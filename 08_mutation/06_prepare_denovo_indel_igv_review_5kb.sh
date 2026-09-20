#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
BAM_DIR="${PROJECT}/02_mapping/bam"
CANDIDATE_ROOT="${PROJECT}/06_mutations/indel_candidates"
COMBINED_CANDIDATES="${CANDIDATE_ROOT}/all_families.preliminary_denovo_INDEL_candidates.tsv.gz"
COMBINED_GENOTYPES="${CANDIDATE_ROOT}/all_families.preliminary_denovo_INDEL_genotypes.long.tsv.gz"
OUT="${PROJECT}/06_mutations/IGV_indel_manual_review"

FAMILIES=(C1 C2 C4)
THREADS=4
FLANK_BP=5000
IGV_VIEW_FLANK_BP=150

for program in samtools bcftools awk gzip; do
    if ! command -v "${program}" >/dev/null 2>&1; then
        echo "ERROR: required program is not available: ${program}" >&2
        exit 1
    fi
done

for input_file in \
    "${REF}" \
    "${REF}.fai" \
    "${COMBINED_CANDIDATES}" \
    "${COMBINED_GENOTYPES}"
do
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: required input is missing or empty: ${input_file}" >&2
        exit 1
    fi
done

mkdir -p "${OUT}/reference" "${OUT}/summary"

# These are links to the exact reference used for mapping/calling.
ln -sfn "${REF}" "${OUT}/reference/Bombus_vestalis.fa"
ln -sfn "${REF}.fai" "${OUT}/reference/Bombus_vestalis.fa.fai"

REVIEW_TEMPLATE="${OUT}/summary/denovo_INDEL_manual_review_template.tsv"
COUNT_SUMMARY="${OUT}/summary/candidate_count_by_family.tsv"

printf '%s\n' \
    $'candidate_id\tfamily_id\tchromosome\tposition\tmutation_type\tindel_length\tmutant_offspring_id\tmutation\treview_decision\treject_reason\tqueen_pattern_ok\tfour_matching_offspring_ok\tmutant_pattern_ok\texact_indel_sequence_ok\tmutant_bidirectional_support_ok\tindependent_read_starts_ok\tmapping_quality_ok\tbase_quality_ok\tread_position_ok\tlocal_alignment_ok\tnearby_soft_clipping_ok\tdepth_CNV_ok\thomopolymer_repeat_ok\trepeat_mappability_ok\treviewer_notes' \
    > "${REVIEW_TEMPLATE}"

printf 'family_id\tn_candidates\n' > "${COUNT_SUMMARY}"

for family in "${FAMILIES[@]}"; do
    echo "Preparing IGV review files for ${family}"

    FAMILY_OUT="${OUT}/${family}"
    mkdir -p \
        "${FAMILY_OUT}/bam" \
        "${FAMILY_OUT}/vcf" \
        "${FAMILY_OUT}/tables" \
        "${FAMILY_OUT}/regions" \
        "${FAMILY_OUT}/snapshots"

    FAMILY_CANDIDATES_GZ="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_INDEL_candidates.tsv.gz"
    FAMILY_CANDIDATE_BED="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_INDEL_candidates.bed"
    SOURCE_VCF="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_INDEL_candidates.vcf.gz"
    REVIEW_VCF="${FAMILY_OUT}/vcf/${family}.preliminary_denovo_INDEL_candidates.vcf.gz"

    for family_input in \
        "${FAMILY_CANDIDATES_GZ}" \
        "${FAMILY_CANDIDATE_BED}" \
        "${SOURCE_VCF}"
    do
        if [[ ! -e "${family_input}" ]]; then
            echo "ERROR: family input is missing: ${family_input}" >&2
            exit 1
        fi
    done

    # Human-readable candidate table.
    gzip -dc "${FAMILY_CANDIDATES_GZ}" \
        > "${FAMILY_OUT}/tables/${family}.preliminary_candidates.tsv"

    cp -f \
        "${FAMILY_CANDIDATE_BED}" \
        "${FAMILY_OUT}/regions/${family}.preliminary_denovo_INDEL_candidates.bed"

    N_CANDIDATES=$(
        gzip -dc "${FAMILY_CANDIDATES_GZ}" |
            awk 'END { print (NR > 0 ? NR - 1 : 0) }'
    )

    if [[ "${N_CANDIDATES}" -eq 0 ]]; then
        echo "${family}: no preliminary INDEL candidates; skipping BAM extraction."

        bcftools view \
            --output-type z \
            --output-file "${REVIEW_VCF}" \
            "${SOURCE_VCF}"

        bcftools index --force --tbi "${REVIEW_VCF}"

        printf "candidate_id\tIGV_locus\tmutation_type\tindel_length\tmutant_offspring_id\tmutation\n" \
            > "${FAMILY_OUT}/regions/${family}.candidate_loci.tsv"

        printf "chromosome\tposition\tREF\tALT\tfamily_samples_GT_DP_GQ_AD_SB\n" \
            > "${FAMILY_OUT}/tables/${family}.VCF_GT_DP_GQ_AD_SB.tsv"

        printf "%s\t%s\n" "${family}" "${N_CANDIDATES}" >> "${COUNT_SUMMARY}"
        continue
    fi

    # Expand each 1-bp candidate interval by +/- FLANK_BP while
    # respecting chromosome boundaries from the FASTA index.
    EXPANDED_BED="${FAMILY_OUT}/regions/${family}.candidate_regions.${FLANK_BP}bp_flank.unmerged.bed"
    REVIEW_REGION_BED="${FAMILY_OUT}/regions/${family}.candidate_regions.${FLANK_BP}bp_flank.merged.bed"

    awk -v OFS='\t' -v flank="${FLANK_BP}" '
        NR == FNR {
            chromosome_length[$1] = $2
            next
        }
        {
            chromosome = $1
            if (!(chromosome in chromosome_length)) {
                print "ERROR: chromosome not found in FASTA index: " chromosome > "/dev/stderr"
                exit 2
            }
            start = $2 - flank
            end = $3 + flank
            if (start < 0) start = 0
            if (end > chromosome_length[chromosome]) end = chromosome_length[chromosome]
            print chromosome, start, end
        }
    ' "${REF}.fai" "${FAMILY_CANDIDATE_BED}" |
        LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
        > "${EXPANDED_BED}"

    # Merge overlapping regions to prevent duplicate alignments.
    awk -v OFS='\t' '
        NR == 1 {
            chromosome = $1
            start = $2
            end = $3
            next
        }
        $1 == chromosome && $2 <= end {
            if ($3 > end) end = $3
            next
        }
        {
            print chromosome, start, end
            chromosome = $1
            start = $2
            end = $3
        }
        END {
            if (NR > 0) print chromosome, start, end
        }
    ' "${EXPANDED_BED}" > "${REVIEW_REGION_BED}"

    if [[ ! -s "${REVIEW_REGION_BED}" ]]; then
        echo "ERROR: no candidate review regions were generated for ${family}" >&2
        exit 1
    fi

    # Keep a small candidate-only VCF containing all six family samples.
    bcftools view \
        --output-type z \
        --output-file "${REVIEW_VCF}" \
        "${SOURCE_VCF}"

    bcftools index \
        --force \
        --tbi \
        "${REVIEW_VCF}"

    # Extract the six-sample quantitative evidence for this family.
    gzip -dc "${COMBINED_GENOTYPES}" |
        awk -F '\t' -v OFS='\t' -v family="${family}" '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "family_id") family_col = i
                }
                if (!family_col) {
                    print "ERROR: family_id column not found" > "/dev/stderr"
                    exit 2
                }
                print
                next
            }
            $family_col == family { print }
        ' > "${FAMILY_OUT}/tables/${family}.candidate_genotypes.long.tsv"

    # Generate a concise list of loci for direct copy/paste into IGV.
    gzip -dc "${FAMILY_CANDIDATES_GZ}" |
        awk -F '\t' -v OFS='\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "candidate_id") candidate_col = i
                    if ($i == "chromosome") chromosome_col = i
                    if ($i == "position") position_col = i
                    if ($i == "mutation_type") type_col = i
                    if ($i == "indel_length") length_col = i
                    if ($i == "mutant_offspring_id") mutant_col = i
                    if ($i == "mutation") mutation_col = i
                }
                if (!candidate_col || !chromosome_col || !position_col || !type_col || !length_col || !mutant_col || !mutation_col) {
                    print "ERROR: required candidate columns not found" > "/dev/stderr"
                    exit 2
                }
                print "candidate_id", "IGV_locus", "mutation_type", "indel_length", "mutant_offspring_id", "mutation"
                next
            }
            {
                print $candidate_col, $chromosome_col ":" $position_col, $type_col, $length_col, $mutant_col, $mutation_col
            }
        ' > "${FAMILY_OUT}/regions/${family}.candidate_loci.tsv"

    # Add blank manual-review rows.
    gzip -dc "${FAMILY_CANDIDATES_GZ}" |
        awk -F '\t' -v OFS='\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "candidate_id") candidate_col = i
                    if ($i == "family_id") family_col = i
                    if ($i == "chromosome") chromosome_col = i
                    if ($i == "position") position_col = i
                    if ($i == "mutation_type") type_col = i
                    if ($i == "indel_length") length_col = i
                    if ($i == "mutant_offspring_id") mutant_col = i
                    if ($i == "mutation") mutation_col = i
                }
                next
            }
            {
                printf "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s", \
                    $candidate_col, OFS, $family_col, OFS, \
                    $chromosome_col, OFS, $position_col, OFS, \
                    $type_col, OFS, $length_col, OFS, \
                    $mutant_col, OFS, $mutation_col
                for (j = 1; j <= 17; j++) {
                    printf "%s", OFS
                }
                printf "\n"
            }
        ' >> "${REVIEW_TEMPLATE}"

    # Save GT/DP/GQ/AD/SB exactly as represented in the candidate VCF.
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT:%DP:%GQ:%AD:%SB]\n' \
        "${REVIEW_VCF}" \
        > "${FAMILY_OUT}/tables/${family}.VCF_GT_DP_GQ_AD_SB.tsv"

    # Extract candidate +/- FLANK_BP regions from the six whole-genome BAMs.
    samples=(
        "${family}_Q"
        "${family}_1"
        "${family}_2"
        "${family}_3"
        "${family}_4"
        "${family}_5"
    )

    for sample in "${samples[@]}"; do
        SOURCE_BAM="${BAM_DIR}/${sample}.markdup.bam"
        REVIEW_BAM="${FAMILY_OUT}/bam/${sample}.candidate_regions.${FLANK_BP}bp_flank.bam"

        if [[ ! -s "${SOURCE_BAM}" ]]; then
            echo "ERROR: BAM is missing: ${SOURCE_BAM}" >&2
            exit 1
        fi

        samtools quickcheck -v "${SOURCE_BAM}"

        if [[ ! -s "${SOURCE_BAM}.bai" && ! -s "${SOURCE_BAM%.bam}.bai" ]]; then
            echo "Creating BAM index: ${SOURCE_BAM}.bai"
            samtools index -@ "${THREADS}" "${SOURCE_BAM}"
        fi

        samtools view \
            -@ "${THREADS}" \
            -bh \
            -L "${REVIEW_REGION_BED}" \
            -o "${REVIEW_BAM}" \
            "${SOURCE_BAM}"

        samtools index \
            -@ "${THREADS}" \
            "${REVIEW_BAM}"

        samtools quickcheck -v "${REVIEW_BAM}"
    done

    # Optional IGV batch file. Run IGV from the review-package root
    # so these relative paths resolve after copying the package locally.
    IGV_BATCH="${FAMILY_OUT}/regions/${family}.igv_batch.txt"
    {
        echo "new"
        echo "genome reference/Bombus_vestalis.fa"
        echo "load ${family}/vcf/${family}.preliminary_denovo_INDEL_candidates.vcf.gz"
        for sample in "${samples[@]}"; do
            echo "load ${family}/bam/${sample}.candidate_regions.${FLANK_BP}bp_flank.bam"
        done
        echo "snapshotDirectory ${family}/snapshots"
        echo "maxPanelHeight 1000"

        gzip -dc "${FAMILY_CANDIDATES_GZ}" |
            awk -F '\t' -v flank="${IGV_VIEW_FLANK_BP}" '
                NR == 1 {
                    for (i = 1; i <= NF; i++) {
                        if ($i == "candidate_id") id_col = i
                        if ($i == "chromosome") chr_col = i
                        if ($i == "position") pos_col = i
                    }
                    next
                }
                {
                    start = $pos_col - flank
                    if (start < 1) start = 1
                    end = $pos_col + flank
                    filename = $id_col
                    gsub(/[^A-Za-z0-9_.-]/ , "_", filename)
                    print "goto " $chr_col ":" start "-" end
                    print "snapshot " filename ".png"
                }
            '
        echo "exit"
    } > "${IGV_BATCH}"

    N_CANDIDATES=$(
        gzip -dc "${FAMILY_CANDIDATES_GZ}" |
            awk 'END { print (NR > 0 ? NR - 1 : 0) }'
    )
    printf '%s\t%s\n' "${family}" "${N_CANDIDATES}" >> "${COUNT_SUMMARY}"

    echo "${family}: ${N_CANDIDATES} candidates"
done

echo
echo "IGV manual-review package created:"
echo "${OUT}"
echo
echo "Local BAM files contain candidate positions +/- ${FLANK_BP} bp."
echo "Load the reference FASTA, one family candidate VCF, and that family's six BAMs into IGV."
