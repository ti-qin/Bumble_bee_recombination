#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
BAM_DIR="${PROJECT}/02_mapping/bam"
CANDIDATE_ROOT="${PROJECT}/06_mutations/snp_candidates"
COMBINED_CANDIDATES="${CANDIDATE_ROOT}/all_families.preliminary_denovo_SNP_candidates.tsv.gz"
COMBINED_GENOTYPES="${CANDIDATE_ROOT}/all_families.preliminary_denovo_SNP_genotypes.long.tsv.gz"
OUT="${PROJECT}/06_mutations/IGV_manual_review"

FAMILIES=(C1 C2 C4)
THREADS=4

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

REVIEW_TEMPLATE="${OUT}/summary/denovo_SNP_manual_review_template.tsv"
COUNT_SUMMARY="${OUT}/summary/candidate_count_by_family.tsv"

printf '%s\n' \
    $'candidate_id\tfamily_id\tchromosome\tposition\tmutant_offspring_id\tmutation\treview_decision\treject_reason\tqueen_pattern_ok\tfour_matching_offspring_ok\tmutant_pattern_ok\tmutant_bidirectional_support_ok\tindependent_read_starts_ok\tmapping_quality_ok\tbase_quality_ok\tread_position_ok\tlocal_alignment_ok\tdepth_CNV_ok\trepeat_mappability_ok\treviewer_notes' \
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

    FAMILY_CANDIDATES_GZ="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_SNP_candidates.tsv.gz"
    FAMILY_CANDIDATE_BED="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_SNP_candidates.bed"
    SOURCE_VCF="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_SNP_candidates.vcf.gz"
    REVIEW_VCF="${FAMILY_OUT}/vcf/${family}.preliminary_denovo_SNP_candidates.vcf.gz"

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
        "${FAMILY_OUT}/regions/${family}.preliminary_denovo_SNP_candidates.bed"

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
                    if ($i == "mutant_offspring_id") mutant_col = i
                    if ($i == "mutation") mutation_col = i
                }
                if (!candidate_col || !chromosome_col || !position_col || !mutant_col || !mutation_col) {
                    print "ERROR: required candidate columns not found" > "/dev/stderr"
                    exit 2
                }
                print "candidate_id", "IGV_locus", "mutant_offspring_id", "mutation"
                next
            }
            {
                print $candidate_col, $chromosome_col ":" $position_col, $mutant_col, $mutation_col
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
                    if ($i == "mutant_offspring_id") mutant_col = i
                    if ($i == "mutation") mutation_col = i
                }
                next
            }
            {
                print $candidate_col, $family_col, $chromosome_col, $position_col, \
                      $mutant_col, $mutation_col, "", "", "", "", "", "", \
                      "", "", "", "", "", "", "", ""
            }
        ' >> "${REVIEW_TEMPLATE}"

    # Save GT/DP/GQ/AD/SB exactly as represented in the candidate VCF.
    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT:%DP:%GQ:%AD:%SB]\n' \
        "${REVIEW_VCF}" \
        > "${FAMILY_OUT}/tables/${family}.VCF_GT_DP_GQ_AD_SB.tsv"

    # Link the original whole-genome BAMs and their indexes.
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
        REVIEW_BAM="${FAMILY_OUT}/bam/${sample}.markdup.bam"

        if [[ ! -s "${SOURCE_BAM}" ]]; then
            echo "ERROR: BAM is missing: ${SOURCE_BAM}" >&2
            exit 1
        fi

        samtools quickcheck -v "${SOURCE_BAM}"

        if [[ -s "${SOURCE_BAM}.bai" ]]; then
            SOURCE_BAI="${SOURCE_BAM}.bai"
        elif [[ -s "${SOURCE_BAM%.bam}.bai" ]]; then
            SOURCE_BAI="${SOURCE_BAM%.bam}.bai"
        else
            echo "Creating BAM index: ${SOURCE_BAM}.bai"
            samtools index -@ "${THREADS}" "${SOURCE_BAM}"
            SOURCE_BAI="${SOURCE_BAM}.bai"
        fi

        ln -sfn "${SOURCE_BAM}" "${REVIEW_BAM}"
        ln -sfn "${SOURCE_BAI}" "${REVIEW_BAM}.bai"
    done

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
echo "Whole-genome BAM files were linked, not copied."
echo "Load the reference FASTA, one family candidate VCF, and that family's six BAMs into IGV."
