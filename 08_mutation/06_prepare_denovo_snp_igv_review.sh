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
FLANK_BP=5000
THREADS=4

# full: link the original whole-genome BAM and index into the
#       review package (recommended when IGV can access the server).
# subset: create smaller BAMs containing only candidate +/- FLANK_BP.
BAM_MODE="${BAM_MODE:-full}"

if [[ "${BAM_MODE}" != "full" && "${BAM_MODE}" != "subset" ]]; then
    echo "ERROR: BAM_MODE must be 'full' or 'subset': ${BAM_MODE}" >&2
    exit 1
fi

for program in python samtools bcftools; do
    command -v "${program}" >/dev/null 2>&1 || {
        echo "ERROR: ${program} is not available in PATH" >&2
        exit 1
    }
done

for input in "${REF}" "${REF}.fai" "${COMBINED_CANDIDATES}" "${COMBINED_GENOTYPES}"; do
    [[ -s "${input}" ]] || {
        echo "ERROR: required input is missing or empty: ${input}" >&2
        exit 1
    }
done

mkdir -p "${OUT}"/{reference,summary}
ln -sfn "${REF}" "${OUT}/reference/Bombus_vestalis.fa"
ln -sfn "${REF}.fai" "${OUT}/reference/Bombus_vestalis.fa.fai"

# Build family-specific review tables, merged +/-5 kb BED files,
# IGV locus lists, and a blank manual-review template.
PROJECT="${PROJECT}" OUT="${OUT}" FLANK_BP="${FLANK_BP}" python - <<'PY'
from pathlib import Path
import os
import pandas as pd

project = Path(os.environ["PROJECT"])
out = Path(os.environ["OUT"])
flank = int(os.environ["FLANK_BP"])
candidate_root = project / "06_mutations/snp_candidates"

candidates = pd.read_csv(
    candidate_root / "all_families.preliminary_denovo_SNP_candidates.tsv.gz",
    sep="\t",
    dtype={"chromosome": str},
)
genotypes = pd.read_csv(
    candidate_root / "all_families.preliminary_denovo_SNP_genotypes.long.tsv.gz",
    sep="\t",
    dtype={"chromosome": str},
)

sizes = {}
with (project / "00_rawdata/reference/Bombus_vestalis.fa.fai").open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        sizes[fields[0]] = int(fields[1])

review_columns = [
    "candidate_id", "family_id", "chromosome", "position",
    "mutant_offspring_id", "mutation", "review_decision",
    "reject_reason", "queen_pattern_ok", "four_matching_offspring_ok",
    "mutant_pattern_ok", "mutant_bidirectional_support_ok",
    "independent_read_starts_ok", "mapping_quality_ok",
    "base_quality_ok", "read_position_ok", "local_alignment_ok",
    "depth_CNV_ok", "repeat_mappability_ok", "reviewer_notes",
]
review_rows = []

for family in ["C1", "C2", "C4"]:
    family_dir = out / family
    for subdir in ["bam", "vcf", "tables", "regions", "snapshots"]:
        (family_dir / subdir).mkdir(parents=True, exist_ok=True)

    fc = candidates.loc[candidates["family_id"].astype(str).eq(family)].copy()
    fg = genotypes.loc[genotypes["family_id"].astype(str).eq(family)].copy()

    fc.to_csv(
        family_dir / "tables" / f"{family}.preliminary_candidates.tsv",
        sep="\t", index=False,
    )
    fg.to_csv(
        family_dir / "tables" / f"{family}.candidate_genotypes.long.tsv",
        sep="\t", index=False,
    )

    intervals = []
    for row in fc.itertuples(index=False):
        chrom = str(row.chromosome)
        if chrom not in sizes:
            raise ValueError(f"{family}: chromosome not in FASTA index: {chrom}")
        position = int(row.position)
        start = max(0, position - 1 - flank)
        end = min(sizes[chrom], position + flank)
        intervals.append((chrom, start, end))

    merged = []
    for chrom, start, end in sorted(intervals):
        if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
            merged.append([chrom, start, end])
        else:
            merged[-1][2] = max(merged[-1][2], end)

    bed = family_dir / "regions" / f"{family}.candidate_regions.{flank}bp_flank.bed"
    with bed.open("w") as handle:
        for chrom, start, end in merged:
            handle.write(f"{chrom}\t{start}\t{end}\n")

    loci = family_dir / "regions" / f"{family}.candidate_loci.tsv"
    fc[["candidate_id", "chromosome", "position", "mutant_offspring_id", "mutation"]].to_csv(
        loci, sep="\t", index=False,
    )

    for row in fc.itertuples(index=False):
        review_rows.append({
            "candidate_id": row.candidate_id,
            "family_id": row.family_id,
            "chromosome": row.chromosome,
            "position": row.position,
            "mutant_offspring_id": row.mutant_offspring_id,
            "mutation": row.mutation,
            **{column: "" for column in review_columns[6:]},
        })

pd.DataFrame(review_rows, columns=review_columns).to_csv(
    out / "summary/denovo_SNP_manual_review_template.tsv",
    sep="\t", index=False,
)
print(f"Prepared review tables for {len(candidates)} candidates")
PY

for family in "${FAMILIES[@]}"; do
    family_out="${OUT}/${family}"
    region_bed="${family_out}/regions/${family}.candidate_regions.${FLANK_BP}bp_flank.bed"
    source_vcf="${CANDIDATE_ROOT}/${family}/${family}.preliminary_denovo_SNP_candidates.vcf.gz"
    review_vcf="${family_out}/vcf/${family}.preliminary_denovo_SNP_candidates.vcf.gz"

    [[ -s "${source_vcf}" ]] || {
        echo "ERROR: family candidate VCF missing: ${source_vcf}" >&2
        exit 1
    }

    bcftools view \
        --output-type z \
        --output-file "${review_vcf}" \
        "${source_vcf}"
    bcftools index --force --tbi "${review_vcf}"

    samples=("${family}_Q" "${family}_1" "${family}_2" "${family}_3" "${family}_4" "${family}_5")
    for sample in "${samples[@]}"; do
        source_bam="${BAM_DIR}/${sample}.markdup.bam"
        [[ -s "${source_bam}" ]] || {
            echo "ERROR: BAM missing: ${source_bam}" >&2
            exit 1
        }

        samtools quickcheck -v "${source_bam}"

        if [[ "${BAM_MODE}" == "full" ]]; then
            review_bam="${family_out}/bam/${sample}.markdup.bam"

            if [[ -s "${source_bam}.bai" ]]; then
                source_bai="${source_bam}.bai"
            elif [[ -s "${source_bam%.bam}.bai" ]]; then
                source_bai="${source_bam%.bam}.bai"
            else
                echo "Creating missing BAM index: ${source_bam}.bai"
                samtools index -@ "${THREADS}" "${source_bam}"
                source_bai="${source_bam}.bai"
            fi

            ln -sfn "${source_bam}" "${review_bam}"
            ln -sfn "${source_bai}" "${review_bam}.bai"
        else
            review_bam="${family_out}/bam/${sample}.candidate_regions.bam"
            samtools view \
                -@ "${THREADS}" \
                -bh \
                -L "${region_bed}" \
                -o "${review_bam}" \
                "${source_bam}"
            samtools index -@ "${THREADS}" "${review_bam}"
            samtools quickcheck -v "${review_bam}"
        fi
    done

    bcftools query \
        -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE:%GT:%DP:%GQ:%AD:%SB]\n' \
        "${review_vcf}" \
        > "${family_out}/tables/${family}.VCF_GT_DP_GQ_AD_SB.tsv"
done

echo
echo "IGV review package created:"
echo "${OUT}"
echo "BAM mode: ${BAM_MODE}"
echo
echo "Load the reference FASTA, one family VCF, and the six BAMs from that family into IGV."
