#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

IN="${PROJECT}/03_variants/joint_calling/all_samples.diploid_for_ploidy.raw.vcf.gz"
OUT="${PROJECT}/03_variants/filtered/all_samples.diploid_for_ploidy.biallelic.snp.qual30.vcf.gz"

THREADS=128

if [[ ! -s "${IN}" ]]; then
    echo "ERROR: input VCF not found: ${IN}" >&2
    exit 1
fi

bcftools view \
    --threads "${THREADS}" \
    -m2 \
    -M2 \
    -v snps \
    "${IN}" |
bcftools filter \
    --threads "${THREADS}" \
    -i 'QUAL>=30' \
    -Oz \
    -o "${OUT}"

bcftools index \
    --threads "${THREADS}" \
    -t "${OUT}"

echo "Filtered VCF saved:"
echo "${OUT}"

bcftools stats "${OUT}" \
    > "${OUT%.vcf.gz}.stats.txt"
