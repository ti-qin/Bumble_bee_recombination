#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

VCF="${PROJECT}/03_variants/filtered/all_samples.diploid_for_ploidy.biallelic.snp.qual30.vcf.gz"
OUT="${PROJECT}/03_variants/genotype_tables/all_samples.gt_dp_ad.tsv.gz"

if [[ ! -s "${VCF}" ]]; then
    echo "ERROR: VCF not found: ${VCF}" >&2
    exit 1
fi

bcftools query \
    -f '[%CHROM\t%POS\t%SAMPLE\t%GT\t%DP\t%AD\n]' \
    "${VCF}" |
gzip -c > "${OUT}"

echo "Saved:"
echo "${OUT}"
