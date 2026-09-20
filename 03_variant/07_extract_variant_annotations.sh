#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

VCF="${PROJECT}/03_variants/gatk_final/filtered/C1_C2_C4.biallelic.snps.raw.vcf.gz"

OUT_DIR="${PROJECT}/03_variants/gatk_final/qc/annotations"
OUT="${OUT_DIR}/biallelic_snp_annotations.tsv.gz"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${VCF}" ]]; then
    echo "ERROR: VCF not found: ${VCF}" >&2
    exit 1
fi

{
    printf "chromosome\tposition\tQUAL\tQD\tMQ\tFS\tSOR\tMQRankSum\tReadPosRankSum\n"

    bcftools query \
        -f '%CHROM\t%POS\t%QUAL\t%INFO/QD\t%INFO/MQ\t%INFO/FS\t%INFO/SOR\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n' \
        "${VCF}"
} | gzip -c > "${OUT}"

echo "Saved:"
echo "${OUT}"

echo
echo "Number of SNPs:"
zcat "${OUT}" | tail -n +2 | wc -l
