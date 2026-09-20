#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
CHR_BED="${PROJECT}/00_rawdata/reference/primary_chromosomes.bed"
BAM_LIST="${PROJECT}/03_variants/joint_calling/bam.list"

OUT_DIR="${PROJECT}/03_variants/joint_calling"
LOG_DIR="${PROJECT}/logs/variants"

THREADS=128

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

RAW_BCF="${OUT_DIR}/all_samples.diploid_for_ploidy.raw.bcf"
RAW_VCF="${OUT_DIR}/all_samples.diploid_for_ploidy.raw.vcf.gz"

if [[ ! -s "${REF}" ]]; then
    echo "ERROR: reference not found: ${REF}" >&2
    exit 1
fi

if [[ ! -s "${CHR_BED}" ]]; then
    echo "ERROR: primary chromosome BED not found: ${CHR_BED}" >&2
    exit 1
fi

if [[ ! -s "${BAM_LIST}" ]]; then
    echo "ERROR: BAM list not found: ${BAM_LIST}" >&2
    exit 1
fi

echo "Running bcftools mpileup/call for ploidy check..."

bcftools mpileup \
    --threads "${THREADS}" \
    -f "${REF}" \
    -R "${CHR_BED}" \
    -q 30 \
    -Q 20 \
    --annotate FORMAT/AD,FORMAT/DP \
    --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -b "${BAM_LIST}" \
    -Ou |
bcftools call \
    --threads "${THREADS}" \
    -m \
    -v \
    -Ob \
    -o "${RAW_BCF}"

bcftools view \
    --threads "${THREADS}" \
    -Oz \
    -o "${RAW_VCF}" \
    "${RAW_BCF}"

bcftools index \
    --threads "${THREADS}" \
    -t "${RAW_VCF}"

echo "Saved:"
echo "${RAW_VCF}"
