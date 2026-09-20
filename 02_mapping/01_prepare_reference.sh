#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"

if [[ ! -s "${REF}" ]]; then
    echo "ERROR: reference genome not found: ${REF}" >&2
    exit 1
fi

echo "Reference genome:"
echo "${REF}"

echo "Building bwa-mem2 index..."
bwa-mem2 index "${REF}"

echo "Building samtools fasta index..."
samtools faidx "${REF}"

if command -v gatk >/dev/null 2>&1; then
    echo "Creating GATK sequence dictionary..."
    gatk CreateSequenceDictionary \
        -R "${REF}" \
        -O "${REF%.fa}.dict"
else
    echo "WARNING: gatk not found, skip sequence dictionary."
    echo "This is okay for bwa/samtools/bcftools, but GATK variant calling will need it later."
fi

echo "Reference preparation completed."
