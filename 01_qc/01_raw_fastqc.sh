#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
RAW="${PROJECT}/00_rawdata/fastq"
OUT="${PROJECT}/01_qc/raw_fastqc"

mkdir -p "${OUT}"

find "${RAW}" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name "*.fq" \
    -print0 |
xargs -0 fastqc \
    --threads 12 \
    --outdir "${OUT}"

echo "Raw FastQC completed."
