#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

MANIFEST="${PROJECT}/00_rawdata/metadata/clean_samples.tsv"
BAM_DIR="${PROJECT}/02_mapping/bam"

STATS_DIR="${PROJECT}/02_mapping/stats"
DEPTH_DIR="${PROJECT}/02_mapping/depth"
COV_DIR="${PROJECT}/02_mapping/coverage"

THREADS=128

mkdir -p "${STATS_DIR}" "${DEPTH_DIR}" "${COV_DIR}"

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r sample_id family_id role r1 r2
do
    bam="${BAM_DIR}/${sample_id}.markdup.bam"

    if [[ ! -s "${bam}" ]]; then
        echo "ERROR: BAM not found for ${sample_id}: ${bam}" >&2
        exit 1
    fi

    echo "Processing mapping statistics for ${sample_id}"

    samtools flagstat \
        -@ "${THREADS}" \
        "${bam}" \
        > "${STATS_DIR}/${sample_id}.flagstat.txt"

    samtools stats \
        -@ "${THREADS}" \
        "${bam}" \
        > "${STATS_DIR}/${sample_id}.samtools.stats.txt"

    samtools idxstats \
        -@ "${THREADS}" \
        "${bam}" \
        > "${STATS_DIR}/${sample_id}.idxstats.txt"

    samtools coverage \
        -q 30 \
        -Q 20 \
        "${bam}" \
        > "${COV_DIR}/${sample_id}.coverage.tsv"

    if command -v mosdepth >/dev/null 2>&1; then
        mosdepth \
            -t "${THREADS}" \
            -Q 30 \
            --fast-mode \
            "${DEPTH_DIR}/${sample_id}" \
            "${bam}"
    else
        echo "WARNING: mosdepth not found, skip mosdepth for ${sample_id}"
    fi
done

echo "Mapping statistics completed."
