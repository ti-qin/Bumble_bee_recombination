#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
CLEAN="${PROJECT}/01_qc/clean_reads"
OUT="${PROJECT}/00_rawdata/metadata/clean_samples.tsv"

mkdir -p "$(dirname "${OUT}")"

printf "sample_id\tfamily_id\trole\tR1\tR2\n" > "${OUT}"

for r1 in "${CLEAN}"/*_R1.clean.fastq.gz; do
    [[ -s "${r1}" ]] || continue

    sample_id=$(basename "${r1}" _R1.clean.fastq.gz)
    r2="${CLEAN}/${sample_id}_R2.clean.fastq.gz"

    if [[ ! -s "${r2}" ]]; then
        echo "ERROR: missing R2 for ${sample_id}" >&2
        exit 1
    fi

    family_id=${sample_id%%_*}

    if [[ "${sample_id}" == *_Q ]]; then
        role="queen"
    else
        role="offspring"
    fi

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "${sample_id}" \
        "${family_id}" \
        "${role}" \
        "${r1}" \
        "${r2}" \
        >> "${OUT}"
done

echo "Created: ${OUT}"
echo
column -t -s $'\t' "${OUT}"
