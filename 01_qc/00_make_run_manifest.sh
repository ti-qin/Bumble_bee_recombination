#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
RAW="${PROJECT}/00_rawdata/fastq"
OUT="${PROJECT}/00_rawdata/metadata/runs.tsv"

mkdir -p "$(dirname "${OUT}")"

printf "sample_id\tfamily_id\trole\trun_id\tR1\tR2\n" > "${OUT}"

for sample_dir in "${RAW}"/C*; do
    [[ -d "${sample_dir}" ]] || continue

    sample_id=$(basename "${sample_dir}")
    family_id=${sample_id%%_*}

    if [[ "${sample_id}" == *_Q ]]; then
        role="queen"
    else
        role="offspring"
    fi

    mapfile -t r1_files < <(
        find "${sample_dir}" \
            -maxdepth 1 \
            -type f \
            -name "*_1.fq" \
            | sort
    )

    if [[ ${#r1_files[@]} -eq 0 ]]; then
        echo "ERROR: no R1 file found in ${sample_dir}" >&2
        exit 1
    fi

    for r1 in "${r1_files[@]}"; do
        r2="${r1%_1.fq}_2.fq"

        if [[ ! -s "${r2}" ]]; then
            echo "ERROR: missing R2 for ${r1}" >&2
            exit 1
        fi

        run_id=$(basename "${r1%_1.fq}")

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${sample_id}" \
            "${family_id}" \
            "${role}" \
            "${run_id}" \
            "${r1}" \
            "${r2}" \
            >> "${OUT}"
    done
done

echo "Created: ${OUT}"
echo
column -t -s $'\t' "${OUT}"
