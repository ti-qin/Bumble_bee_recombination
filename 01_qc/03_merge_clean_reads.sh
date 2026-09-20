#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
MANIFEST="${PROJECT}/00_rawdata/metadata/runs.tsv"

RUN_CLEAN="${PROJECT}/01_qc/fastp_per_run/clean"
FINAL="${PROJECT}/01_qc/clean_reads"

THREADS=32

mkdir -p "${FINAL}"

if ! command -v pigz >/dev/null 2>&1; then
    echo "ERROR: pigz was not found." >&2
    echo "Install with: conda install -c conda-forge pigz" >&2
    exit 1
fi

mapfile -t samples < <(
    awk -F '\t' 'NR > 1 {print $1}' "${MANIFEST}" |
    sort -u
)

for sample_id in "${samples[@]}"; do

    mapfile -t run_ids < <(
        awk -F '\t' -v sample="${sample_id}" \
            'NR > 1 && $1 == sample {print $4}' \
            "${MANIFEST}" |
        sort
    )

    r1_files=()
    r2_files=()

    for run_id in "${run_ids[@]}"; do
        r1="${RUN_CLEAN}/${sample_id}/${run_id}_R1.clean.fastq.gz"
        r2="${RUN_CLEAN}/${sample_id}/${run_id}_R2.clean.fastq.gz"

        if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
            echo "ERROR: missing clean FASTQ for ${sample_id}, ${run_id}" >&2
            exit 1
        fi

        r1_files+=("${r1}")
        r2_files+=("${r2}")
    done

    final_r1="${FINAL}/${sample_id}_R1.clean.fastq.gz"
    final_r2="${FINAL}/${sample_id}_R2.clean.fastq.gz"

    echo "${sample_id}: ${#run_ids[@]} sequencing run(s)"

    if [[ ${#run_ids[@]} -eq 1 ]]; then
        ln -sfn "$(realpath "${r1_files[0]}")" "${final_r1}"
        ln -sfn "$(realpath "${r2_files[0]}")" "${final_r2}"
    else
        tmp_r1="${final_r1}.tmp"
        tmp_r2="${final_r2}.tmp"

        rm -f "${tmp_r1}" "${tmp_r2}"

        pigz -dc "${r1_files[@]}" |
            pigz -p "${THREADS}" > "${tmp_r1}"

        pigz -dc "${r2_files[@]}" |
            pigz -p "${THREADS}" > "${tmp_r2}"

        mv "${tmp_r1}" "${final_r1}"
        mv "${tmp_r2}" "${final_r2}"
    fi
done

echo "Final clean reads:"
ls -lh "${FINAL}"
