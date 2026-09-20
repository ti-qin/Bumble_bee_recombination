#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"
MANIFEST="${PROJECT}/00_rawdata/metadata/runs.tsv"

CLEAN="${PROJECT}/01_qc/fastp_per_run/clean"
HTML="${PROJECT}/01_qc/fastp_per_run/html"
JSON="${PROJECT}/01_qc/fastp_per_run/json"
LOG="${PROJECT}/01_qc/fastp_per_run/logs"

THREADS=8

mkdir -p "${CLEAN}" "${HTML}" "${JSON}" "${LOG}"

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r \
    sample_id family_id role run_id r1 r2
do
    sample_clean_dir="${CLEAN}/${sample_id}"
    sample_html_dir="${HTML}/${sample_id}"
    sample_json_dir="${JSON}/${sample_id}"
    sample_log_dir="${LOG}/${sample_id}"

    mkdir -p \
        "${sample_clean_dir}" \
        "${sample_html_dir}" \
        "${sample_json_dir}" \
        "${sample_log_dir}"

    out_r1="${sample_clean_dir}/${run_id}_R1.clean.fastq.gz"
    out_r2="${sample_clean_dir}/${run_id}_R2.clean.fastq.gz"

    html="${sample_html_dir}/${run_id}.fastp.html"
    json="${sample_json_dir}/${run_id}.fastp.json"
    log="${sample_log_dir}/${run_id}.fastp.log"

    echo "Processing ${sample_id}: ${run_id}"

    fastp \
        --in1 "${r1}" \
        --in2 "${r2}" \
        --out1 "${out_r1}" \
        --out2 "${out_r2}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --n_base_limit 5 \
        --length_required 50 \
        --thread "${THREADS}" \
        --html "${html}" \
        --json "${json}" \
        > "${log}" 2>&1

    if [[ ! -s "${out_r1}" || ! -s "${out_r2}" ]]; then
        echo "ERROR: empty output for ${sample_id}: ${run_id}" >&2
        exit 1
    fi
done

echo "All fastp runs completed."
