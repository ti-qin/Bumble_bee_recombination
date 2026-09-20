#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
MANIFEST="${PROJECT}/00_rawdata/metadata/clean_samples.tsv"

BAM_DIR="${PROJECT}/02_mapping/bam"
TMP_DIR="${PROJECT}/02_mapping/tmp"
LOG_DIR="${PROJECT}/logs/mapping"

THREADS=128

mkdir -p "${BAM_DIR}" "${TMP_DIR}" "${LOG_DIR}"

if [[ ! -s "${REF}" ]]; then
    echo "ERROR: reference genome not found: ${REF}" >&2
    exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: clean sample manifest not found: ${MANIFEST}" >&2
    exit 1
fi

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r sample_id family_id role r1 r2
do
    echo "=========================================="
    echo "Mapping sample: ${sample_id}"
    echo "Family: ${family_id}"
    echo "Role: ${role}"
    echo "R1: ${r1}"
    echo "R2: ${r2}"
    echo "=========================================="

    final_bam="${BAM_DIR}/${sample_id}.markdup.bam"

    if [[ -s "${final_bam}" && -s "${final_bam}.bai" ]]; then
        echo "Skip ${sample_id}: BAM already exists."
        continue
    fi

    sample_tmp="${TMP_DIR}/${sample_id}"
    mkdir -p "${sample_tmp}"

    name_bam="${sample_tmp}/${sample_id}.name.bam"
    fixmate_bam="${sample_tmp}/${sample_id}.fixmate.bam"
    coord_bam="${sample_tmp}/${sample_id}.coord.bam"

    rg="@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA\tLB:${sample_id}\tPU:${sample_id}"

    log_file="${LOG_DIR}/${sample_id}.mapping.log"

    {
        echo "[1/5] bwa-mem2 mem + name sort"

        bwa-mem2 mem \
            -t "${THREADS}" \
            -R "${rg}" \
            "${REF}" \
            "${r1}" \
            "${r2}" |
        samtools sort \
            -n \
            -@ "${THREADS}" \
            -o "${name_bam}" \
            -

        echo "[2/5] samtools fixmate"

        samtools fixmate \
            -m \
            -@ "${THREADS}" \
            "${name_bam}" \
            "${fixmate_bam}"

        echo "[3/5] coordinate sort"

        samtools sort \
            -@ "${THREADS}" \
            -o "${coord_bam}" \
            "${fixmate_bam}"

        echo "[4/5] mark duplicates"

        samtools markdup \
            -@ "${THREADS}" \
            "${coord_bam}" \
            "${final_bam}"

        echo "[5/5] index BAM"

        samtools index \
            -@ "${THREADS}" \
            "${final_bam}"

        echo "Cleaning temporary files"

        rm -f \
            "${name_bam}" \
            "${fixmate_bam}" \
            "${coord_bam}"

        rmdir "${sample_tmp}" || true

        echo "Completed: ${sample_id}"
    } > "${log_file}" 2>&1

done

echo "All mapping jobs completed."
