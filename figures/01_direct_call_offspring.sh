#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# 1. Paths
# ============================================================

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
INTERVALS="${PROJECT}/00_rawdata/reference/primary_chromosomes.bed"

BAM_DIR="${PROJECT}/02_mapping/bam"

OUT_ROOT="${PROJECT}/figure/offspring_NJ"
CALL_DIR="${OUT_ROOT}/01_direct_call"
TMP_DIR="${OUT_ROOT}/tmp"

LOG_DIR="${PROJECT}/logs/figure1_offspring_NJ"

mkdir -p \
    "${CALL_DIR}" \
    "${TMP_DIR}" \
    "${LOG_DIR}"


# ============================================================
# 2. Samples
# ============================================================

SAMPLES=(
    C1_1 C1_2 C1_3 C1_4 C1_5
    C2_1 C2_2 C2_3 C2_4 C2_5
    C4_1 C4_2 C4_3 C4_4 C4_5
)


# ============================================================
# 3. Parameters
# ============================================================

# Fifteen BAMs are processed in the same HaplotypeCaller run.
# Increase this if the server has sufficient memory.
JAVA_MEM="64g"

PAIRHMM_THREADS=96

STAND_CALL_CONF=30.0
MIN_BASE_QUALITY=15

WORK_VCF="${CALL_DIR}/offspring15.direct_diploid_SB.raw.vcf.gz"

GATK_LOG="${LOG_DIR}/offspring15.direct_diploid_SB.HaplotypeCaller.log"


# ============================================================
# 4. Check programs and inputs
# ============================================================

for program in gatk samtools bcftools
do
    if ! command -v "${program}" >/dev/null 2>&1
    then
        echo "ERROR: ${program} was not found in PATH." >&2
        exit 1
    fi
done

if [[ ! -f "${REF}" ]]
then
    echo "ERROR: reference not found: ${REF}" >&2
    exit 1
fi

if [[ ! -f "${INTERVALS}" ]]
then
    echo "ERROR: interval file not found: ${INTERVALS}" >&2
    exit 1
fi

if [[ ! -d "${BAM_DIR}" ]]
then
    echo "ERROR: BAM directory not found: ${BAM_DIR}" >&2
    exit 1
fi


# ============================================================
# 5. Locate BAM file
# ============================================================

find_sample_bam()
{
    local sample="$1"
    local bam

    # Modify these candidates if your BAM naming differs.
    local candidates=(
        "${BAM_DIR}/${sample}.markdup.bam"
        "${BAM_DIR}/${sample}.sorted.markdup.bam"
        "${BAM_DIR}/${sample}.markdup.sorted.bam"
        "${BAM_DIR}/${sample}.dedup.bam"
        "${BAM_DIR}/${sample}.bam"
    )

    for bam in "${candidates[@]}"
    do
        if [[ -f "${bam}" ]]
        then
            printf '%s\n' "${bam}"
            return 0
        fi
    done

    local -a hits=()

    mapfile -t hits < <(
        find "${BAM_DIR}" \
            -maxdepth 1 \
            -type f \
            -name "${sample}*markdup*.bam" \
            | sort
    )

    if (( ${#hits[@]} == 1 ))
    then
        printf '%s\n' "${hits[0]}"
        return 0
    fi

    if (( ${#hits[@]} > 1 ))
    then
        echo "ERROR: multiple BAM files found for ${sample}:" >&2
        printf '  %s\n' "${hits[@]}" >&2
        return 1
    fi

    echo "ERROR: no BAM found for ${sample}" >&2
    return 1
}


# ============================================================
# 6. Build BAM argument array
# ============================================================

bam_args=()

MANIFEST="${CALL_DIR}/offspring15.sample_bam.tsv"

printf 'sample_id\tbam\tbam_SM\n' > "${MANIFEST}"

for sample in "${SAMPLES[@]}"
do
    bam="$(find_sample_bam "${sample}")"

    echo "Checking ${sample}: ${bam}"

    # Basic BAM integrity check
    if ! samtools quickcheck "${bam}"
    then
        echo "ERROR: BAM failed samtools quickcheck: ${bam}" >&2
        exit 1
    fi

    # Ensure BAM index exists
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]
    then
        echo "Creating BAM index: ${bam}"
        samtools index \
            -@ "${PAIRHMM_THREADS}" \
            "${bam}"
    fi

    # Extract unique SM values from read groups
    mapfile -t bam_sm_values < <(
        samtools view -H "${bam}" |
        awk -F '\t' '
            $1 == "@RG" {
                for (i = 1; i <= NF; i++) {
                    if ($i ~ /^SM:/) {
                        sub(/^SM:/, "", $i)
                        print $i
                    }
                }
            }
        ' |
        sort -u
    )

    if (( ${#bam_sm_values[@]} == 0 ))
    then
        echo "ERROR: no SM tag found in BAM read groups:" >&2
        echo "${bam}" >&2
        exit 1
    fi

    if (( ${#bam_sm_values[@]} > 1 ))
    then
        echo "ERROR: multiple SM values found in one BAM:" >&2
        echo "${bam}" >&2
        printf '  %s\n' "${bam_sm_values[@]}" >&2
        exit 1
    fi

    bam_sm="${bam_sm_values[0]}"

    if [[ "${bam_sm}" != "${sample}" ]]
    then
        echo "ERROR: BAM SM does not match expected sample." >&2
        echo "Expected sample: ${sample}" >&2
        echo "Observed SM:     ${bam_sm}" >&2
        echo "BAM:             ${bam}" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\n' \
        "${sample}" \
        "${bam}" \
        "${bam_sm}" \
        >> "${MANIFEST}"

    # Add one -I argument for each BAM
    bam_args+=(
        -I "${bam}"
    )
done


# ============================================================
# 7. Confirm 15 unique samples
# ============================================================

n_samples="$(
    tail -n +2 "${MANIFEST}" |
    cut -f1 |
    sort -u |
    wc -l
)"

n_sm="$(
    tail -n +2 "${MANIFEST}" |
    cut -f3 |
    sort -u |
    wc -l
)"

if [[ "${n_samples}" -ne 15 ]]
then
    echo "ERROR: expected 15 sample IDs, found ${n_samples}." >&2
    exit 1
fi

if [[ "${n_sm}" -ne 15 ]]
then
    echo "ERROR: expected 15 unique BAM SM values, found ${n_sm}." >&2
    exit 1
fi


# ============================================================
# 8. Direct multi-sample HaplotypeCaller
# ============================================================

echo
echo "Running direct HaplotypeCaller on 15 offspring BAMs..."
echo "Output: ${WORK_VCF}"
echo

gatk \
    --java-options \
    "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP_DIR}" \
    HaplotypeCaller \
    -R "${REF}" \
    "${bam_args[@]}" \
    -L "${INTERVALS}" \
    -O "${WORK_VCF}" \
    -A StrandBiasBySample \
    --standard-min-confidence-threshold-for-calling \
    "${STAND_CALL_CONF}" \
    --min-base-quality-score \
    "${MIN_BASE_QUALITY}" \
    --native-pair-hmm-threads \
    "${PAIRHMM_THREADS}" \
    --tmp-dir "${TMP_DIR}" \
    > "${GATK_LOG}" 2>&1


# ============================================================
# 9. Validate output
# ============================================================

if [[ ! -s "${WORK_VCF}" ]]
then
    echo "ERROR: output VCF was not generated." >&2
    echo "See log: ${GATK_LOG}" >&2
    exit 1
fi

if [[ ! -s "${WORK_VCF}.tbi" ]]
then
    gatk IndexFeatureFile \
        -I "${WORK_VCF}"
fi

echo
echo "Samples in output VCF:"
bcftools query -l "${WORK_VCF}"

vcf_sample_count="$(
    bcftools query -l "${WORK_VCF}" |
    wc -l
)"

if [[ "${vcf_sample_count}" -ne 15 ]]
then
    echo "ERROR: output VCF contains ${vcf_sample_count} samples." >&2
    exit 1
fi

variant_count="$(
    bcftools view \
        -H \
        "${WORK_VCF}" |
    wc -l
)"

echo
echo "Direct multi-sample calling completed."
echo "Samples:  ${vcf_sample_count}"
echo "Variants: ${variant_count}"
echo "VCF:      ${WORK_VCF}"
echo "Log:      ${GATK_LOG}"
