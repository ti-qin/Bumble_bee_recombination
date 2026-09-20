#!/usr/bin/env bash

set -Eeuo pipefail


# ============================================================
# 1. Paths and resource parameters
# ============================================================

PROJECT="${PROJECT:-/localdata/qinti/Project/Bee}"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
INTERVALS="${PROJECT}/00_rawdata/reference/primary_chromosomes.bed"
META="${PROJECT}/00_rawdata/metadata/analysis_samples.noC3.tsv"

BAM_DIR="${PROJECT}/02_mapping/bam"

OUT_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/families"
SUMMARY_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/summary"
TMP_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/tmp"
LOG_ROOT="${PROJECT}/logs/gatk_direct_diploid_SB"

FAMILIES=(C1 C2 C4)

# Three families run concurrently by default.  These values may be
# overridden without editing this file, for example:
#
#   JAVA_MEM_PER_FAMILY=40g \
#   PAIRHMM_THREADS_PER_FAMILY=16 \
#   bash this_script.sh
#
# Total maximum resources with the defaults below are approximately
# 3 x 60 GB Java heap and 3 x 32 PairHMM threads.
JAVA_MEM_PER_FAMILY="${JAVA_MEM_PER_FAMILY:-60g}"
PAIRHMM_THREADS_PER_FAMILY="${PAIRHMM_THREADS_PER_FAMILY:-32}"
INDEX_THREADS_PER_FAMILY="${INDEX_THREADS_PER_FAMILY:-16}"

STAND_CALL_CONF="${STAND_CALL_CONF:-30.0}"
MIN_BASE_QUALITY="${MIN_BASE_QUALITY:-15}"

# false: reuse an existing family VCF if it passes every validation.
# true:  archive the existing VCF and call that family again.
FORCE_RERUN="${FORCE_RERUN:-true}"


# ============================================================
# 2. Shared-input validation
# ============================================================

for program in gatk samtools bcftools awk sort diff grep
do
    if ! command -v "${program}" >/dev/null 2>&1
    then
        echo "ERROR: ${program} was not found in PATH." >&2
        exit 1
    fi
done

REF_FAI="${REF}.fai"
REF_DICT="${REF%.*}.dict"

for file in \
    "${REF}" \
    "${REF_FAI}" \
    "${REF_DICT}" \
    "${INTERVALS}" \
    "${META}"
do
    if [[ ! -s "${file}" ]]
    then
        echo "ERROR: required file not found or empty: ${file}" >&2
        exit 1
    fi
done

mkdir -p \
    "${OUT_ROOT}" \
    "${SUMMARY_ROOT}" \
    "${TMP_ROOT}" \
    "${LOG_ROOT}"


# ============================================================
# 3. Helper functions
# ============================================================

get_family_samples() {
    local family="$1"

    awk -F '\t' -v family="${family}" '
    NR == 1 {
        next
    }
    {
        sub(/\r$/, "", $NF)
    }
    $2 == family && $3 == "queen" {
        print 0 "\t" $1
    }
    $2 == family && $3 == "offspring" {
        print 1 "\t" $1
    }
    ' "${META}" |
    sort -k1,1n -k2,2 |
    cut -f2
}


validate_family_vcf() {
    local family="$1"
    local vcf="$2"
    local expected_samples="$3"
    local summary_dir="$4"

    local observed_samples="${summary_dir}/${family}.raw_vcf.samples.list"
    local header_file="${summary_dir}/${family}.raw_vcf.header.txt"
    local tag

    if [[ ! -s "${vcf}" ]]
    then
        echo "ERROR: VCF not found or empty: ${vcf}" >&2
        return 1
    fi

    if [[ ! -s "${vcf}.tbi" && ! -s "${vcf}.csi" ]]
    then
        bcftools index \
            --threads "${INDEX_THREADS_PER_FAMILY}" \
            -t \
            "${vcf}"
    fi

    bcftools query -l "${vcf}" > "${observed_samples}"

    if ! diff -q \
        <(sort "${expected_samples}") \
        <(sort "${observed_samples}") \
        >/dev/null
    then
        echo "ERROR: ${family} VCF samples differ from metadata." >&2
        diff \
            <(sort "${expected_samples}") \
            <(sort "${observed_samples}") \
            >&2 || true
        return 1
    fi

    # Save the complete header first.  Do not use
    # `bcftools view -h | grep -q` with pipefail: grep may close the
    # pipe early and make bcftools return SIGPIPE (141).
    bcftools view -h "${vcf}" > "${header_file}"

    for tag in GT AD DP GQ SB
    do
        if ! grep -q "^##FORMAT=<ID=${tag}," "${header_file}"
        then
            echo "ERROR: FORMAT/${tag} is absent from ${family} VCF." >&2
            return 1
        fi
    done

    for tag in QD FS SOR MQ MQRankSum ReadPosRankSum
    do
        if ! grep -q "^##INFO=<ID=${tag}," "${header_file}"
        then
            echo "ERROR: INFO/${tag} is absent from ${family} VCF." >&2
            return 1
        fi
    done

    # Read the pipeline to completion.  The former implementation used
    # `awk ... { exit }`; under `set -o pipefail` that gave bcftools a
    # SIGPIPE and silently terminated the parent script after C1.
    if ! bcftools query -f '[%SB\n]' "${vcf}" |
        awk '
        $0 != "." && $0 != "" {
            found = 1
        }
        END {
            exit(found ? 0 : 1)
        }
        '
    then
        echo "ERROR: no non-missing FORMAT/SB value was found for ${family}." >&2
        return 1
    fi

    return 0
}


archive_existing_vcf() {
    local family="$1"
    local vcf="$2"
    local summary_dir="$3"
    local timestamp
    local archive_dir
    local path

    if [[ ! -e "${vcf}" && ! -e "${vcf}.tbi" && ! -e "${vcf}.csi" ]]
    then
        return 0
    fi

    timestamp="$(date +%Y%m%d_%H%M%S)"
    archive_dir="${summary_dir}/previous_call_${timestamp}"
    mkdir -p "${archive_dir}"

    for path in "${vcf}" "${vcf}.tbi" "${vcf}.csi"
    do
        if [[ -e "${path}" ]]
        then
            mv "${path}" "${archive_dir}/"
        fi
    done

    echo "Archived the previous ${family} VCF under: ${archive_dir}"
}


run_family() (
    set -Eeuo pipefail

    local family="$1"
    local family_root="${OUT_ROOT}/${family}"
    local raw_dir="${family_root}/raw"
    local family_summary_dir="${SUMMARY_ROOT}/${family}"
    local family_tmp="${TMP_ROOT}/${family}/haplotypecaller"
    local family_log_dir="${LOG_ROOT}/${family}"

    local raw_vcf="${raw_dir}/${family}.direct_diploid_SB.raw.vcf.gz"
    local work_vcf="${raw_dir}/${family}.direct_diploid_SB.raw.inprogress.vcf.gz"
    local sample_list="${family_summary_dir}/${family}.samples.list"
    local bam_list="${family_summary_dir}/${family}.bams.list"
    local parameter_file="${family_summary_dir}/${family}.calling_parameters.tsv"
    local gatk_log="${family_log_dir}/01_haplotypecaller_direct_diploid_SB.log"
    local success_file="${family_summary_dir}/${family}.01_calling.SUCCESS"
    local failed_file="${family_summary_dir}/${family}.01_calling.FAILED"

    local queen_count
    local offspring_count
    local sample_id
    local bam
    local bam_sm
    local failure_status=1

    mkdir -p \
        "${raw_dir}" \
        "${family_summary_dir}" \
        "${family_tmp}" \
        "${family_log_dir}"

    rm -f "${failed_file}"

    trap '
        failure_status=$?
        printf "family\t%s\nexit_status\t%s\ntime\t%s\n" \
            "${family}" "${failure_status}" "$(date -Is)" \
            > "${failed_file}"
        echo "ERROR: ${family} failed with status ${failure_status}." >&2
        echo "GATK log: ${gatk_log}" >&2
        exit "${failure_status}"
    ' ERR

    mapfile -t samples < <(get_family_samples "${family}")

    queen_count="$(
        awk -F '\t' -v family="${family}" '
        NR > 1 && $2 == family && $3 == "queen" {n++}
        END {print n + 0}
        ' "${META}"
    )"

    offspring_count="$(
        awk -F '\t' -v family="${family}" '
        NR > 1 && $2 == family && $3 == "offspring" {n++}
        END {print n + 0}
        ' "${META}"
    )"

    if [[ "${#samples[@]}" -ne 6 || "${queen_count}" -ne 1 || "${offspring_count}" -ne 5 ]]
    then
        echo "ERROR: ${family} must contain exactly 1 queen and 5 offspring." >&2
        false
    fi

    printf "%s\n" "${samples[@]}" > "${sample_list}"
    : > "${bam_list}"

    bam_args=()

    for sample_id in "${samples[@]}"
    do
        bam="${BAM_DIR}/${sample_id}.markdup.bam"

        if [[ ! -s "${bam}" ]]
        then
            echo "ERROR: BAM not found: ${bam}" >&2
            false
        fi

        samtools quickcheck -v "${bam}"

        if [[ ! -s "${bam}.bai" && ! -s "${bam%.bam}.bai" ]]
        then
            samtools index \
                -@ "${INDEX_THREADS_PER_FAMILY}" \
                "${bam}"
        fi

        bam_sm="$(
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
        )"

        if [[ "${bam_sm}" != "${sample_id}" ]]
        then
            echo "ERROR: ${bam} has SM=${bam_sm}; expected ${sample_id}." >&2
            false
        fi

        printf "%s\n" "${bam}" >> "${bam_list}"
        bam_args+=( -I "${bam}" )
    done

    if [[ "${FORCE_RERUN}" != "true" ]] &&
        validate_family_vcf \
            "${family}" \
            "${raw_vcf}" \
            "${sample_list}" \
            "${family_summary_dir}"
    then
        printf "family\t%s\nstatus\treused_existing_valid_vcf\ntime\t%s\n" \
            "${family}" "$(date -Is)" > "${success_file}"
        rm -f "${failed_file}"
        trap - ERR
        echo "${family}: existing validated VCF reused."
        exit 0
    fi

    archive_existing_vcf \
        "${family}" \
        "${raw_vcf}" \
        "${family_summary_dir}"

    rm -f \
        "${work_vcf}" \
        "${work_vcf}.tbi" \
        "${work_vcf}.csi" \
        "${success_file}"

    echo "${family}: starting HaplotypeCaller with ${#samples[@]} BAMs."
    echo "${family}: GATK log: ${gatk_log}"

    gatk \
        --java-options \
        "-Xmx${JAVA_MEM_PER_FAMILY} -Djava.io.tmpdir=${family_tmp}" \
        HaplotypeCaller \
        -R "${REF}" \
        "${bam_args[@]}" \
        -L "${INTERVALS}" \
        -O "${work_vcf}" \
        -A StrandBiasBySample \
        --standard-min-confidence-threshold-for-calling \
        "${STAND_CALL_CONF}" \
        --min-base-quality-score \
        "${MIN_BASE_QUALITY}" \
        --native-pair-hmm-threads \
        "${PAIRHMM_THREADS_PER_FAMILY}" \
        > "${gatk_log}" 2>&1

    validate_family_vcf \
        "${family}" \
        "${work_vcf}" \
        "${sample_list}" \
        "${family_summary_dir}"

    mv "${work_vcf}" "${raw_vcf}"

    if [[ -e "${work_vcf}.tbi" ]]
    then
        mv "${work_vcf}.tbi" "${raw_vcf}.tbi"
    elif [[ -e "${work_vcf}.csi" ]]
    then
        mv "${work_vcf}.csi" "${raw_vcf}.csi"
    else
        bcftools index \
            --threads "${INDEX_THREADS_PER_FAMILY}" \
            -t \
            "${raw_vcf}"
    fi

    bcftools stats "${raw_vcf}" \
        > "${family_summary_dir}/${family}.direct_diploid_SB.raw.stats.txt"

    {
        printf "parameter\tvalue\n"
        printf "family_id\t%s\n" "${family}"
        printf "n_samples\t%s\n" "${#samples[@]}"
        printf "sample_ploidy\tdefault_2_for_all_samples\n"
        printf "emit_ref_confidence\tNONE\n"
        printf "StrandBiasBySample\tenabled\n"
        printf "stand_call_conf\t%s\n" "${STAND_CALL_CONF}"
        printf "min_base_quality_score\t%s\n" "${MIN_BASE_QUALITY}"
        printf "native_pair_hmm_threads\t%s\n" "${PAIRHMM_THREADS_PER_FAMILY}"
        printf "java_heap_per_family\t%s\n" "${JAVA_MEM_PER_FAMILY}"
        printf "raw_vcf\t%s\n" "${raw_vcf}"
    } > "${parameter_file}"

    printf "family\t%s\nstatus\tcompleted\ntime\t%s\nvcf\t%s\n" \
        "${family}" "$(date -Is)" "${raw_vcf}" > "${success_file}"

    rm -f "${failed_file}"
    trap - ERR

    echo "${family}: completed successfully."
    echo "${family}: ${raw_vcf}"
)


# ============================================================
# 4. Launch all families concurrently
# ============================================================

declare -A FAMILY_PID=()
declare -A FAMILY_WRAPPER_LOG=()

echo "Launching family-specific HaplotypeCaller jobs in parallel."
echo "Families: ${FAMILIES[*]}"
echo "Java heap per family: ${JAVA_MEM_PER_FAMILY}"
echo "PairHMM threads per family: ${PAIRHMM_THREADS_PER_FAMILY}"
echo

for family in "${FAMILIES[@]}"
do
    wrapper_log="${LOG_ROOT}/${family}/01_family_parallel_wrapper.log"
    mkdir -p "$(dirname "${wrapper_log}")"

    run_family "${family}" > "${wrapper_log}" 2>&1 &

    FAMILY_PID["${family}"]=$!
    FAMILY_WRAPPER_LOG["${family}"]="${wrapper_log}"

    echo "Launched ${family}: PID ${FAMILY_PID[${family}]}"
    echo "  wrapper log: ${wrapper_log}"
done


# ============================================================
# 5. Wait for every family and report all outcomes
# ============================================================

failed_families=()

for family in "${FAMILIES[@]}"
do
    pid="${FAMILY_PID[${family}]}"

    # `wait` is deliberately inside an if statement.  A failed family
    # therefore does not trigger the parent shell's `set -e`; the other
    # background jobs continue and all outcomes are collected.
    if wait "${pid}"
    then
        echo "SUCCESS: ${family}"
    else
        status=$?
        echo "FAILED: ${family} (exit status ${status})" >&2
        echo "  wrapper log: ${FAMILY_WRAPPER_LOG[${family}]}" >&2
        failed_families+=("${family}")
    fi
done

echo

if [[ "${#failed_families[@]}" -gt 0 ]]
then
    echo "ERROR: failed families: ${failed_families[*]}" >&2
    echo "Successful families were retained and do not need to be rerun." >&2
    exit 1
fi

echo "All family-specific direct HaplotypeCaller runs completed successfully."
