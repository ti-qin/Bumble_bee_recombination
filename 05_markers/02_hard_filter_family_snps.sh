#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# 1. Paths and parameters
# ============================================================

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"

OUT_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/families"
SUMMARY_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/summary"
TMP_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/tmp"
LOG_ROOT="${PROJECT}/logs/gatk_direct_diploid_SB"

FAMILIES=(C1 C2 C4)

THREADS=64
JAVA_MEM="32g"

COMBINED_COUNTS="${SUMMARY_ROOT}/family_variant_counts_by_stage.tsv"

mkdir -p \
    "${SUMMARY_ROOT}" \
    "${TMP_ROOT}" \
    "${LOG_ROOT}"


# ============================================================
# 2. Locate programs and validate reference
# ============================================================

for program in gatk bcftools
do
    if ! command -v "${program}" >/dev/null 2>&1
    then
        echo "ERROR: ${program} was not found in PATH." >&2
        exit 1
    fi
done

if [[ ! -s "${REF}" ]]
then
    echo "ERROR: reference not found: ${REF}" >&2
    exit 1
fi

printf "family_id\tstage\tn_records\tvcf\n" \
    > "${COMBINED_COUNTS}"


# ============================================================
# 3. Filter each family independently
# ============================================================

count_records() {
    bcftools view -H "$1" | wc -l
}

for family in "${FAMILIES[@]}"
do
    FAMILY_ROOT="${OUT_ROOT}/${family}"
    RAW="${FAMILY_ROOT}/raw/${family}.direct_diploid_SB.raw.vcf.gz"

    FILTER_DIR="${FAMILY_ROOT}/filtered"
    FAMILY_SUMMARY_DIR="${SUMMARY_ROOT}/${family}"
    FAMILY_TMP="${TMP_ROOT}/${family}/filtering"
    FAMILY_LOG_DIR="${LOG_ROOT}/${family}"

    BIALLELIC="${FILTER_DIR}/${family}.direct_diploid_SB.biallelic.snps.raw.vcf.gz"
    FILTERED="${FILTER_DIR}/${family}.direct_diploid_SB.biallelic.snps.filtered.vcf.gz"
    PASS="${FILTER_DIR}/${family}.direct_diploid_SB.biallelic.snps.PASS.vcf.gz"

    LOG="${FAMILY_LOG_DIR}/02_hard_filter_family_snps.log"

    mkdir -p \
        "${FILTER_DIR}" \
        "${FAMILY_SUMMARY_DIR}" \
        "${FAMILY_TMP}" \
        "${FAMILY_LOG_DIR}"

    if [[ ! -s "${RAW}" ]]
    then
        echo "ERROR: ${family} raw VCF not found: ${RAW}" >&2
        exit 1
    fi

    if ! bcftools view -h "${RAW}" |
        grep -q '^##FORMAT=<ID=SB,'
    then
        echo "ERROR: ${family} raw VCF does not contain FORMAT/SB." >&2
        exit 1
    fi

    echo
    echo "Filtering family ${family}"

    # --------------------------------------------------------
    # 3a. Keep biallelic SNPs only
    # --------------------------------------------------------

    rm -f \
        "${BIALLELIC}" \
        "${BIALLELIC}.tbi" \
        "${BIALLELIC}.csi"

    bcftools view \
        --threads "${THREADS}" \
        -m2 \
        -M2 \
        -v snps \
        -Oz \
        -o "${BIALLELIC}" \
        "${RAW}"

    bcftools index \
        --threads "${THREADS}" \
        -t \
        "${BIALLELIC}"

    # --------------------------------------------------------
    # 3b. Apply the existing GATK SNP hard filters
    # --------------------------------------------------------

    rm -f \
        "${FILTERED}" \
        "${FILTERED}.tbi" \
        "${FILTERED}.csi"

    gatk \
        --java-options \
        "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${FAMILY_TMP}" \
        VariantFiltration \
        -R "${REF}" \
        -V "${BIALLELIC}" \
        -O "${FILTERED}" \
        --filter-name "QD2" \
        --filter-expression "QD < 2.0" \
        --filter-name "FS60" \
        --filter-expression "FS > 60.0" \
        --filter-name "SOR3" \
        --filter-expression "SOR > 3.0" \
        --filter-name "MQ40" \
        --filter-expression "MQ < 40.0" \
        --filter-name "MQRankSum-12.5" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "ReadPosRankSum-8" \
        --filter-expression "ReadPosRankSum < -8.0" \
        > "${LOG}" 2>&1

    if [[ ! -s "${FILTERED}.tbi" && ! -s "${FILTERED}.csi" ]]
    then
        bcftools index \
            --threads "${THREADS}" \
            -t \
            "${FILTERED}"
    fi

    # --------------------------------------------------------
    # 3c. Extract PASS records
    # --------------------------------------------------------

    rm -f \
        "${PASS}" \
        "${PASS}.tbi" \
        "${PASS}.csi"

    bcftools view \
        --threads "${THREADS}" \
        -f PASS \
        -Oz \
        -o "${PASS}" \
        "${FILTERED}"

    bcftools index \
        --threads "${THREADS}" \
        -t \
        "${PASS}"

    # --------------------------------------------------------
    # 3d. Validate the final family VCF
    # --------------------------------------------------------

    for tag in GT AD DP GQ SB
    do
        if ! bcftools view -h "${PASS}" |
            grep -q "^##FORMAT=<ID=${tag},"
        then
            echo "ERROR: FORMAT/${tag} is absent from ${family} PASS VCF." >&2
            exit 1
        fi
    done

    n_raw="$(count_records "${RAW}")"
    n_biallelic="$(count_records "${BIALLELIC}")"
    n_filtered="$(count_records "${FILTERED}")"
    n_pass="$(count_records "${PASS}")"

    {
        printf "stage\tn_records\tvcf\n"
        printf "raw_direct_multisample\t%s\t%s\n" \
            "${n_raw}" "${RAW}"
        printf "biallelic_snps\t%s\t%s\n" \
            "${n_biallelic}" "${BIALLELIC}"
        printf "hard_filtered_all\t%s\t%s\n" \
            "${n_filtered}" "${FILTERED}"
        printf "biallelic_snps_PASS\t%s\t%s\n" \
            "${n_pass}" "${PASS}"
    } > "${FAMILY_SUMMARY_DIR}/${family}.variant_counts_by_stage.tsv"

    printf "%s\traw_direct_multisample\t%s\t%s\n" \
        "${family}" "${n_raw}" "${RAW}" \
        >> "${COMBINED_COUNTS}"

    printf "%s\tbiallelic_snps\t%s\t%s\n" \
        "${family}" "${n_biallelic}" "${BIALLELIC}" \
        >> "${COMBINED_COUNTS}"

    printf "%s\thard_filtered_all\t%s\t%s\n" \
        "${family}" "${n_filtered}" "${FILTERED}" \
        >> "${COMBINED_COUNTS}"

    printf "%s\tbiallelic_snps_PASS\t%s\t%s\n" \
        "${family}" "${n_pass}" "${PASS}" \
        >> "${COMBINED_COUNTS}"

    {
        printf "filter\tn_records\n"

        bcftools query \
            -f '%FILTER\n' \
            "${FILTERED}" |
        awk '
        {
            n = split($0, values, ";")
            for (i = 1; i <= n; i++) {
                count[values[i]]++
            }
        }
        END {
            for (key in count) {
                print key "\t" count[key]
            }
        }
        ' |
        sort -k1,1
    } > "${FAMILY_SUMMARY_DIR}/${family}.hard_filter_counts.tsv"

    bcftools query \
        -f '[%SAMPLE\t%GT\n]' \
        "${PASS}" |
    awk '
    BEGIN {
        OFS = "\t"
    }
    {
        count[$1, $2]++
    }
    END {
        print "sample_id", "GT", "n_calls"
        for (key in count) {
            split(key, values, SUBSEP)
            print values[1], values[2], count[key]
        }
    }
    ' |
    sort -k1,1 -k2,2 \
    > "${FAMILY_SUMMARY_DIR}/${family}.PASS_genotype_counts.tsv"

    bcftools stats \
        "${PASS}" \
        > "${FAMILY_SUMMARY_DIR}/${family}.direct_diploid_SB.biallelic.snps.PASS.stats.txt"

    echo "${family} final PASS SNP VCF:"
    echo "${PASS}"
    echo "PASS SNP count: ${n_pass}"
done

echo
echo "All family-specific SNP filtering completed."
echo "Combined count summary:"
echo "${COMBINED_COUNTS}"

