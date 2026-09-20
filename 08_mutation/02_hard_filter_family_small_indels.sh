#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

REF="${PROJECT}/00_rawdata/reference/Bombus_vestalis.fa"
ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/families"

FAMILIES=(C1 C2 C4)

JAVA_MEM="64g"
THREADS=128

for family in "${FAMILIES[@]}"
do
    FAMILY_ROOT="${ROOT}/${family}"
    RAW="${FAMILY_ROOT}/raw/${family}.direct_diploid_SB.raw.vcf.gz"

    OUT="${FAMILY_ROOT}/filtered"
    TMP="${PROJECT}/03_variants/gatk_direct_diploid_SB/tmp/${family}/indel_filter"

    SELECTED="${OUT}/${family}.direct_diploid_SB.biallelic.small_indels.raw.vcf.gz"
    NORMALIZED="${OUT}/${family}.direct_diploid_SB.biallelic.small_indels.normalized.vcf.gz"
    FILTERED="${OUT}/${family}.direct_diploid_SB.biallelic.small_indels.filtered.vcf.gz"
    PASS="${OUT}/${family}.direct_diploid_SB.biallelic.small_indels.PASS.vcf.gz"

    mkdir -p "${OUT}" "${TMP}"

    if [[ ! -s "${RAW}" ]]
    then
        echo "ERROR: missing raw VCF: ${RAW}" >&2
        exit 1
    fi

    echo
    echo "Processing ${family}"

    # --------------------------------------------------------
    # 1. Select biallelic INDELs of length 1–19 bp
    # --------------------------------------------------------

    gatk \
        --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP}" \
        SelectVariants \
        -R "${REF}" \
        -V "${RAW}" \
        --select-type-to-include INDEL \
        --restrict-alleles-to BIALLELIC \
        --min-indel-size 1 \
        --max-indel-size 19 \
        -O "${SELECTED}"

    # --------------------------------------------------------
    # 2. Left-align and trim INDEL representation
    # --------------------------------------------------------

    gatk \
        --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP}" \
        LeftAlignAndTrimVariants \
        -R "${REF}" \
        -V "${SELECTED}" \
        -O "${NORMALIZED}"

    # --------------------------------------------------------
    # 3. INDEL-specific hard filtering
    # --------------------------------------------------------

    gatk \
        --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${TMP}" \
        VariantFiltration \
        -R "${REF}" \
        -V "${NORMALIZED}" \
        --filter-name "QD2" \
        --filter-expression "QD < 2.0" \
        --filter-name "QUAL30" \
        --filter-expression "QUAL < 30.0" \
        --filter-name "FS200" \
        --filter-expression "FS > 200.0" \
        --filter-name "ReadPosRankSum-20" \
        --filter-expression "ReadPosRankSum < -20.0" \
        -O "${FILTERED}"

    # --------------------------------------------------------
    # 4. Extract PASS
    # --------------------------------------------------------

    bcftools view \
        --threads "${THREADS}" \
        -f PASS \
        -Oz \
        -o "${PASS}" \
        "${FILTERED}"

    bcftools index \
        --threads "${THREADS}" \
        -f \
        -t \
        "${PASS}"

    # --------------------------------------------------------
    # 5. Validate FORMAT fields
    # --------------------------------------------------------

    for tag in GT AD DP GQ SB
    do
        if ! bcftools view -h "${PASS}" |
            grep -q "^##FORMAT=<ID=${tag},"
        then
            echo "ERROR: FORMAT/${tag} missing from ${PASS}" >&2
            exit 1
        fi
    done

    echo "${family} counts:"

    printf "selected\t"
    bcftools view -H "${SELECTED}" | wc -l

    printf "PASS\t"
    bcftools view -H "${PASS}" | wc -l

    echo "Final VCF: ${PASS}"
done

echo
echo "All family small-INDEL filtering completed."
