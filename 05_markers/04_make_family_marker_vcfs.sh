#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# 1. Paths and parameters
# ============================================================

PROJECT="/localdata/qinti/Project/Bee"

VCF_ROOT="${PROJECT}/03_variants/gatk_direct_diploid_SB/families"
META="${PROJECT}/00_rawdata/metadata/analysis_samples.noC3.tsv"

FAMILIES=(C1 C2 C4)
THREADS=32


# ============================================================
# 2. Locate bcftools and validate metadata
# ============================================================

if ! command -v bcftools >/dev/null 2>&1
then
    echo "ERROR: bcftools was not found in PATH." >&2
    exit 1
fi

if [[ ! -s "${META}" ]]
then
    echo "ERROR: metadata not found: ${META}" >&2
    exit 1
fi


# ============================================================
# 3. Create one marker-only VCF per family
# ============================================================

for family in "${FAMILIES[@]}"
do
    FAMILY_DIR="${PROJECT}/05_markers/family_markers/${family}"

    SOURCE_VCF="${VCF_ROOT}/${family}/filtered/${family}.direct_diploid_SB.biallelic.snps.PASS.vcf.gz"

    BED="${FAMILY_DIR}/${family}.high_confidence_markers.bed"
    SAMPLE_LIST="${FAMILY_DIR}/${family}.samples.txt"
    OUT="${FAMILY_DIR}/${family}.high_confidence_markers.vcf.gz"
    STATS="${FAMILY_DIR}/${family}.high_confidence_markers.stats.txt"

    mkdir -p "${FAMILY_DIR}"

    if [[ ! -s "${SOURCE_VCF}" ]]
    then
        echo "ERROR: source family VCF not found: ${SOURCE_VCF}" >&2
        exit 1
    fi

    if [[ ! -s "${BED}" ]]
    then
        echo "ERROR: marker BED not found or empty: ${BED}" >&2
        exit 1
    fi

    awk -F '\t' -v family="${family}" '
    NR > 1 && $2 == family && $3 == "queen" {
        print 0 "\t" $1
    }
    NR > 1 && $2 == family && $3 == "offspring" {
        print 1 "\t" $1
    }
    ' "${META}" |
    sort -k1,1n -k2,2 |
    cut -f2 \
    > "${SAMPLE_LIST}"

    if [[ "$(wc -l < "${SAMPLE_LIST}")" -ne 6 ]]
    then
        echo "ERROR: ${family} sample list should contain six samples." >&2
        exit 1
    fi

    rm -f \
        "${OUT}" \
        "${OUT}.tbi" \
        "${OUT}.csi"

    bcftools view \
        --threads "${THREADS}" \
        -S "${SAMPLE_LIST}" \
        -R "${BED}" \
        -Oz \
        -o "${OUT}" \
        "${SOURCE_VCF}"

    bcftools index \
        --threads "${THREADS}" \
        -t \
        "${OUT}"

    expected_markers="$(
        awk '
        NF >= 3 {
            key = $1 ":" $2 ":" $3
            seen[key] = 1
        }
        END {
            print length(seen)
        }
        ' "${BED}"
    )"

    observed_markers="$(
        bcftools view -H "${OUT}" |
        wc -l
    )"

    observed_samples="$(
        bcftools query -l "${OUT}" |
        wc -l
    )"

    if [[ "${observed_samples}" -ne 6 ]]
    then
        echo "ERROR: ${family} marker VCF contains ${observed_samples} samples instead of six." >&2
        exit 1
    fi

    if [[ "${observed_markers}" -ne "${expected_markers}" ]]
    then
        echo "ERROR: ${family} marker count mismatch." >&2
        echo "  expected from BED: ${expected_markers}" >&2
        echo "  observed in VCF: ${observed_markers}" >&2
        exit 1
    fi

    for tag in GT AD DP GQ SB
    do
        if ! bcftools view -h "${OUT}" |
            grep -q "^##FORMAT=<ID=${tag},"
        then
            echo "ERROR: FORMAT/${tag} is absent from ${family} marker VCF." >&2
            exit 1
        fi
    done

    bcftools stats \
        "${OUT}" \
        > "${STATS}"

    echo
    echo "${family}:"
    echo "  source VCF = ${SOURCE_VCF}"
    echo "  samples    = ${observed_samples}"
    echo "  markers    = ${observed_markers}"
    echo "  output     = ${OUT}"
done

echo
echo "All family marker VCFs completed."

