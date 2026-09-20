#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

VCF="${PROJECT}/03_variants/filtered/all_samples.diploid_for_ploidy.biallelic.snp.qual30.vcf.gz"
META="${PROJECT}/00_rawdata/metadata/clean_samples.tsv"

OUT_DIR="${PROJECT}/04_pedigree/intermediate"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${VCF}" ]]; then
    echo "ERROR: VCF not found: ${VCF}" >&2
    exit 1
fi

if [[ ! -s "${META}" ]]; then
    echo "ERROR: metadata not found: ${META}" >&2
    exit 1
fi

bcftools query -l "${VCF}" \
    > "${OUT_DIR}/vcf_samples.txt"

tail -n +2 "${META}" |
cut -f1 |
sort \
    > "${OUT_DIR}/metadata_samples.sorted.txt"

sort "${OUT_DIR}/vcf_samples.txt" \
    > "${OUT_DIR}/vcf_samples.sorted.txt"

echo "VCF sample number:"
wc -l "${OUT_DIR}/vcf_samples.txt"

echo "Metadata sample number:"
wc -l "${OUT_DIR}/metadata_samples.sorted.txt"

if ! diff -u \
    "${OUT_DIR}/metadata_samples.sorted.txt" \
    "${OUT_DIR}/vcf_samples.sorted.txt"
then
    echo "ERROR: VCF and metadata sample names are inconsistent." >&2
    exit 1
fi

echo
echo "Sample composition:"

awk -F '\t' '
NR > 1 {
    count[$3]++
}
END {
    for (role in count) {
        print role, count[role]
    }
}
' "${META}"

echo
echo "Samples by family and role:"

awk -F '\t' '
NR > 1 {
    count[$2 FS $3]++
}
END {
    for (key in count) {
        split(key, parts, FS)
        print parts[1], parts[2], count[key]
    }
}
' "${META}" |
sort

echo
echo "Pedigree input validation completed."
