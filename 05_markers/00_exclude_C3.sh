#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

META_IN="${PROJECT}/00_rawdata/metadata/clean_samples.tsv"
META_OUT="${PROJECT}/00_rawdata/metadata/analysis_samples.noC3.tsv"
EXCLUDED="${PROJECT}/05_markers/exclusion/excluded_samples.tsv"

if [[ ! -s "${META_IN}" ]]; then
    echo "ERROR: metadata not found: ${META_IN}" >&2
    exit 1
fi

# 保留 C1、C2 和 C4
awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    print
    next
}
$2 != "C3" {
    print
}
' "${META_IN}" > "${META_OUT}"

# 记录排除样本和原因
printf "sample_id\tfamily_id\trole\texclusion_stage\treason\n" \
    > "${EXCLUDED}"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR > 1 && $2 == "C3" {
    print $1, $2, $3, "pedigree_validation", \
    "C3_Q was genetically incompatible with the C3 offspring cluster"
}
' "${META_IN}" >> "${EXCLUDED}"

echo "Retained samples:"
tail -n +2 "${META_OUT}" | wc -l

echo "Retained queens:"
awk -F '\t' 'NR > 1 && $3 == "queen"' "${META_OUT}" | wc -l

echo "Retained offspring:"
awk -F '\t' 'NR > 1 && $3 == "offspring"' "${META_OUT}" | wc -l

echo
echo "Excluded samples:"
column -t -s $'\t' "${EXCLUDED}"
