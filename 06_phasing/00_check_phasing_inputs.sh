#!/usr/bin/env bash

set -euo pipefail

PROJECT="/localdata/qinti/Project/Bee"

META="${PROJECT}/00_rawdata/metadata/analysis_samples.noC3.tsv"
PRIMARY_BED="${PROJECT}/00_rawdata/reference/primary_chromosomes.bed"

if [[ ! -s "${META}" ]]; then
    echo "ERROR: metadata not found: ${META}" >&2
    exit 1
fi

if [[ ! -s "${PRIMARY_BED}" ]]; then
    echo "ERROR: primary chromosome BED not found: ${PRIMARY_BED}" >&2
    exit 1
fi

echo "========================================"
echo "Phasing input validation"
echo "========================================"

for family in C1 C2 C4
do
    marker_file="${PROJECT}/05_markers/family_markers/${family}/${family}.high_confidence_markers.tsv.gz"

    if [[ ! -s "${marker_file}" ]]; then
        echo "ERROR: marker file not found:"
        echo "${marker_file}"
        exit 1
    fi

    echo
    echo "${family}"

    echo -n "Samples in metadata: "
    awk -F '\t' -v fam="${family}" '
    NR > 1 && $2 == fam {
        count++
    }
    END {
        print count + 0
    }
    ' "${META}"

    echo -n "Offspring in metadata: "
    awk -F '\t' -v fam="${family}" '
    NR > 1 && $2 == fam && $3 == "offspring" {
        count++
    }
    END {
        print count + 0
    }
    ' "${META}"

    echo -n "Marker number: "
    zcat "${marker_file}" |
        tail -n +2 |
        wc -l

    echo "Header:"
    zcat "${marker_file}" |
        head -1 |
        tr '\t' '\n' |
        nl -ba
done

echo
echo "Input validation completed."
