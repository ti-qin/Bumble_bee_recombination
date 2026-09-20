from pathlib import Path
import gzip
import re

import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path("/localdata/qinti/Project/Bee")

REFERENCE = (
    PROJECT
    / "00_rawdata/reference/"
      "Bombus_vestalis.fa"
)

PRIMARY_LIST = (
    PROJECT
    / "00_rawdata/reference/"
      "primary_chromosomes.txt"
)

OUT_DIR = (
    PROJECT
    / "07_recombination/reference_masks"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True
)

GAP_BED = OUT_DIR / "reference_N_gaps.bed"
GAP_SUMMARY = OUT_DIR / "reference_N_gap_summary.tsv"


# ============================================================
# 2. Helpers
# ============================================================

def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")

    return path.open("r")


primary_chromosomes = {
    line.strip()
    for line in PRIMARY_LIST.open()
    if line.strip()
}


# ============================================================
# 3. Parse FASTA
# ============================================================

gap_rows = []

current_name = None
sequence_parts = []


def process_sequence(
    sequence_name,
    parts,
):
    if sequence_name is None:
        return

    if sequence_name not in primary_chromosomes:
        return

    sequence = "".join(parts)

    for match in re.finditer(
        r"[Nn]+",
        sequence
    ):
        # BED is zero-based and half-open
        gap_rows.append(
            {
                "chromosome": sequence_name,
                "start": match.start(),
                "end": match.end(),
                "length": (
                    match.end()
                    - match.start()
                ),
            }
        )


with open_text(REFERENCE) as handle:
    for line in handle:
        line = line.rstrip("\n")

        if line.startswith(">"):
            process_sequence(
                current_name,
                sequence_parts,
            )

            current_name = (
                line[1:]
                .split()[0]
            )

            sequence_parts = []
        else:
            if (
                current_name
                in primary_chromosomes
            ):
                sequence_parts.append(
                    line.strip()
                )

process_sequence(
    current_name,
    sequence_parts,
)


# ============================================================
# 4. Save BED and summary
# ============================================================

gaps = pd.DataFrame(gap_rows)

if gaps.empty:
    GAP_BED.write_text("")

    summary = pd.DataFrame(
        columns=[
            "chromosome",
            "n_gaps",
            "total_gap_bp",
            "maximum_gap_bp",
        ]
    )
else:
    gaps = gaps.sort_values(
        [
            "chromosome",
            "start",
        ]
    )

    gaps[
        [
            "chromosome",
            "start",
            "end",
        ]
    ].to_csv(
        GAP_BED,
        sep="\t",
        header=False,
        index=False,
    )

    summary = (
        gaps.groupby(
            "chromosome",
            as_index=False
        )
        .agg(
            n_gaps=(
                "length",
                "size"
            ),
            total_gap_bp=(
                "length",
                "sum"
            ),
            maximum_gap_bp=(
                "length",
                "max"
            ),
        )
    )

summary.to_csv(
    GAP_SUMMARY,
    sep="\t",
    index=False,
)

print(f"Reference gaps: {len(gaps):,}")
print(f"Saved: {GAP_BED}")
print(f"Saved: {GAP_SUMMARY}")
