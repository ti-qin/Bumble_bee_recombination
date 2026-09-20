from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path("/localdata/qinti/Project/Bee")

INPUT = (
    PROJECT
    / "03_variants/gatk_final/qc/annotations/"
      "biallelic_snp_annotations.tsv.gz"
)

OUT_DIR = (
    PROJECT
    / "03_variants/gatk_final/qc/annotations"
)

SUMMARY_OUT = OUT_DIR / "variant_annotation_summary.tsv"
THRESHOLD_OUT = OUT_DIR / "rank_sum_threshold_summary.tsv"


# ============================================================
# 2. Read annotations
# ============================================================

data = pd.read_csv(
    INPUT,
    sep="\t",
    na_values=[".", "nan", "NaN"]
)

annotation_columns = [
    "QUAL",
    "QD",
    "MQ",
    "FS",
    "SOR",
    "MQRankSum",
    "ReadPosRankSum",
]

for column in annotation_columns:
    data[column] = pd.to_numeric(
        data[column],
        errors="coerce"
    )

n_total = len(data)

print(f"Total biallelic SNPs: {n_total:,}")


# ============================================================
# 3. Quantile summary
# ============================================================

summary_rows = []

quantiles = [
    0.001,
    0.005,
    0.01,
    0.025,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.975,
    0.99,
    0.995,
    0.999,
]

for annotation in annotation_columns:
    values = data[annotation].dropna()

    row = {
        "annotation": annotation,
        "n_total_sites": n_total,
        "n_non_missing": len(values),
        "n_missing": n_total - len(values),
        "missing_fraction": 1 - len(values) / n_total,
        "mean": values.mean(),
        "standard_deviation": values.std(),
        "minimum": values.min(),
        "maximum": values.max(),
    }

    for quantile in quantiles:
        name = f"q{quantile:g}"
        row[name] = values.quantile(quantile)

    summary_rows.append(row)

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False
)

print("\nAnnotation summary:")
print(summary.to_string(index=False))


# ============================================================
# 4. RankSum threshold summary
# ============================================================

mq_values = data["MQRankSum"].dropna()
readpos_values = data["ReadPosRankSum"].dropna()

threshold_rows = [
    {
        "annotation": "MQRankSum",
        "threshold": -12.5,
        "direction": "less_than",
        "n_non_missing": len(mq_values),
        "n_failing": int((mq_values < -12.5).sum()),
        "fraction_failing": float((mq_values < -12.5).mean()),
    },
    {
        "annotation": "ReadPosRankSum",
        "threshold": -8.0,
        "direction": "less_than",
        "n_non_missing": len(readpos_values),
        "n_failing": int((readpos_values < -8.0).sum()),
        "fraction_failing": float(
            (readpos_values < -8.0).mean()
        ),
    },
]

threshold_summary = pd.DataFrame(threshold_rows)

threshold_summary.to_csv(
    THRESHOLD_OUT,
    sep="\t",
    index=False
)

print("\nThreshold summary:")
print(threshold_summary.to_string(index=False))

print(f"\nSaved: {SUMMARY_OUT}")
print(f"Saved: {THRESHOLD_OUT}")
