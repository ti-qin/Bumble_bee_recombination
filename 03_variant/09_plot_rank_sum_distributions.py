from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path("/localdata/qinti/Project/Bee")

INPUT = (
    PROJECT
    / "03_variants/gatk_final/qc/annotations/"
      "biallelic_snp_annotations.tsv.gz"
)

FIG_DIR = (
    PROJECT
    / "03_variants/gatk_final/qc/figures"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. Read data
# ============================================================

data = pd.read_csv(
    INPUT,
    sep="\t",
    na_values=[".", "nan", "NaN"]
)

data["MQRankSum"] = pd.to_numeric(
    data["MQRankSum"],
    errors="coerce"
)

data["ReadPosRankSum"] = pd.to_numeric(
    data["ReadPosRankSum"],
    errors="coerce"
)

mq = data["MQRankSum"].dropna().to_numpy()
readpos = data["ReadPosRankSum"].dropna().to_numpy()

print(f"MQRankSum non-missing sites: {len(mq):,}")
print(
    f"ReadPosRankSum non-missing sites: "
    f"{len(readpos):,}"
)


# ============================================================
# 3. MQRankSum histogram
# ============================================================

mq_lower = np.quantile(mq, 0.001)
mq_upper = np.quantile(mq, 0.999)

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.hist(
    mq,
    bins=100,
    range=(mq_lower, mq_upper),
    density=True
)

ax.axvline(
    -12.5,
    linestyle="--",
    linewidth=1,
    label="GATK starting threshold: −12.5"
)

ax.axvline(
    0,
    linestyle=":",
    linewidth=0.8
)

ax.set_xlabel("MQRankSum")
ax.set_ylabel("Density")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "MQRankSum_histogram.pdf"
)

plt.close(fig)


# ============================================================
# 4. ReadPosRankSum histogram
# ============================================================

readpos_lower = np.quantile(
    readpos,
    0.001
)

readpos_upper = np.quantile(
    readpos,
    0.999
)

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.hist(
    readpos,
    bins=100,
    range=(readpos_lower, readpos_upper),
    density=True
)

ax.axvline(
    -8.0,
    linestyle="--",
    linewidth=1,
    label="GATK starting threshold: −8"
)

ax.axvline(
    0,
    linestyle=":",
    linewidth=0.8
)

ax.set_xlabel("ReadPosRankSum")
ax.set_ylabel("Density")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "ReadPosRankSum_histogram.pdf"
)

plt.close(fig)


# ============================================================
# 5. MQRankSum ECDF
# ============================================================

mq_sorted = np.sort(mq)

mq_ecdf = (
    np.arange(1, len(mq_sorted) + 1)
    / len(mq_sorted)
)

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.plot(
    mq_sorted,
    mq_ecdf,
    linewidth=1
)

ax.axvline(
    -12.5,
    linestyle="--",
    linewidth=1,
    label="−12.5"
)

ax.axvline(
    0,
    linestyle=":",
    linewidth=0.8
)

ax.set_xlim(
    mq_lower,
    mq_upper
)

ax.set_xlabel("MQRankSum")
ax.set_ylabel("Cumulative proportion")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "MQRankSum_ECDF.pdf"
)

plt.close(fig)


# ============================================================
# 6. ReadPosRankSum ECDF
# ============================================================

readpos_sorted = np.sort(
    readpos
)

readpos_ecdf = (
    np.arange(1, len(readpos_sorted) + 1)
    / len(readpos_sorted)
)

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.plot(
    readpos_sorted,
    readpos_ecdf,
    linewidth=1
)

ax.axvline(
    -8.0,
    linestyle="--",
    linewidth=1,
    label="−8"
)

ax.axvline(
    0,
    linestyle=":",
    linewidth=0.8
)

ax.set_xlim(
    readpos_lower,
    readpos_upper
)

ax.set_xlabel("ReadPosRankSum")
ax.set_ylabel("Cumulative proportion")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "ReadPosRankSum_ECDF.pdf"
)

plt.close(fig)


# ============================================================
# 7. Left-tail zoom for MQRankSum
# ============================================================

mq_left_limit = np.quantile(
    mq,
    0.10
)

mq_left = mq[mq <= mq_left_limit]

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.hist(
    mq_left,
    bins=80,
    density=True
)

ax.axvline(
    -12.5,
    linestyle="--",
    linewidth=1,
    label="−12.5"
)

ax.set_xlabel("MQRankSum, lower 10%")
ax.set_ylabel("Density")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "MQRankSum_left_tail.pdf"
)

plt.close(fig)


# ============================================================
# 8. Left-tail zoom for ReadPosRankSum
# ============================================================

readpos_left_limit = np.quantile(
    readpos,
    0.10
)

readpos_left = readpos[
    readpos <= readpos_left_limit
]

fig, ax = plt.subplots(
    figsize=(3.5, 2.8),
    dpi=150
)

ax.hist(
    readpos_left,
    bins=80,
    density=True
)

ax.axvline(
    -8.0,
    linestyle="--",
    linewidth=1,
    label="−8"
)

ax.set_xlabel("ReadPosRankSum, lower 10%")
ax.set_ylabel("Density")
ax.legend(
    frameon=False,
    fontsize=6
)

fig.tight_layout()

fig.savefig(
    FIG_DIR / "ReadPosRankSum_left_tail.pdf"
)

plt.close(fig)

print(f"Figures saved under: {FIG_DIR}")
