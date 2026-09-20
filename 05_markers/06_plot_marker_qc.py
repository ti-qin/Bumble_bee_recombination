from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


PROJECT = Path("/localdata/qinti/Project/Bee")

SUMMARY_DIR = PROJECT / "05_markers/summary"
FIG_DIR = PROJECT / "05_markers/figures"

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)

chromosome_summary = pd.read_csv(
    SUMMARY_DIR
    / "family_marker_chromosome_summary.tsv",
    sep="\t"
)

missingness = pd.read_csv(
    SUMMARY_DIR
    / "offspring_marker_missingness.tsv",
    sep="\t"
)


# ============================================================
# 1. Marker number per chromosome
# ============================================================

marker_matrix = chromosome_summary.pivot(
    index="chromosome",
    columns="family_id",
    values="n_markers"
)

fig, ax = plt.subplots(
    figsize=(7.2, 3.2),
    dpi=150
)

marker_matrix.plot(
    kind="bar",
    ax=ax,
    width=0.8
)

ax.set_xlabel("Chromosome")
ax.set_ylabel("Number of markers")
ax.legend(
    title="Family",
    frameon=False
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "marker_number_per_chromosome.pdf"
)

plt.close(fig)


# ============================================================
# 2. Marker density per chromosome
# ============================================================

density_matrix = chromosome_summary.pivot(
    index="chromosome",
    columns="family_id",
    values="markers_per_Mb"
)

fig, ax = plt.subplots(
    figsize=(7.2, 3.2),
    dpi=150
)

density_matrix.plot(
    kind="bar",
    ax=ax,
    width=0.8
)

ax.set_xlabel("Chromosome")
ax.set_ylabel("Markers per Mb")
ax.legend(
    title="Family",
    frameon=False
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "marker_density_per_chromosome.pdf"
)

plt.close(fig)


# ============================================================
# 3. Maximum marker gap per chromosome
# ============================================================

gap_matrix = chromosome_summary.pivot(
    index="chromosome",
    columns="family_id",
    values="max_marker_gap_bp"
) / 1_000_000

fig, ax = plt.subplots(
    figsize=(7.2, 3.2),
    dpi=150
)

gap_matrix.plot(
    kind="bar",
    ax=ax,
    width=0.8
)

ax.set_xlabel("Chromosome")
ax.set_ylabel("Maximum marker gap (Mb)")
ax.legend(
    title="Family",
    frameon=False
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "maximum_marker_gap_per_chromosome.pdf"
)

plt.close(fig)


# ============================================================
# 4. Missingness per offspring
# ============================================================

plot_data = missingness.sort_values(
    ["family_id", "sample_id"]
)

fig, ax = plt.subplots(
    figsize=(6.4, 3.0),
    dpi=150
)

ax.bar(
    plot_data["sample_id"],
    plot_data["missing_rate"]
)

ax.set_xlabel("")
ax.set_ylabel("Marker missing rate")
ax.tick_params(
    axis="x",
    rotation=90
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "offspring_marker_missing_rate.pdf"
)

plt.close(fig)
