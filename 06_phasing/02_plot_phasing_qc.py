from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


PROJECT = Path("/localdata/qinti/Project/Bee")

PHASING_ROOT = (
    PROJECT / "06_phasing/families"
)

SUMMARY_DIR = (
    PROJECT / "06_phasing/summary"
)

FIG_DIR = (
    PROJECT / "06_phasing/figures"
)

METADATA = (
    PROJECT
    / "00_rawdata/metadata/"
      "analysis_samples.noC3.tsv"
)

PRIMARY_BED = (
    PROJECT
    / "00_rawdata/reference/"
      "primary_chromosomes.bed"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str
)

primary_bed = pd.read_csv(
    PRIMARY_BED,
    sep="\t",
    header=None,
    names=[
        "chromosome",
        "start",
        "end",
    ]
)

chromosome_order = (
    primary_bed["chromosome"]
    .astype(str)
    .tolist()
)

families = sorted(
    metadata["family_id"].unique()
)

links = pd.read_csv(
    SUMMARY_DIR
    / "all_families_adjacent_phase_links.tsv.gz",
    sep="\t"
)

switch_summary = pd.read_csv(
    SUMMARY_DIR
    / "all_families_raw_switch_summary.tsv",
    sep="\t"
)

chromosome_summary = pd.read_csv(
    SUMMARY_DIR
    / "phasing_chromosome_summary.tsv",
    sep="\t"
)


# ============================================================
# 1. Phase-support histogram
# ============================================================

values = links[
    "phase_support"
].dropna()

fig, ax = plt.subplots(
    figsize=(3.8, 2.8),
    dpi=150
)

ax.hist(
    values,
    bins=np.linspace(
        0.5,
        1.0,
        11
    )
)

ax.axvline(
    0.8,
    linestyle="--",
    linewidth=1,
    label="4/5 support"
)

ax.set_xlabel(
    "Adjacent-marker phase support"
)

ax.set_ylabel(
    "Number of marker links"
)

ax.legend(
    frameon=False,
    fontsize=7
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "adjacent_phase_support_histogram.pdf"
)

plt.close(fig)


# ============================================================
# 2. Weak links per chromosome
# ============================================================

weak_matrix = (
    chromosome_summary
    .pivot(
        index="chromosome",
        columns="family_id",
        values="n_weak_majority_links"
    )
    .reindex(chromosome_order)
    .fillna(0)
)

fig, ax = plt.subplots(
    figsize=(7.2, 3.2),
    dpi=150
)

weak_matrix.plot(
    kind="bar",
    ax=ax,
    width=0.8
)

ax.set_xlabel("Chromosome")
ax.set_ylabel(
    "Number of weak 3:2 phase links"
)

ax.legend(
    title="Family",
    frameon=False
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "weak_phase_links_per_chromosome.pdf"
)

plt.close(fig)


# ============================================================
# 3. Raw switches per offspring
# ============================================================

offspring_total = (
    switch_summary
    .groupby(
        [
            "family_id",
            "sample_id",
        ],
        as_index=False
    )
    .agg(
        n_raw_haplotype_switches=(
            "n_raw_haplotype_switches",
            "sum"
        ),
        n_comparable_marker_links=(
            "n_comparable_marker_links",
            "sum"
        ),
    )
)

offspring_total[
    "raw_switch_fraction"
] = (
    offspring_total[
        "n_raw_haplotype_switches"
    ]
    / offspring_total[
        "n_comparable_marker_links"
    ]
)

offspring_total = (
    offspring_total
    .sort_values(
        [
            "family_id",
            "sample_id",
        ]
    )
)

fig, ax = plt.subplots(
    figsize=(6.4, 3.0),
    dpi=150
)

ax.bar(
    offspring_total["sample_id"],
    offspring_total[
        "n_raw_haplotype_switches"
    ]
)

ax.set_xlabel("")
ax.set_ylabel(
    "Raw H0/H1 switches"
)

ax.tick_params(
    axis="x",
    rotation=90
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "raw_haplotype_switches_per_offspring.pdf"
)

plt.close(fig)


# ============================================================
# 4. Multi-page haplotype mosaic for each family
# ============================================================

default_colors = (
    plt.rcParams[
        "axes.prop_cycle"
    ]
    .by_key()["color"]
)

haplotype_colors = {
    0: default_colors[0],
    1: default_colors[1],
}

for family in families:
    phased_file = (
        PHASING_ROOT
        / family
        / f"{family}.phased_markers.tsv.gz"
    )

    phased = pd.read_csv(
        phased_file,
        sep="\t"
    )

    offspring = (
        metadata.loc[
            (
                metadata["family_id"]
                == family
            )
            & (
                metadata["role"]
                == "offspring"
            ),
            "sample_id"
        ]
        .sort_values()
        .tolist()
    )

    output_pdf = (
        FIG_DIR
        / f"{family}.phased_haplotype_mosaic.pdf"
    )

    with PdfPages(output_pdf) as pdf:
        for chromosome in chromosome_order:
            chromosome_data = phased.loc[
                phased["chromosome"]
                == chromosome
            ].sort_values("position")

            if chromosome_data.empty:
                continue

            fig, ax = plt.subplots(
                figsize=(8.0, 2.8),
                dpi=150
            )

            for child_index, child in enumerate(
                offspring
            ):
                state_column = (
                    f"{child}_haplotype"
                )

                positions = (
                    chromosome_data[
                        "position"
                    ].to_numpy(
                        dtype=float
                    )
                    / 1_000_000
                )

                states = pd.to_numeric(
                    chromosome_data[
                        state_column
                    ],
                    errors="coerce"
                ).to_numpy(
                    dtype=float
                )

                valid = np.isfinite(states)

                positions = positions[valid]
                states = states[valid].astype(int)

                if len(states) == 0:
                    continue

                change_points = (
                    np.where(
                        states[1:]
                        != states[:-1]
                    )[0]
                    + 1
                )

                block_starts = np.r_[
                    0,
                    change_points
                ]

                block_ends = np.r_[
                    change_points - 1,
                    len(states) - 1
                ]

                for start, end in zip(
                    block_starts,
                    block_ends
                ):
                    state = int(
                        states[start]
                    )

                    x_start = positions[start]
                    x_end = positions[end]

                    if x_start == x_end:
                        ax.scatter(
                            x_start,
                            child_index,
                            s=5,
                            color=(
                                haplotype_colors[
                                    state
                                ]
                            ),
                        )
                    else:
                        ax.hlines(
                            child_index,
                            x_start,
                            x_end,
                            linewidth=4,
                            color=(
                                haplotype_colors[
                                    state
                                ]
                            ),
                        )

            ax.set_yticks(
                np.arange(
                    len(offspring)
                )
            )

            ax.set_yticklabels(
                offspring
            )

            ax.set_xlabel(
                "Chromosomal position (Mb)"
            )

            ax.set_ylabel("Offspring")

            ax.set_title(
                f"{family} — {chromosome}"
            )

            ax.set_ylim(
                -0.75,
                len(offspring) - 0.25
            )

            fig.tight_layout()

            pdf.savefig(fig)

            plt.close(fig)

print(
    f"Phasing QC figures saved under: "
    f"{FIG_DIR}"
)
