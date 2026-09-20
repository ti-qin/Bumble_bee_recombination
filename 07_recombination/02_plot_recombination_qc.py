from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PROJECT = Path("/localdata/qinti/Project/Bee")

ROOT = (
    PROJECT
    / "07_recombination/classified/"
      "threshold_10000bp"
)

FIG_DIR = (
    PROJECT
    / "07_recombination/figures/"
      "threshold_10000bp"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)

events = pd.read_csv(
    ROOT
    / "classified_recombination_events.tsv.gz",
    sep="\t",
)

co = pd.read_csv(
    ROOT
    / "CO_breakpoints.tsv.gz",
    sep="\t",
)

summary = pd.read_csv(
    ROOT
    / "summary/"
      "recombination_counts_by_offspring.tsv",
    sep="\t",
)


# ============================================================
# 1. CO number per offspring
# ============================================================

plot_data = summary.sort_values(
    [
        "family_id",
        "sample_id",
    ]
)

fig, ax = plt.subplots(
    figsize=(6.4, 3.0),
    dpi=150,
)

ax.bar(
    plot_data["sample_id"],
    plot_data["n_CO_events"],
)

ax.set_xlabel("")
ax.set_ylabel("Number of CO events")

ax.tick_params(
    axis="x",
    rotation=90,
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "CO_events_per_offspring.pdf"
)

plt.close(fig)


# ============================================================
# 2. High-confidence CO number
# ============================================================

fig, ax = plt.subplots(
    figsize=(6.4, 3.0),
    dpi=150,
)

ax.bar(
    plot_data["sample_id"],
    plot_data[
        "n_high_confidence_CO"
    ],
)

ax.set_xlabel("")
ax.set_ylabel(
    "Candidate high-confidence COs"
)

ax.tick_params(
    axis="x",
    rotation=90,
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "high_confidence_CO_per_offspring.pdf"
)

plt.close(fig)


# ============================================================
# 3. CO breakpoint resolution
# ============================================================

co_width = pd.to_numeric(
    co["breakpoint_interval_bp"],
    errors="coerce",
).dropna()

fig, ax = plt.subplots(
    figsize=(3.8, 2.8),
    dpi=150,
)

ax.hist(
    co_width,
    bins=60,
)

ax.set_xlabel(
    "CO breakpoint interval (bp)"
)

ax.set_ylabel("Number of COs")

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "CO_breakpoint_interval_distribution.pdf"
)

plt.close(fig)


# ============================================================
# 4. NCO tract spans
# ============================================================

nco = events.loc[
    events["event_class"]
    == "NCO_CANDIDATE"
].copy()

nco_span = pd.to_numeric(
    nco["core_span_bp"],
    errors="coerce",
).dropna()

if len(nco_span) > 0:
    fig, ax = plt.subplots(
        figsize=(3.8, 2.8),
        dpi=150,
    )

    ax.hist(
        nco_span,
        bins=50,
    )

    ax.axvline(
        10_000,
        linestyle="--",
        linewidth=1,
        label="10 kb threshold",
    )

    ax.set_xlabel(
        "Minimum NCO tract span (bp)"
    )

    ax.set_ylabel(
        "Number of candidates"
    )

    ax.legend(
        frameon=False,
        fontsize=7,
    )

    fig.tight_layout()

    fig.savefig(
        FIG_DIR
        / "NCO_candidate_tract_span.pdf"
    )

    plt.close(fig)


# ============================================================
# 5. Event review-status counts
# ============================================================

status_counts = (
    events["review_status"]
    .value_counts()
)

fig, ax = plt.subplots(
    figsize=(6.0, 3.0),
    dpi=150,
)

ax.bar(
    status_counts.index,
    status_counts.values,
)

ax.set_xlabel("")
ax.set_ylabel("Number of event records")

ax.tick_params(
    axis="x",
    rotation=45,
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "recombination_event_review_status.pdf"
)

plt.close(fig)

print(
    f"Figures saved under: {FIG_DIR}"
)
