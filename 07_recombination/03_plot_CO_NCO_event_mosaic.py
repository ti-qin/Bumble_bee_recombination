from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_rgb


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path(
    "/localdata/qinti/Project/Bee"
)

METADATA_FILE = (
    PROJECT
    / "00_rawdata/metadata/"
      "analysis_samples.noC3.tsv"
)

PRIMARY_BED = (
    PROJECT
    / "00_rawdata/reference/"
      "primary_chromosomes.bed"
)

PHASING_ROOT = (
    PROJECT
    / "06_phasing/families"
)

EVENT_ROOT = (
    PROJECT
    / "07_recombination/classified/"
      "majority_phase_5of5/"
      "detailed_review_10000bp"
)

EVENT_FILE = (
    EVENT_ROOT
    / "CO_NCO_events.detailed_review.tsv.gz"
)

FIG_DIR = (
    EVENT_ROOT
    / "figures/all_event_mosaics"
)

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. Plotting parameters
# ============================================================

# True:
# C1:C1_1:chromosome_1:E0000001
#
# False:
# E0000001
LABEL_FULL_EVENT_ID = False

LABEL_FONT_SIZE = 7

# Vertical distance between offspring tracks.
ROW_STEP = 3

# Multiple event-ID levels reduce overlap between nearby
# events from the same offspring.
LABEL_LANE_OFFSETS = [
    0.38,
    0.72,
    1.06,
    1.40,
    1.74,
    2.08,
]

# Two labels closer than this fraction of chromosome length
# are assigned to different label lanes where possible.
MIN_LABEL_GAP_FRACTION = 0.018

# Absolute lower bound for label separation in Mb.
MIN_LABEL_GAP_MB = 0.08

HAPLOTYPE_LINEWIDTH = 4.0

# True: chromosomes without identified events are omitted.
# False: all chromosomes are included, as in the original
# phasing QC PDF.
ONLY_CHROMOSOMES_WITH_EVENTS = False


# ============================================================
# 3. Visual encodings
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

# review_status controls only the event-ID label color.
REVIEW_STATUS_COLORS = {
    "HIGH": "#4DAF4A",
    "HIGH_AFTER_AUTO_REVIEW": "#FFB000",
    "LOW_MANUAL_REVIEW": "#D62728",
    "MISSING": "#969696",
}

UNKNOWN_STATUS_COLOR = "#969696"


# ============================================================
# 4. Helper functions
# ============================================================

def contrasting_text_color(
    background_color: str,
) -> str:
    """
    Return black or white depending on background brightness.
    """

    red, green, blue = to_rgb(
        background_color
    )

    luminance = (
        0.2126 * red
        + 0.7152 * green
        + 0.0722 * blue
    )

    if luminance > 0.58:
        return "black"

    return "white"


def format_event_id(
    event_id: str,
) -> str:
    """
    Format the event label shown on the plot.
    """

    event_id = str(event_id)

    if LABEL_FULL_EVENT_ID:
        return event_id

    return event_id.rsplit(
        ":",
        1
    )[-1]


def assign_label_lanes(
    event_subset: pd.DataFrame,
    chromosome_length_mb: float,
) -> dict:
    """
    Greedily assign nearby event labels to different
    vertical lanes.

    Returns
    -------
    dict
        DataFrame index -> label-lane index.
    """

    if event_subset.empty:
        return {}

    event_subset = (
        event_subset
        .sort_values(
            [
                "event_midpoint",
                "event_interval_start",
                "event_id",
            ]
        )
    )

    minimum_gap = max(
        chromosome_length_mb
        * MIN_LABEL_GAP_FRACTION,
        MIN_LABEL_GAP_MB,
    )

    n_lanes = len(
        LABEL_LANE_OFFSETS
    )

    last_x_by_lane = np.full(
        n_lanes,
        -np.inf,
        dtype=float,
    )

    lane_assignment = {}

    for index, event in event_subset.iterrows():
        x_position = (
            float(
                event["event_midpoint"]
            )
            / 1_000_000
        )

        available_lanes = [
            lane
            for lane in range(n_lanes)
            if (
                x_position
                - last_x_by_lane[lane]
                >= minimum_gap
            )
        ]

        if available_lanes:
            lane = available_lanes[0]

        else:
            # If every lane is occupied nearby, reuse the
            # lane whose previous label is farthest away.
            lane = int(
                np.argmin(
                    last_x_by_lane
                )
            )

        lane_assignment[index] = lane

        last_x_by_lane[lane] = (
            x_position
        )

    return lane_assignment


def draw_haplotype_blocks(
    ax,
    chromosome_data: pd.DataFrame,
    haplotype_column: str,
    y_position: float,
):
    """
    Draw consecutive H0/H1 haplotype blocks for one offspring.
    """

    positions = (
        chromosome_data[
            "position"
        ]
        .to_numpy(dtype=float)
        / 1_000_000
    )

    states = pd.to_numeric(
        chromosome_data[
            haplotype_column
        ],
        errors="coerce",
    ).to_numpy(dtype=float)

    valid = np.isfinite(
        states
    )

    positions = positions[
        valid
    ]

    states = states[
        valid
    ].astype(int)

    if len(states) == 0:
        return

    invalid_states = ~np.isin(
        states,
        [0, 1],
    )

    if invalid_states.any():
        raise ValueError(
            f"{haplotype_column}: "
            "haplotype values other than 0/1 detected: "
            f"{np.unique(states[invalid_states])}"
        )

    change_points = (
        np.where(
            states[1:]
            != states[:-1]
        )[0]
        + 1
    )

    block_starts = np.r_[
        0,
        change_points,
    ]

    block_ends = np.r_[
        change_points - 1,
        len(states) - 1,
    ]

    for start_index, end_index in zip(
        block_starts,
        block_ends,
    ):
        state = int(
            states[start_index]
        )

        x_start = float(
            positions[start_index]
        )

        x_end = float(
            positions[end_index]
        )

        if np.isclose(
            x_start,
            x_end,
        ):
            ax.scatter(
                x_start,
                y_position,
                s=6,
                color=haplotype_colors[
                    state
                ],
                zorder=2,
            )

        else:
            ax.hlines(
                y_position,
                x_start,
                x_end,
                linewidth=(
                    HAPLOTYPE_LINEWIDTH
                ),
                color=haplotype_colors[
                    state
                ],
                zorder=2,
            )


def draw_event(
    ax,
    event: pd.Series,
    y_position: float,
    label_lane: int,
):
    """
    Draw one CO or NCO event and its colored event-ID label.
    """

    event_type = str(
        event["event_type"]
    )

    review_status = str(
        event["review_status"]
    )

    event_id = format_event_id(
        event["event_id"]
    )

    x_start = (
        float(
            event[
                "event_interval_start"
            ]
        )
        / 1_000_000
    )

    x_end = (
        float(
            event[
                "event_interval_end"
            ]
        )
        / 1_000_000
    )

    x_midpoint = (
        float(
            event["event_midpoint"]
        )
        / 1_000_000
    )

    if x_end < x_start:
        x_start, x_end = (
            x_end,
            x_start,
        )

    # --------------------------------------------------------
    # Possible event/breakpoint interval
    # --------------------------------------------------------

    bracket_y = (
        y_position
        + 0.13
    )

    ax.hlines(
        bracket_y,
        x_start,
        x_end,
        linewidth=0.55,
        color="0.25",
        zorder=4,
    )

    ax.vlines(
        [
            x_start,
            x_end,
        ],
        bracket_y - 0.055,
        bracket_y + 0.055,
        linewidth=0.45,
        color="0.25",
        zorder=4,
    )

    # --------------------------------------------------------
    # Event-type encoding
    # --------------------------------------------------------

    if event_type == "CO":
        # CO = vertical switching marker.
        ax.vlines(
            x_midpoint,
            y_position - 0.23,
            y_position + 0.23,
            linewidth=0.9,
            color="black",
            zorder=5,
        )

        ax.scatter(
            x_midpoint,
            y_position + 0.23,
            marker="v",
            s=11,
            facecolor="black",
            edgecolor="black",
            linewidth=0.3,
            zorder=6,
        )

    elif event_type == "NCO":
        # NCO = diamond at event midpoint.
        ax.scatter(
            x_midpoint,
            y_position,
            marker="D",
            s=16,
            facecolor="white",
            edgecolor="black",
            linewidth=0.7,
            zorder=6,
        )

    else:
        ax.scatter(
            x_midpoint,
            y_position,
            marker="o",
            s=14,
            facecolor="white",
            edgecolor="black",
            linewidth=0.7,
            zorder=6,
        )

    # --------------------------------------------------------
    # Core tract, where available
    # --------------------------------------------------------

    core_start = pd.to_numeric(
        pd.Series(
            [event.get("core_start", np.nan)]
        ),
        errors="coerce",
    ).iloc[0]

    core_end = pd.to_numeric(
        pd.Series(
            [event.get("core_end", np.nan)]
        ),
        errors="coerce",
    ).iloc[0]

    if (
        pd.notna(core_start)
        and
        pd.notna(core_end)
    ):
        core_x_start = (
            float(core_start)
            / 1_000_000
        )

        core_x_end = (
            float(core_end)
            / 1_000_000
        )

        if core_x_end < core_x_start:
            core_x_start, core_x_end = (
                core_x_end,
                core_x_start,
            )

        # The short tract may be invisible at chromosome
        # scale, so both a segment and a midpoint marker
        # are retained.
        ax.hlines(
            y_position - 0.15,
            core_x_start,
            core_x_end,
            linewidth=2.2,
            color="black",
            zorder=5,
        )

    # --------------------------------------------------------
    # Event-ID label colored by review_status
    # --------------------------------------------------------

    status_color = (
        REVIEW_STATUS_COLORS.get(
            review_status,
            UNKNOWN_STATUS_COLOR,
        )
    )

    label_text_color = (
        contrasting_text_color(
            status_color
        )
    )

    label_y = (
        y_position
        + LABEL_LANE_OFFSETS[
            label_lane
        ]
    )

    # Connector from event marker to ID label.
    ax.plot(
        [
            x_midpoint,
            x_midpoint,
        ],
        [
            y_position + 0.24,
            label_y - 0.035,
        ],
        linewidth=0.35,
        color="0.35",
        zorder=7,
    )

    ax.text(
        x_midpoint,
        label_y,
        event_id,
        rotation=30,
        rotation_mode="anchor",
        ha="left",
        va="center",
        fontsize=LABEL_FONT_SIZE,
        color=label_text_color,
        clip_on=False,
        zorder=8,
        bbox={
            "boxstyle": (
                "round,pad=0.16"
            ),
            "facecolor": (
                status_color
            ),
            "edgecolor": "black",
            "linewidth": 0.25,
            "alpha": 0.94,
        },
    )


# ============================================================
# 5. Read metadata, chromosomes and event table
# ============================================================

metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str,
)

primary_bed = pd.read_csv(
    PRIMARY_BED,
    sep="\t",
    header=None,
    names=[
        "chromosome",
        "start",
        "end",
    ],
)

primary_bed["chromosome"] = (
    primary_bed[
        "chromosome"
    ].astype(str)
)

primary_bed["start"] = pd.to_numeric(
    primary_bed["start"],
    errors="raise",
).astype(int)

primary_bed["end"] = pd.to_numeric(
    primary_bed["end"],
    errors="raise",
).astype(int)

chromosome_order = (
    primary_bed[
        "chromosome"
    ].tolist()
)

chromosome_lengths = dict(
    zip(
        primary_bed[
            "chromosome"
        ],
        (
            primary_bed["end"]
            - primary_bed["start"]
        ),
    )
)

events = pd.read_csv(
    EVENT_FILE,
    sep="\t",
    dtype={
        "family_id": str,
        "sample_id": str,
        "chromosome": str,
        "event_id": str,
        "event_type": str,
        "event_context": str,
        "review_status": str,
    },
)


# ============================================================
# 6. Validate the event table
# ============================================================

required_event_columns = [
    "family_id",
    "sample_id",
    "chromosome",
    "event_id",
    "event_type",
    "event_context",
    "review_status",
    "event_interval_start",
    "event_interval_end",
    "event_midpoint",
    "core_start",
    "core_end",
]

missing_columns = [
    column
    for column in required_event_columns
    if column not in events.columns
]

if missing_columns:
    raise KeyError(
        "Missing columns in event table: "
        + ", ".join(missing_columns)
    )

if events["event_id"].duplicated().any():
    duplicated_ids = (
        events.loc[
            events[
                "event_id"
            ].duplicated(
                keep=False
            ),
            "event_id",
        ]
        .tolist()
    )

    raise ValueError(
        "Duplicated event_id values detected:\n"
        + "\n".join(
            duplicated_ids[:20]
        )
    )

numeric_event_columns = [
    "event_interval_start",
    "event_interval_end",
    "event_midpoint",
    "core_start",
    "core_end",
]

for column in numeric_event_columns:
    events[column] = pd.to_numeric(
        events[column],
        errors=(
            "coerce"
            if column in [
                "core_start",
                "core_end",
            ]
            else "raise"
        ),
    )

events["review_status"] = (
    events["review_status"]
    .fillna("MISSING")
    .astype(str)
)

known_statuses = set(
    REVIEW_STATUS_COLORS
)

observed_statuses = set(
    events["review_status"].unique()
)

unknown_statuses = sorted(
    observed_statuses
    - known_statuses
)

if unknown_statuses:
    warnings.warn(
        "The following review_status values are not "
        "defined in REVIEW_STATUS_COLORS and will be "
        "shown in grey: "
        + ", ".join(unknown_statuses)
    )


# ============================================================
# 7. Family and sample validation
# ============================================================

families = sorted(
    metadata[
        "family_id"
    ].dropna().unique()
)

metadata_samples = set(
    metadata["sample_id"]
)

unknown_event_samples = sorted(
    set(events["sample_id"])
    - metadata_samples
)

if unknown_event_samples:
    raise ValueError(
        "Event samples absent from metadata:\n"
        + "\n".join(
            unknown_event_samples
        )
    )


# ============================================================
# 8. Legend handles
# ============================================================

haplotype_handles = [
    Line2D(
        [0],
        [0],
        color=haplotype_colors[0],
        linewidth=4,
        label="H0",
    ),
    Line2D(
        [0],
        [0],
        color=haplotype_colors[1],
        linewidth=4,
        label="H1",
    ),
]

event_type_handles = [
    Line2D(
        [0],
        [0],
        marker="|",
        linestyle="None",
        color="black",
        markeredgewidth=1.2,
        markersize=10,
        label="CO",
    ),
    Line2D(
        [0],
        [0],
        marker="D",
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="black",
        markersize=4.5,
        label="NCO",
    ),
]

status_order = [
    "HIGH",
    "HIGH_AFTER_AUTO_REVIEW",
    "LOW_MANUAL_REVIEW",
]

status_order.extend(
    status
    for status in sorted(
        observed_statuses
    )
    if status not in status_order
)

status_handles = []

for status in status_order:
    if status not in observed_statuses:
        continue

    status_handles.append(
        Patch(
            facecolor=(
                REVIEW_STATUS_COLORS.get(
                    status,
                    UNKNOWN_STATUS_COLOR,
                )
            ),
            edgecolor="black",
            linewidth=0.3,
            label=status,
        )
    )


# ============================================================
# 9. Plot all family/chromosome mosaics
# ============================================================

for family in families:
    phased_file = (
        PHASING_ROOT
        / family
        / f"{family}.phased_markers.tsv.gz"
    )

    if not phased_file.exists():
        warnings.warn(
            f"Missing phased-marker file: "
            f"{phased_file}"
        )
        continue

    phased = pd.read_csv(
        phased_file,
        sep="\t",
    )

    phased["chromosome"] = (
        phased[
            "chromosome"
        ].astype(str)
    )

    phased["position"] = pd.to_numeric(
        phased["position"],
        errors="raise",
    ).astype(int)

    offspring = (
        metadata.loc[
            (
                metadata["family_id"]
                == family
            )
            &
            (
                metadata["role"]
                == "offspring"
            ),
            "sample_id",
        ]
        .sort_values()
        .tolist()
    )

    family_events = events.loc[
        events["family_id"]
        == family
    ].copy()

    output_pdf = (
        FIG_DIR
        / (
            f"{family}."
            "all_CO_NCO_events."
            "review_status.pdf"
        )
    )

    with PdfPages(
        output_pdf
    ) as pdf:
        for chromosome in chromosome_order:
            chromosome_data = (
                phased.loc[
                    phased[
                        "chromosome"
                    ]
                    == chromosome
                ]
                .sort_values(
                    "position"
                )
                .reset_index(
                    drop=True
                )
            )

            if chromosome_data.empty:
                continue

            chromosome_events = (
                family_events.loc[
                    family_events[
                        "chromosome"
                    ]
                    == chromosome
                ]
                .copy()
            )

            if (
                ONLY_CHROMOSOMES_WITH_EVENTS
                and chromosome_events.empty
            ):
                continue

            chromosome_length_bp = (
                chromosome_lengths.get(
                    chromosome,
                    chromosome_data[
                        "position"
                    ].max(),
                )
            )

            chromosome_length_mb = (
                chromosome_length_bp
                / 1_000_000
            )

            figure_height = max(
                4.6,
                1.1
                + len(offspring)
                * 1.15,
            )

            fig, ax = plt.subplots(
                figsize=(
                    10.5,
                    figure_height,
                ),
                dpi=150,
            )

            row_positions = {
                child: (
                    child_index
                    * ROW_STEP
                )
                for child_index, child
                in enumerate(offspring)
            }

            # ------------------------------------------------
            # Haplotype mosaic
            # ------------------------------------------------

            for child in offspring:
                haplotype_column = (
                    f"{child}_haplotype"
                )

                if (
                    haplotype_column
                    not in chromosome_data.columns
                ):
                    warnings.warn(
                        f"Missing column "
                        f"{haplotype_column} in "
                        f"{family} {chromosome}"
                    )
                    continue

                draw_haplotype_blocks(
                    ax=ax,
                    chromosome_data=(
                        chromosome_data
                    ),
                    haplotype_column=(
                        haplotype_column
                    ),
                    y_position=(
                        row_positions[child]
                    ),
                )

            # ------------------------------------------------
            # CO/NCO event overlays
            # ------------------------------------------------

            for child in offspring:
                child_events = (
                    chromosome_events.loc[
                        chromosome_events[
                            "sample_id"
                        ]
                        == child
                    ]
                    .sort_values(
                        [
                            "event_midpoint",
                            "event_interval_start",
                        ]
                    )
                )

                lane_assignment = (
                    assign_label_lanes(
                        event_subset=(
                            child_events
                        ),
                        chromosome_length_mb=(
                            chromosome_length_mb
                        ),
                    )
                )

                for event_index, event in (
                    child_events.iterrows()
                ):
                    draw_event(
                        ax=ax,
                        event=event,
                        y_position=(
                            row_positions[
                                child
                            ]
                        ),
                        label_lane=(
                            lane_assignment[
                                event_index
                            ]
                        ),
                    )

            # ------------------------------------------------
            # Axes
            # ------------------------------------------------

            y_tick_positions = [
                row_positions[child]
                for child in offspring
            ]

            ax.set_yticks(
                y_tick_positions
            )

            ax.set_yticklabels(
                offspring
            )

            ax.set_xlim(
                0,
                chromosome_length_mb,
            )

            final_row_y = (
                max(y_tick_positions)
                if y_tick_positions
                else 0
            )

            ax.set_ylim(
                -0.65,
                final_row_y
                + max(
                    LABEL_LANE_OFFSETS
                )
                + 0.65,
            )

            ax.set_xlabel(
                "Chromosomal position (Mb)"
            )

            ax.set_ylabel(
                "Offspring"
            )

            n_co = int(
                (
                    chromosome_events[
                        "event_type"
                    ]
                    == "CO"
                ).sum()
            )

            n_nco = int(
                (
                    chromosome_events[
                        "event_type"
                    ]
                    == "NCO"
                ).sum()
            )

            ax.set_title(
                f"{family} — {chromosome}  "
                f"(CO={n_co}, NCO={n_nco})"
            )

            ax.grid(
                axis="x",
                linewidth=0.35,
                alpha=0.25,
            )

            ax.tick_params(
                axis="both",
                labelsize=7,
            )

            # ------------------------------------------------
            # Legends
            # ------------------------------------------------

            event_legend = ax.legend(
                handles=(
                    haplotype_handles
                    + event_type_handles
                ),
                loc="upper left",
                bbox_to_anchor=(
                    0,
                    -0.14,
                ),
                ncol=4,
                frameon=False,
                fontsize=6.5,
                handlelength=1.8,
                columnspacing=1.2,
            )

            ax.add_artist(
                event_legend
            )

            if status_handles:
                ax.legend(
                    handles=status_handles,
                    title="review_status",
                    loc="upper right",
                    bbox_to_anchor=(
                        1,
                        -0.13,
                    ),
                    ncol=min(
                        3,
                        len(status_handles),
                    ),
                    frameon=False,
                    fontsize=6.2,
                    title_fontsize=6.5,
                    handlelength=1.1,
                    columnspacing=1.0,
                )

            fig.subplots_adjust(
                left=0.09,
                right=0.985,
                top=0.91,
                bottom=0.22,
            )

            pdf.savefig(
                fig,
                bbox_inches="tight",
            )

            plt.close(
                fig
            )

    print(
        f"Saved: {output_pdf}"
    )


# ============================================================
# 10. Summary
# ============================================================

print(
    "\nAll-event haplotype mosaics completed."
)

print(
    f"\nFigures:\n{FIG_DIR}"
)

print(
    "\nreview_status counts:"
)

print(
    events[
        "review_status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)
