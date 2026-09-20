from pathlib import Path
from collections import defaultdict
from bisect import bisect_right
import argparse

import numpy as np
import pandas as pd


# ============================================================
# 1. Arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Identify and classify CO, NCO and "
        "CO-associated gene-conversion candidates "
        "from phased haploid offspring markers."
    )
)

parser.add_argument(
    "--tract-threshold-bp",
    type=int,
    default=10_000,
    help=(
        "Maximum internal tract span classified "
        "as a gene-conversion candidate."
    ),
)

parser.add_argument(
    "--large-breakpoint-bp",
    type=int,
    default=100_000,
    help=(
        "Breakpoint intervals wider than this "
        "are flagged as low resolution."
    ),
)

parser.add_argument(
    "--shared-tolerance-bp",
    type=int,
    default=0,
    help=(
        "Distance tolerance for clustering "
        "events across offspring."
    ),
)

args = parser.parse_args()

TRACT_THRESHOLD_BP = (
    args.tract_threshold_bp
)

LARGE_BREAKPOINT_BP = (
    args.large_breakpoint_bp
)

SHARED_TOLERANCE_BP = (
    args.shared_tolerance_bp
)


# ============================================================
# 2. Paths
# ============================================================

PROJECT = Path(
    "/localdata/qinti/Project/Bee"
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

PHASING_ROOT = (
    PROJECT
    / "06_phasing/families"
)

GAP_BED = (
    PROJECT
    / "07_recombination/reference_masks/"
      "reference_N_gaps.bed"
)

OUTPUT_ROOT = (
    PROJECT
    / "07_recombination/classified"
    / f"threshold_{TRACT_THRESHOLD_BP}bp"
)

SUMMARY_DIR = (
    OUTPUT_ROOT / "summary"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 3. Metadata
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
)

families = sorted(
    metadata["family_id"].unique()
)

family_offspring = {}

for family in families:
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
            "sample_id",
        ]
        .sort_values()
        .tolist()
    )

    family_offspring[family] = offspring


# ============================================================
# 4. Chromosome information
# ============================================================

primary = pd.read_csv(
    PRIMARY_BED,
    sep="\t",
    header=None,
    names=[
        "chromosome",
        "start",
        "end",
    ],
)

primary["chromosome"] = (
    primary["chromosome"].astype(str)
)

chromosome_order = (
    primary["chromosome"].tolist()
)

chromosome_rank = {
    chromosome: index
    for index, chromosome
    in enumerate(chromosome_order)
}

chromosome_lengths = dict(
    zip(
        primary["chromosome"],
        primary["end"]
        - primary["start"],
    )
)


# ============================================================
# 5. Reference gap lookup
# ============================================================

gap_intervals = defaultdict(list)
gap_starts = {}

if GAP_BED.exists() and GAP_BED.stat().st_size > 0:
    gaps = pd.read_csv(
        GAP_BED,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "start",
            "end",
        ],
    )

    gaps["chromosome"] = (
        gaps["chromosome"].astype(str)
    )

    for chromosome, subset in gaps.groupby(
        "chromosome"
    ):
        intervals = list(
            zip(
                subset["start"].astype(int),
                subset["end"].astype(int),
            )
        )

        intervals.sort()

        gap_intervals[chromosome] = (
            intervals
        )

        gap_starts[chromosome] = [
            interval[0]
            for interval in intervals
        ]


def overlaps_reference_gap(
    chromosome: str,
    start_1based: int,
    end_1based: int,
) -> bool:
    """
    Test overlap between a one-based inclusive
    interval and zero-based half-open gap BED.
    """

    intervals = gap_intervals.get(
        chromosome,
        [],
    )

    if not intervals:
        return False

    start0 = max(
        0,
        int(start_1based) - 1,
    )

    end0 = int(end_1based)

    starts = gap_starts[chromosome]

    index = (
        bisect_right(
            starts,
            end0 - 1,
        )
        - 1
    )

    for candidate in [
        index - 1,
        index,
        index + 1,
    ]:
        if (
            candidate < 0
            or candidate >= len(intervals)
        ):
            continue

        gap_start, gap_end = intervals[
            candidate
        ]

        if (
            gap_start < end0
            and gap_end > start0
        ):
            return True

    return False


# ============================================================
# 6. Build haplotype blocks
# ============================================================

def build_blocks(
    chromosome_data: pd.DataFrame,
    state_column: str,
):
    states = pd.to_numeric(
        chromosome_data[state_column],
        errors="coerce",
    ).to_numpy(dtype=float)

    called_mask = np.isfinite(states)

    called = chromosome_data.loc[
        called_mask
    ].copy()

    called["_row_index"] = np.where(
        called_mask
    )[0]

    called["_state"] = (
        states[called_mask]
        .astype(np.int8)
    )

    called = called.reset_index(
        drop=True
    )

    if called.empty:
        return []

    blocks = []
    start_index = 0

    for index in range(
        1,
        len(called),
    ):
        if (
            called.loc[index, "_state"]
            != called.loc[
                index - 1,
                "_state"
            ]
        ):
            block_data = called.iloc[
                start_index:index
            ]

            blocks.append(
                make_block(
                    block_data,
                    len(blocks),
                )
            )

            start_index = index

    block_data = called.iloc[
        start_index:
    ]

    blocks.append(
        make_block(
            block_data,
            len(blocks),
        )
    )

    return blocks


def make_block(
    block_data: pd.DataFrame,
    block_index: int,
):
    start_position = int(
        block_data["position"].iloc[0]
    )

    end_position = int(
        block_data["position"].iloc[-1]
    )

    first_row_index = int(
        block_data["_row_index"].iloc[0]
    )

    last_row_index = int(
        block_data["_row_index"].iloc[-1]
    )

    n_markers = len(block_data)

    return {
        "block_index": block_index,
        "state": int(
            block_data["_state"].iloc[0]
        ),
        "start_position": start_position,
        "end_position": end_position,
        "span_bp": (
            end_position
            - start_position
        ),
        "n_markers": n_markers,
        "first_row_index": (
            first_row_index
        ),
        "last_row_index": (
            last_row_index
        ),
        "n_missing_markers_inside": (
            last_row_index
            - first_row_index
            + 1
            - n_markers
        ),
    }


# ============================================================
# 7. Adjacent block boundaries
# ============================================================

BAD_PHASE_CATEGORIES = {
    "weak_majority",
    "tie",
    "unresolved",
}


def build_boundaries(
    chromosome_data: pd.DataFrame,
    blocks: list[dict],
    chromosome: str,
):
    boundaries = []

    for boundary_index in range(
        len(blocks) - 1
    ):
        left_block = blocks[
            boundary_index
        ]

        right_block = blocks[
            boundary_index + 1
        ]

        phase_slice = chromosome_data.iloc[
            (
                left_block[
                    "last_row_index"
                ]
                + 1
            ):
            (
                right_block[
                    "first_row_index"
                ]
                + 1
            )
        ]

        supports = pd.to_numeric(
            phase_slice[
                "phase_support_from_previous"
            ],
            errors="coerce",
        ).dropna()

        categories = (
            phase_slice[
                "phase_category_from_previous"
            ]
            .dropna()
            .astype(str)
            .tolist()
        )

        n_bad_links = sum(
            category
            in BAD_PHASE_CATEGORIES
            for category in categories
        )

        left_position = int(
            left_block["end_position"]
        )

        right_position = int(
            right_block["start_position"]
        )

        boundaries.append(
            {
                "boundary_index": (
                    boundary_index
                ),
                "left_block_index": (
                    boundary_index
                ),
                "right_block_index": (
                    boundary_index + 1
                ),
                "left_state": (
                    left_block["state"]
                ),
                "right_state": (
                    right_block["state"]
                ),
                "left_marker_position": (
                    left_position
                ),
                "right_marker_position": (
                    right_position
                ),
                "breakpoint_interval_bp": (
                    right_position
                    - left_position
                ),
                "breakpoint_midpoint": (
                    (
                        left_position
                        + right_position
                    )
                    / 2
                ),
                "minimum_phase_support": (
                    supports.min()
                    if len(supports) > 0
                    else np.nan
                ),
                "median_phase_support": (
                    supports.median()
                    if len(supports) > 0
                    else np.nan
                ),
                "n_phase_links_crossed": (
                    len(phase_slice)
                ),
                "n_weak_or_unresolved_links": (
                    n_bad_links
                ),
                "phase_issue": (
                    n_bad_links > 0
                ),
                "phase_categories_crossed": (
                    ",".join(categories)
                ),
                "overlaps_reference_gap": (
                    overlaps_reference_gap(
                        chromosome,
                        left_position,
                        right_position,
                    )
                ),
            }
        )

    return boundaries


def summarize_boundaries(
    boundaries: list[dict],
    boundary_indices: list[int],
):
    selected = [
        boundaries[index]
        for index in boundary_indices
    ]

    supports = [
        boundary[
            "minimum_phase_support"
        ]
        for boundary in selected
        if np.isfinite(
            boundary[
                "minimum_phase_support"
            ]
        )
    ]

    return {
        "minimum_phase_support": (
            min(supports)
            if supports
            else np.nan
        ),
        "n_weak_or_unresolved_links": (
            sum(
                boundary[
                    "n_weak_or_unresolved_links"
                ]
                for boundary in selected
            )
        ),
        "phase_issue": any(
            boundary["phase_issue"]
            for boundary in selected
        ),
        "overlaps_reference_gap": any(
            boundary[
                "overlaps_reference_gap"
            ]
            for boundary in selected
        ),
        "maximum_breakpoint_interval_bp": (
            max(
                boundary[
                    "breakpoint_interval_bp"
                ]
                for boundary in selected
            )
        ),
        "boundary_indices": ",".join(
            str(index)
            for index in boundary_indices
        ),
    }


# ============================================================
# 8. Output containers
# ============================================================

block_rows = []
boundary_rows = []
event_rows = []
co_rows = []

event_serial = 0
co_serial = 0
double_co_serial = 0


# ============================================================
# 9. Event classification
# ============================================================

for family in families:
    phased_file = (
        PHASING_ROOT
        / family
        / f"{family}.phased_markers.tsv.gz"
    )

    phased = pd.read_csv(
        phased_file,
        sep="\t",
    )

    phased["chromosome"] = (
        phased["chromosome"].astype(str)
    )

    phased["position"] = pd.to_numeric(
        phased["position"],
        errors="raise",
    ).astype(int)

    offspring = family_offspring[
        family
    ]

    for chromosome in chromosome_order:
        chromosome_data = (
            phased.loc[
                phased["chromosome"]
                == chromosome
            ]
            .sort_values("position")
            .reset_index(drop=True)
            .copy()
        )

        if chromosome_data.empty:
            continue

        for child in offspring:
            state_column = (
                f"{child}_haplotype"
            )

            blocks = build_blocks(
                chromosome_data,
                state_column,
            )

            if len(blocks) == 0:
                continue

            boundaries = build_boundaries(
                chromosome_data,
                blocks,
                chromosome,
            )

            for block in blocks:
                block_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        **block,
                    }
                )

            for boundary in boundaries:
                boundary_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        **boundary,
                    }
                )

            if len(blocks) < 2:
                continue

            consumed_boundaries = set()

            # ------------------------------------------------
            # 9A. CO-associated gene-conversion pattern:
            #
            # A | B-short | A-short | B
            # ------------------------------------------------

            for block_index in range(
                1,
                len(blocks) - 2,
            ):
                middle_1 = blocks[
                    block_index
                ]

                middle_2 = blocks[
                    block_index + 1
                ]

                boundary_indices = [
                    block_index - 1,
                    block_index,
                    block_index + 1,
                ]

                if any(
                    index in consumed_boundaries
                    for index in boundary_indices
                ):
                    continue

                if not (
                    middle_1["span_bp"]
                    <= TRACT_THRESHOLD_BP
                    and
                    middle_2["span_bp"]
                    <= TRACT_THRESHOLD_BP
                ):
                    continue

                left_flank = blocks[
                    block_index - 1
                ]

                right_flank = blocks[
                    block_index + 2
                ]

                # Overall haplotype must change
                # from the left flank to the right flank.
                if (
                    left_flank["state"]
                    == right_flank["state"]
                ):
                    continue

                event_serial += 1
                co_serial += 1

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_serial:07d}"
                )

                co_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:CO"
                    f"{co_serial:07d}"
                )

                qc = summarize_boundaries(
                    boundaries,
                    boundary_indices,
                )

                event_start = int(
                    left_flank[
                        "end_position"
                    ]
                )

                event_end = int(
                    right_flank[
                        "start_position"
                    ]
                )

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_class": (
                            "CO_ASSOCIATED_GC_CANDIDATE"
                        ),
                        "n_CO_events": 1,
                        "n_NCO_events": 0,
                        "n_GC_events": 1,
                        "state_before": (
                            left_flank["state"]
                        ),
                        "state_after": (
                            right_flank["state"]
                        ),
                        "event_interval_start": (
                            event_start
                        ),
                        "event_interval_end": (
                            event_end
                        ),
                        "event_interval_bp": (
                            event_end
                            - event_start
                        ),
                        "core_start": (
                            middle_1[
                                "start_position"
                            ]
                        ),
                        "core_end": (
                            middle_2[
                                "end_position"
                            ]
                        ),
                        "core_span_bp": (
                            middle_2[
                                "end_position"
                            ]
                            - middle_1[
                                "start_position"
                            ]
                        ),
                        "component_block_spans_bp": (
                            f"{middle_1['span_bp']},"
                            f"{middle_2['span_bp']}"
                        ),
                        "n_core_markers": (
                            middle_1["n_markers"]
                            + middle_2["n_markers"]
                        ),
                        "single_marker_tract": (
                            middle_1[
                                "n_markers"
                            ] == 1
                            or
                            middle_2[
                                "n_markers"
                            ] == 1
                        ),
                        "double_co_group_id": "",
                        **qc,
                    }
                )

                co_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "co_id": co_id,
                        "event_id": event_id,
                        "co_context": (
                            "CO_ASSOCIATED_GC_UNRESOLVED"
                        ),
                        "left_marker_position": (
                            event_start
                        ),
                        "right_marker_position": (
                            event_end
                        ),
                        "breakpoint_interval_bp": (
                            event_end
                            - event_start
                        ),
                        "breakpoint_midpoint": (
                            (
                                event_start
                                + event_end
                            )
                            / 2
                        ),
                        "double_co_group_id": "",
                        **qc,
                    }
                )

                consumed_boundaries.update(
                    boundary_indices
                )

            # ------------------------------------------------
            # 9B. Internal short tract:
            #
            # A | B-short | A
            # ------------------------------------------------

            for block_index in range(
                1,
                len(blocks) - 1,
            ):
                tract = blocks[
                    block_index
                ]

                boundary_indices = [
                    block_index - 1,
                    block_index,
                ]

                if any(
                    index in consumed_boundaries
                    for index in boundary_indices
                ):
                    continue

                if (
                    tract["span_bp"]
                    > TRACT_THRESHOLD_BP
                ):
                    continue

                left_flank = blocks[
                    block_index - 1
                ]

                right_flank = blocks[
                    block_index + 1
                ]

                if (
                    left_flank["state"]
                    != right_flank["state"]
                ):
                    continue

                event_serial += 1

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_serial:07d}"
                )

                qc = summarize_boundaries(
                    boundaries,
                    boundary_indices,
                )

                event_start = int(
                    left_flank[
                        "end_position"
                    ]
                )

                event_end = int(
                    right_flank[
                        "start_position"
                    ]
                )

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_class": (
                            "NCO_CANDIDATE"
                        ),
                        "n_CO_events": 0,
                        "n_NCO_events": 1,
                        "n_GC_events": 1,
                        "state_before": (
                            left_flank["state"]
                        ),
                        "state_after": (
                            right_flank["state"]
                        ),
                        "event_interval_start": (
                            event_start
                        ),
                        "event_interval_end": (
                            event_end
                        ),
                        "event_interval_bp": (
                            event_end
                            - event_start
                        ),
                        "core_start": (
                            tract[
                                "start_position"
                            ]
                        ),
                        "core_end": (
                            tract[
                                "end_position"
                            ]
                        ),
                        "core_span_bp": (
                            tract["span_bp"]
                        ),
                        "component_block_spans_bp": (
                            str(
                                tract["span_bp"]
                            )
                        ),
                        "n_core_markers": (
                            tract["n_markers"]
                        ),
                        "single_marker_tract": (
                            tract[
                                "n_markers"
                            ] == 1
                        ),
                        "double_co_group_id": "",
                        **qc,
                    }
                )

                consumed_boundaries.update(
                    boundary_indices
                )

            # ------------------------------------------------
            # 9C. Identify long internal blocks that
            #     produce paired CO breakpoints.
            # ------------------------------------------------

            double_co_boundary_map = {}

            for block_index in range(
                1,
                len(blocks) - 1,
            ):
                tract = blocks[
                    block_index
                ]

                left_boundary = (
                    block_index - 1
                )

                right_boundary = (
                    block_index
                )

                if (
                    tract["span_bp"]
                    <= TRACT_THRESHOLD_BP
                ):
                    continue

                if (
                    left_boundary
                    in consumed_boundaries
                    or
                    right_boundary
                    in consumed_boundaries
                ):
                    continue

                double_co_serial += 1

                double_group_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:DCO"
                    f"{double_co_serial:06d}"
                )

                double_co_boundary_map[
                    left_boundary
                ] = (
                    double_group_id,
                    "left_boundary",
                )

                double_co_boundary_map[
                    right_boundary
                ] = (
                    double_group_id,
                    "right_boundary",
                )

            # ------------------------------------------------
            # 9D. All unconsumed transitions are COs.
            # ------------------------------------------------

            for boundary_index, boundary in enumerate(
                boundaries
            ):
                if (
                    boundary_index
                    in consumed_boundaries
                ):
                    continue

                event_serial += 1
                co_serial += 1

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_serial:07d}"
                )

                co_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:CO"
                    f"{co_serial:07d}"
                )

                if (
                    boundary_index
                    in double_co_boundary_map
                ):
                    (
                        double_group_id,
                        co_context,
                    ) = (
                        double_co_boundary_map[
                            boundary_index
                        ]
                    )
                else:
                    double_group_id = ""
                    co_context = "SIMPLE_CO"

                qc = summarize_boundaries(
                    boundaries,
                    [boundary_index],
                )

                event_start = int(
                    boundary[
                        "left_marker_position"
                    ]
                )

                event_end = int(
                    boundary[
                        "right_marker_position"
                    ]
                )

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_class": "CO",
                        "n_CO_events": 1,
                        "n_NCO_events": 0,
                        "n_GC_events": 0,
                        "state_before": (
                            boundary[
                                "left_state"
                            ]
                        ),
                        "state_after": (
                            boundary[
                                "right_state"
                            ]
                        ),
                        "event_interval_start": (
                            event_start
                        ),
                        "event_interval_end": (
                            event_end
                        ),
                        "event_interval_bp": (
                            event_end
                            - event_start
                        ),
                        "core_start": np.nan,
                        "core_end": np.nan,
                        "core_span_bp": np.nan,
                        "component_block_spans_bp": "",
                        "n_core_markers": np.nan,
                        "single_marker_tract": False,
                        "double_co_group_id": (
                            double_group_id
                        ),
                        **qc,
                    }
                )

                co_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "co_id": co_id,
                        "event_id": event_id,
                        "co_context": co_context,
                        "left_marker_position": (
                            event_start
                        ),
                        "right_marker_position": (
                            event_end
                        ),
                        "breakpoint_interval_bp": (
                            event_end
                            - event_start
                        ),
                        "breakpoint_midpoint": (
                            (
                                event_start
                                + event_end
                            )
                            / 2
                        ),
                        "double_co_group_id": (
                            double_group_id
                        ),
                        **qc,
                    }
                )


# ============================================================
# 10. Convert to DataFrames
# ============================================================

blocks_df = pd.DataFrame(
    block_rows
)

boundaries_df = pd.DataFrame(
    boundary_rows
)

events_df = pd.DataFrame(
    event_rows
)

co_df = pd.DataFrame(
    co_rows
)


# ============================================================
# 11. Detect shared events across offspring
# ============================================================

def add_shared_clusters(
    data: pd.DataFrame,
    group_columns: list[str],
    start_column: str,
    end_column: str,
    prefix: str,
):
    result = data.copy()

    result[
        f"{prefix}_cluster_id"
    ] = ""

    result[
        f"{prefix}_n_events"
    ] = 1

    result[
        f"{prefix}_n_offspring"
    ] = 1

    result[
        f"{prefix}_shared_across_offspring"
    ] = False

    cluster_serial = 0

    for _, subset in result.groupby(
        group_columns,
        sort=False,
    ):
        subset = subset.sort_values(
            [
                start_column,
                end_column,
            ]
        )

        current_indices = []
        current_end = None

        def assign_cluster(indices):
            nonlocal cluster_serial

            if not indices:
                return

            cluster_serial += 1

            cluster_id = (
                f"{prefix.upper()}"
                f"{cluster_serial:06d}"
            )

            n_events = len(indices)

            n_offspring = (
                result.loc[
                    indices,
                    "sample_id"
                ]
                .nunique()
            )

            result.loc[
                indices,
                f"{prefix}_cluster_id"
            ] = cluster_id

            result.loc[
                indices,
                f"{prefix}_n_events"
            ] = n_events

            result.loc[
                indices,
                f"{prefix}_n_offspring"
            ] = n_offspring

            result.loc[
                indices,
                f"{prefix}_shared_across_offspring"
            ] = (
                n_offspring > 1
            )

        for index, row in subset.iterrows():
            start = int(
                row[start_column]
            )

            end = int(
                row[end_column]
            )

            if not current_indices:
                current_indices = [index]
                current_end = end
                continue

            if (
                start
                <= current_end
                + SHARED_TOLERANCE_BP
            ):
                current_indices.append(
                    index
                )

                current_end = max(
                    current_end,
                    end,
                )
            else:
                assign_cluster(
                    current_indices
                )

                current_indices = [index]
                current_end = end

        assign_cluster(
            current_indices
        )

    return result


if not events_df.empty:
    events_df["event_group"] = np.where(
        events_df["n_CO_events"] > 0,
        "CO_CONTAINING",
        "GC_ONLY",
    )

    events_df = add_shared_clusters(
        events_df,
        group_columns=[
            "family_id",
            "chromosome",
            "event_group",
        ],
        start_column=(
            "event_interval_start"
        ),
        end_column=(
            "event_interval_end"
        ),
        prefix="event",
    )

if not co_df.empty:
    co_df = add_shared_clusters(
        co_df,
        group_columns=[
            "family_id",
            "chromosome",
        ],
        start_column=(
            "left_marker_position"
        ),
        end_column=(
            "right_marker_position"
        ),
        prefix="co",
    )


# ============================================================
# 12. Assign review status
# ============================================================

def assign_review_status(row):
    if (
        bool(row["phase_issue"])
        or bool(
            row[
                "overlaps_reference_gap"
            ]
        )
    ):
        return (
            "EXCLUDE_OR_RESOLVE_PHASE_OR_GAP"
        )

    if bool(
        row.get(
            "event_shared_across_offspring",
            False,
        )
    ):
        return (
            "MANUAL_REVIEW_SHARED_EVENT"
        )

    if (
        row[
            "maximum_breakpoint_interval_bp"
        ]
        > LARGE_BREAKPOINT_BP
    ):
        return (
            "MANUAL_REVIEW_LOW_RESOLUTION"
        )

    if (
        row["event_class"]
        != "CO"
        and bool(
            row["single_marker_tract"]
        )
    ):
        return (
            "MANUAL_REVIEW_SINGLE_MARKER_GC"
        )

    return "CANDIDATE_HIGH_CONFIDENCE"


if not events_df.empty:
    events_df["review_status"] = (
        events_df.apply(
            assign_review_status,
            axis=1,
        )
    )

if not co_df.empty and not events_df.empty:
    event_annotations = (
        events_df[
            [
                "event_id",
                "event_shared_across_offspring",
                "review_status",
            ]
        ]
        .drop_duplicates("event_id")
    )

    co_df = co_df.merge(
        event_annotations,
        on="event_id",
        how="left",
    )


# ============================================================
# 13. Summaries
# ============================================================

sample_summary = (
    events_df.groupby(
        [
            "family_id",
            "sample_id",
        ],
        as_index=False,
    )
    .agg(
        n_CO_events=(
            "n_CO_events",
            "sum"
        ),
        n_NCO_candidates=(
            "n_NCO_events",
            "sum"
        ),
        n_GC_candidates=(
            "n_GC_events",
            "sum"
        ),
        n_event_records=(
            "event_id",
            "size"
        ),
        n_gap_or_phase_flagged=(
            "review_status",
            lambda values: (
                values
                == (
                    "EXCLUDE_OR_RESOLVE_"
                    "PHASE_OR_GAP"
                )
            ).sum(),
        ),
        n_shared_event_records=(
            "event_shared_across_offspring",
            "sum"
        ),
    )
)

high_confidence = events_df.loc[
    events_df["review_status"]
    == "CANDIDATE_HIGH_CONFIDENCE"
]

high_confidence_summary = (
    high_confidence.groupby(
        [
            "family_id",
            "sample_id",
        ],
        as_index=False,
    )
    .agg(
        n_high_confidence_CO=(
            "n_CO_events",
            "sum"
        ),
        n_high_confidence_NCO=(
            "n_NCO_events",
            "sum"
        ),
        n_high_confidence_GC=(
            "n_GC_events",
            "sum"
        ),
    )
)

sample_summary = sample_summary.merge(
    high_confidence_summary,
    on=[
        "family_id",
        "sample_id",
    ],
    how="left",
)

for column in [
    "n_high_confidence_CO",
    "n_high_confidence_NCO",
    "n_high_confidence_GC",
]:
    sample_summary[column] = (
        sample_summary[column]
        .fillna(0)
        .astype(int)
    )

chromosome_summary = (
    events_df.groupby(
        [
            "family_id",
            "chromosome",
        ],
        as_index=False,
    )
    .agg(
        n_CO_events=(
            "n_CO_events",
            "sum"
        ),
        n_NCO_candidates=(
            "n_NCO_events",
            "sum"
        ),
        n_GC_candidates=(
            "n_GC_events",
            "sum"
        ),
    )
)

chromosome_summary[
    "chromosome_length"
] = chromosome_summary[
    "chromosome"
].map(
    chromosome_lengths
)

chromosome_summary[
    "CO_events_per_Mb_per_offspring"
] = (
    chromosome_summary["n_CO_events"]
    / chromosome_summary[
        "chromosome_length"
    ]
    * 1_000_000
    / chromosome_summary[
        "family_id"
    ].map(
        {
            family: len(
                family_offspring[family]
            )
            for family in families
        }
    )
)


# ============================================================
# 14. Save outputs
# ============================================================

blocks_df.to_csv(
    OUTPUT_ROOT
    / "haplotype_blocks.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

boundaries_df.to_csv(
    OUTPUT_ROOT
    / "raw_haplotype_switch_boundaries.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

events_df.to_csv(
    OUTPUT_ROOT
    / "classified_recombination_events.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

co_df.to_csv(
    OUTPUT_ROOT
    / "CO_breakpoints.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

events_df.loc[
    events_df["n_GC_events"] > 0
].to_csv(
    OUTPUT_ROOT
    / "gene_conversion_candidates.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

events_df.loc[
    events_df["review_status"]
    == "CANDIDATE_HIGH_CONFIDENCE"
].to_csv(
    OUTPUT_ROOT
    / "candidate_high_confidence_events.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

sample_summary.to_csv(
    SUMMARY_DIR
    / "recombination_counts_by_offspring.tsv",
    sep="\t",
    index=False,
)

chromosome_summary.to_csv(
    SUMMARY_DIR
    / "recombination_counts_by_chromosome.tsv",
    sep="\t",
    index=False,
)

events_df.loc[
    events_df[
        "event_shared_across_offspring"
    ]
].to_csv(
    SUMMARY_DIR
    / "shared_recombination_event_clusters.tsv",
    sep="\t",
    index=False,
)

print(
    f"Tract threshold: "
    f"{TRACT_THRESHOLD_BP:,} bp"
)

print(
    f"Classified event records: "
    f"{len(events_df):,}"
)

print(
    f"CO breakpoints: "
    f"{len(co_df):,}"
)

print(
    f"NCO candidates: "
    f"{events_df['n_NCO_events'].sum():,}"
)

print(
    f"CO-associated GC candidates:"
    f"{(events_df['event_class']== 'CO_ASSOCIATED_GC_CANDIDATE').sum():,}"
)

print("\nCounts by offspring:")
print(
    sample_summary.to_string(
        index=False
    )
)

print(
    f"\nOutputs saved under:\n"
    f"{OUTPUT_ROOT}"
)
