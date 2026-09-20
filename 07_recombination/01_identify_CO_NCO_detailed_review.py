from pathlib import Path
from collections import defaultdict
from bisect import bisect_right
from itertools import combinations, count

import numpy as np
import pandas as pd


# ============================================================
# 1. Parameters
# ============================================================

# Short internal haplotype blocks are classified as
# NCO/gene-conversion candidates.
NCO_MAX_TRACT_BP = 10_000

# Breakpoint intervals wider than this are flagged as
# low-resolution. This does not automatically make a
# simple CO low confidence.
LOW_RESOLUTION_BP = 100_000

# Exact or near-shared breakpoints are distinguished.
NEAR_SHARED_DISTANCE_BP = 10_000

# Minimum number of markers required in each stable
# flanking haplotype block.
MIN_FLANK_MARKERS = 3

# NCO or short CO-associated tracts should preferably
# contain at least two markers.
MIN_CORE_MARKERS = 2

# Phase support:
# 5:0 = 1.0
# 4:1 = 0.8
# 3:2 = 0.6
MIN_STRONG_PHASE_SUPPORT = 0.80

# Shared paired COs within this distance are particularly
# suspicious when accompanied by gaps or weak phase.
SUSPICIOUS_PAIRED_CO_MAX_SPAN_BP = 700_000

# Number of markers shown on each side in the
# event-marker context table.
CONTEXT_FLANK_MARKERS = 5


# ============================================================
# 2. Paths
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

REFERENCE_GAP_BED = (
    PROJECT
    / "07_recombination/reference_masks/"
      "reference_N_gaps.bed"
)

OUTPUT_ROOT = (
    PROJECT
    / "07_recombination/classified/"
      "majority_phase_5of5/"
      "detailed_review_10000bp"
)

SUMMARY_DIR = (
    OUTPUT_ROOT / "summary"
)

REVIEW_DIR = (
    OUTPUT_ROOT / "manual_review"
)

OUTPUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True
)

REVIEW_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 3. Metadata
# ============================================================

metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str
)

families = sorted(
    metadata["family_id"].unique()
)

family_offspring = {}

for family in families:
    offspring = (
        metadata.loc[
            (
                metadata["family_id"] == family
            )
            &
            (
                metadata["role"] == "offspring"
            ),
            "sample_id"
        ]
        .sort_values()
        .tolist()
    )

    if len(offspring) != 5:
        raise ValueError(
            f"{family}: expected 5 offspring, "
            f"found {len(offspring)}"
        )

    family_offspring[family] = offspring


# ============================================================
# 4. Chromosome order
# ============================================================

primary = pd.read_csv(
    PRIMARY_BED,
    sep="\t",
    header=None,
    names=[
        "chromosome",
        "start",
        "end",
    ]
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
        - primary["start"]
    )
)


# ============================================================
# 5. Reference-gap lookup
# ============================================================

gap_intervals = defaultdict(list)
gap_starts = {}

if (
    REFERENCE_GAP_BED.exists()
    and REFERENCE_GAP_BED.stat().st_size > 0
):
    gaps = pd.read_csv(
        REFERENCE_GAP_BED,
        sep="\t",
        header=None,
        names=[
            "chromosome",
            "start",
            "end",
        ]
    )

    gaps["chromosome"] = (
        gaps["chromosome"].astype(str)
    )

    for chromosome, subset in gaps.groupby(
        "chromosome"
    ):
        intervals = sorted(
            zip(
                subset["start"].astype(int),
                subset["end"].astype(int),
            )
        )

        gap_intervals[chromosome] = intervals

        gap_starts[chromosome] = [
            start
            for start, end in intervals
        ]


def overlaps_reference_gap(
    chromosome: str,
    start_1based: int,
    end_1based: int,
) -> bool:
    """
    Test overlap between a one-based inclusive interval
    and zero-based half-open BED intervals.
    """

    intervals = gap_intervals.get(
        chromosome,
        []
    )

    if not intervals:
        return False

    start0 = max(
        int(start_1based) - 1,
        0
    )

    end0 = int(end_1based)

    starts = gap_starts[chromosome]

    index = (
        bisect_right(
            starts,
            end0 - 1
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
# 6. Haplotype block construction
# ============================================================

def make_block(
    chromosome_data: pd.DataFrame,
    states: np.ndarray,
    start_index: int,
    end_index: int,
    block_index: int,
) -> dict:
    start_position = int(
        chromosome_data.loc[
            start_index,
            "position"
        ]
    )

    end_position = int(
        chromosome_data.loc[
            end_index,
            "position"
        ]
    )

    return {
        "block_index": int(block_index),
        "state": int(states[start_index]),
        "first_marker_index": int(start_index),
        "last_marker_index": int(end_index),
        "start_position": start_position,
        "end_position": end_position,
        "span_bp": (
            end_position
            - start_position
        ),
        "n_markers": (
            end_index
            - start_index
            + 1
        ),
    }


def build_haplotype_blocks(
    chromosome_data: pd.DataFrame,
    haplotype_column: str,
) -> list[dict]:
    states = pd.to_numeric(
        chromosome_data[haplotype_column],
        errors="coerce"
    )

    valid = states.isin(
        [0, 1]
    )

    if not valid.all():
        invalid = chromosome_data.loc[
            ~valid,
            [
                "chromosome",
                "position",
                haplotype_column,
            ]
        ]

        raise ValueError(
            "Invalid haplotype values detected:\n"
            f"{invalid.head(20)}"
        )

    states = states.astype(
        int
    ).to_numpy()

    blocks = []
    start_index = 0

    for index in range(
        1,
        len(chromosome_data)
    ):
        if states[index] != states[index - 1]:
            blocks.append(
                make_block(
                    chromosome_data,
                    states,
                    start_index,
                    index - 1,
                    len(blocks),
                )
            )

            start_index = index

    blocks.append(
        make_block(
            chromosome_data,
            states,
            start_index,
            len(chromosome_data) - 1,
            len(blocks),
        )
    )

    return blocks


# ============================================================
# 7. Haplotype boundaries
# ============================================================

def build_boundaries(
    chromosome_data: pd.DataFrame,
    blocks: list[dict],
    chromosome: str,
) -> list[dict]:
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

        left_marker_index = (
            left_block["last_marker_index"]
        )

        right_marker_index = (
            right_block["first_marker_index"]
        )

        left_position = int(
            chromosome_data.loc[
                left_marker_index,
                "position"
            ]
        )

        right_position = int(
            chromosome_data.loc[
                right_marker_index,
                "position"
            ]
        )

        phase_support = pd.to_numeric(
            pd.Series(
                [
                    chromosome_data.loc[
                        right_marker_index,
                        "phase_support_from_previous"
                    ]
                ]
            ),
            errors="coerce"
        ).iloc[0]

        phase_category = str(
            chromosome_data.loc[
                right_marker_index,
                "phase_category_from_previous"
            ]
        )

        boundaries.append(
            {
                "boundary_index": boundary_index,
                "left_block_index": boundary_index,
                "right_block_index": boundary_index + 1,
                "left_state": left_block["state"],
                "right_state": right_block["state"],
                "left_marker_position": left_position,
                "right_marker_position": right_position,
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
                "phase_support": phase_support,
                "phase_category": phase_category,
                "weak_phase": (
                    pd.isna(phase_support)
                    or
                    phase_support
                    < MIN_STRONG_PHASE_SUPPORT
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
) -> dict:
    selected = [
        boundaries[index]
        for index in boundary_indices
    ]

    supports = [
        boundary["phase_support"]
        for boundary in selected
        if pd.notna(
            boundary["phase_support"]
        )
    ]

    signatures = [
        (
            f"{boundary['left_marker_position']}-"
            f"{boundary['right_marker_position']}"
        )
        for boundary in selected
    ]

    return {
        "boundary_indices": ",".join(
            str(index)
            for index in boundary_indices
        ),
        "boundary_signature": ";".join(
            signatures
        ),
        "boundary_phase_supports": ",".join(
            str(boundary["phase_support"])
            for boundary in selected
        ),
        "minimum_phase_support": (
            min(supports)
            if supports
            else np.nan
        ),
        "n_weak_phase_boundaries": sum(
            boundary["weak_phase"]
            for boundary in selected
        ),
        "maximum_breakpoint_interval_bp": max(
            boundary["breakpoint_interval_bp"]
            for boundary in selected
        ),
        "boundary_overlaps_reference_gap": any(
            boundary["overlaps_reference_gap"]
            for boundary in selected
        ),
    }


# ============================================================
# 8. Event identification
# ============================================================

block_rows = []
boundary_rows = []
event_rows = []

event_counter = count(1)
paired_co_counter = count(1)

# Keep phased data for marker-context output.
phased_lookup = {}


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

    phased["chromosome"] = (
        phased["chromosome"].astype(str)
    )

    phased["position"] = pd.to_numeric(
        phased["position"],
        errors="raise"
    ).astype(int)

    offspring = family_offspring[
        family
    ]

    # --------------------------------------------------------
    # Validate strict 5/5 phase-support values
    # --------------------------------------------------------

    for chromosome, subset in phased.groupby(
        "chromosome"
    ):
        subset = (
            subset
            .sort_values("position")
            .reset_index(drop=True)
        )

        support = pd.to_numeric(
            subset[
                "phase_support_from_previous"
            ],
            errors="coerce"
        )

        values = (
            support.iloc[1:]
            .dropna()
            .to_numpy()
        )

        if len(values) > 0:
            expected = np.isclose(
                values[:, None],
                np.array(
                    [0.6, 0.8, 1.0]
                )[None, :]
            ).any(axis=1)

            if not expected.all():
                raise ValueError(
                    f"{family} {chromosome}: "
                    "unexpected phase support values: "
                    f"{np.unique(values[~expected])}"
                )

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

        phased_lookup[
            (
                family,
                chromosome,
            )
        ] = chromosome_data

        for child in offspring:
            haplotype_column = (
                f"{child}_haplotype"
            )

            blocks = build_haplotype_blocks(
                chromosome_data,
                haplotype_column,
            )

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

            consumed_boundaries = set()

            # =================================================
            # 8A. CO-associated complex short-tract pattern
            #
            # A | B-short | A-short | B
            #
            # This is counted as one CO, not one NCO plus CO.
            # =================================================

            for block_index in range(
                1,
                len(blocks) - 2
            ):
                left_flank = blocks[
                    block_index - 1
                ]

                short_1 = blocks[
                    block_index
                ]

                short_2 = blocks[
                    block_index + 1
                ]

                right_flank = blocks[
                    block_index + 2
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
                    short_1["span_bp"]
                    <= NCO_MAX_TRACT_BP
                    and
                    short_2["span_bp"]
                    <= NCO_MAX_TRACT_BP
                ):
                    continue

                if (
                    left_flank["state"]
                    == right_flank["state"]
                ):
                    continue

                event_number = next(
                    event_counter
                )

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_number:07d}"
                )

                event_start = int(
                    left_flank["end_position"]
                )

                event_end = int(
                    right_flank["start_position"]
                )

                boundary_summary = (
                    summarize_boundaries(
                        boundaries,
                        boundary_indices,
                    )
                )

                component_markers = [
                    short_1["n_markers"],
                    short_2["n_markers"],
                ]

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_type": "CO",
                        "event_context": (
                            "CO_ASSOCIATED_COMPLEX_TRACT"
                        ),
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
                        "event_midpoint": (
                            (
                                event_start
                                + event_end
                            )
                            / 2
                        ),
                        "core_start": (
                            short_1[
                                "start_position"
                            ]
                        ),
                        "core_end": (
                            short_2[
                                "end_position"
                            ]
                        ),
                        "core_span_bp": (
                            short_2[
                                "end_position"
                            ]
                            - short_1[
                                "start_position"
                            ]
                        ),
                        "component_core_spans_bp": (
                            f"{short_1['span_bp']},"
                            f"{short_2['span_bp']}"
                        ),
                        "component_core_markers": (
                            f"{short_1['n_markers']},"
                            f"{short_2['n_markers']}"
                        ),
                        "minimum_component_core_markers": (
                            min(component_markers)
                        ),
                        "n_core_markers": (
                            short_1["n_markers"]
                            + short_2["n_markers"]
                        ),
                        "left_flank_n_markers": (
                            left_flank[
                                "n_markers"
                            ]
                        ),
                        "right_flank_n_markers": (
                            right_flank[
                                "n_markers"
                            ]
                        ),
                        "left_flank_span_bp": (
                            left_flank[
                                "span_bp"
                            ]
                        ),
                        "right_flank_span_bp": (
                            right_flank[
                                "span_bp"
                            ]
                        ),
                        "minimum_flank_markers": min(
                            left_flank["n_markers"],
                            right_flank["n_markers"],
                        ),
                        "minimum_flank_span_bp": min(
                            left_flank["span_bp"],
                            right_flank["span_bp"],
                        ),
                        "paired_co_group_id": "",
                        "paired_co_tract_span_bp": (
                            np.nan
                        ),
                        "paired_co_boundary_role": "",
                        "event_overlaps_reference_gap": (
                            boundary_summary[
                                "boundary_overlaps_reference_gap"
                            ]
                            or
                            overlaps_reference_gap(
                                chromosome,
                                event_start,
                                event_end,
                            )
                        ),
                        **boundary_summary,
                    }
                )

                consumed_boundaries.update(
                    boundary_indices
                )

            # =================================================
            # 8B. NCO candidate
            #
            # A | B-short | A
            # =================================================

            for block_index in range(
                1,
                len(blocks) - 1
            ):
                left_flank = blocks[
                    block_index - 1
                ]

                core_block = blocks[
                    block_index
                ]

                right_flank = blocks[
                    block_index + 1
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
                    left_flank["state"]
                    != right_flank["state"]
                ):
                    continue

                if (
                    core_block["span_bp"]
                    > NCO_MAX_TRACT_BP
                ):
                    continue

                event_number = next(
                    event_counter
                )

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_number:07d}"
                )

                event_start = int(
                    left_flank["end_position"]
                )

                event_end = int(
                    right_flank["start_position"]
                )

                boundary_summary = (
                    summarize_boundaries(
                        boundaries,
                        boundary_indices,
                    )
                )

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_type": "NCO",
                        "event_context": (
                            "INTERNAL_SHORT_TRACT"
                        ),
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
                        "event_midpoint": (
                            (
                                event_start
                                + event_end
                            )
                            / 2
                        ),
                        "core_start": (
                            core_block[
                                "start_position"
                            ]
                        ),
                        "core_end": (
                            core_block[
                                "end_position"
                            ]
                        ),
                        "core_span_bp": (
                            core_block[
                                "span_bp"
                            ]
                        ),
                        "component_core_spans_bp": (
                            str(
                                core_block[
                                    "span_bp"
                                ]
                            )
                        ),
                        "component_core_markers": (
                            str(
                                core_block[
                                    "n_markers"
                                ]
                            )
                        ),
                        "minimum_component_core_markers": (
                            core_block[
                                "n_markers"
                            ]
                        ),
                        "n_core_markers": (
                            core_block[
                                "n_markers"
                            ]
                        ),
                        "left_flank_n_markers": (
                            left_flank[
                                "n_markers"
                            ]
                        ),
                        "right_flank_n_markers": (
                            right_flank[
                                "n_markers"
                            ]
                        ),
                        "left_flank_span_bp": (
                            left_flank[
                                "span_bp"
                            ]
                        ),
                        "right_flank_span_bp": (
                            right_flank[
                                "span_bp"
                            ]
                        ),
                        "minimum_flank_markers": min(
                            left_flank["n_markers"],
                            right_flank["n_markers"],
                        ),
                        "minimum_flank_span_bp": min(
                            left_flank["span_bp"],
                            right_flank["span_bp"],
                        ),
                        "paired_co_group_id": "",
                        "paired_co_tract_span_bp": (
                            np.nan
                        ),
                        "paired_co_boundary_role": "",
                        "event_overlaps_reference_gap": (
                            boundary_summary[
                                "boundary_overlaps_reference_gap"
                            ]
                            or
                            overlaps_reference_gap(
                                chromosome,
                                event_start,
                                event_end,
                            )
                        ),
                        **boundary_summary,
                    }
                )

                consumed_boundaries.update(
                    boundary_indices
                )

            # =================================================
            # 8C. Identify paired CO boundaries surrounding
            #     a long internal haplotype block.
            #
            # A | B-long | A
            # =================================================

            paired_boundary_map = {}

            for block_index in range(
                1,
                len(blocks) - 1
            ):
                central_block = blocks[
                    block_index
                ]

                left_flank = blocks[
                    block_index - 1
                ]

                right_flank = blocks[
                    block_index + 1
                ]

                if (
                    central_block["span_bp"]
                    <= NCO_MAX_TRACT_BP
                ):
                    continue

                if (
                    left_flank["state"]
                    != right_flank["state"]
                ):
                    continue

                left_boundary = (
                    block_index - 1
                )

                right_boundary = (
                    block_index
                )

                if (
                    left_boundary
                    in consumed_boundaries
                    or
                    right_boundary
                    in consumed_boundaries
                ):
                    continue

                if (
                    left_boundary
                    in paired_boundary_map
                    or
                    right_boundary
                    in paired_boundary_map
                ):
                    continue

                paired_number = next(
                    paired_co_counter
                )

                group_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:DCO"
                    f"{paired_number:06d}"
                )

                paired_boundary_map[
                    left_boundary
                ] = {
                    "group_id": group_id,
                    "role": "left_boundary",
                    "tract_span_bp": (
                        central_block[
                            "span_bp"
                        ]
                    ),
                }

                paired_boundary_map[
                    right_boundary
                ] = {
                    "group_id": group_id,
                    "role": "right_boundary",
                    "tract_span_bp": (
                        central_block[
                            "span_bp"
                        ]
                    ),
                }

            # =================================================
            # 8D. Remaining boundaries are CO candidates
            # =================================================

            for boundary_index, boundary in enumerate(
                boundaries
            ):
                if (
                    boundary_index
                    in consumed_boundaries
                ):
                    continue

                event_number = next(
                    event_counter
                )

                event_id = (
                    f"{family}:{child}:"
                    f"{chromosome}:E"
                    f"{event_number:07d}"
                )

                left_block = blocks[
                    boundary[
                        "left_block_index"
                    ]
                ]

                right_block = blocks[
                    boundary[
                        "right_block_index"
                    ]
                ]

                paired_info = (
                    paired_boundary_map.get(
                        boundary_index
                    )
                )

                if paired_info is None:
                    event_context = "SIMPLE_CO"
                    paired_group_id = ""
                    paired_role = ""
                    paired_span = np.nan
                else:
                    event_context = "PAIRED_CO"
                    paired_group_id = (
                        paired_info["group_id"]
                    )
                    paired_role = (
                        paired_info["role"]
                    )
                    paired_span = (
                        paired_info[
                            "tract_span_bp"
                        ]
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

                boundary_summary = (
                    summarize_boundaries(
                        boundaries,
                        [boundary_index],
                    )
                )

                event_rows.append(
                    {
                        "family_id": family,
                        "sample_id": child,
                        "chromosome": chromosome,
                        "event_id": event_id,
                        "event_type": "CO",
                        "event_context": (
                            event_context
                        ),
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
                        "event_midpoint": (
                            (
                                event_start
                                + event_end
                            )
                            / 2
                        ),
                        "core_start": np.nan,
                        "core_end": np.nan,
                        "core_span_bp": np.nan,
                        "component_core_spans_bp": "",
                        "component_core_markers": "",
                        "minimum_component_core_markers": (
                            np.nan
                        ),
                        "n_core_markers": np.nan,
                        "left_flank_n_markers": (
                            left_block[
                                "n_markers"
                            ]
                        ),
                        "right_flank_n_markers": (
                            right_block[
                                "n_markers"
                            ]
                        ),
                        "left_flank_span_bp": (
                            left_block[
                                "span_bp"
                            ]
                        ),
                        "right_flank_span_bp": (
                            right_block[
                                "span_bp"
                            ]
                        ),
                        "minimum_flank_markers": min(
                            left_block["n_markers"],
                            right_block["n_markers"],
                        ),
                        "minimum_flank_span_bp": min(
                            left_block["span_bp"],
                            right_block["span_bp"],
                        ),
                        "paired_co_group_id": (
                            paired_group_id
                        ),
                        "paired_co_tract_span_bp": (
                            paired_span
                        ),
                        "paired_co_boundary_role": (
                            paired_role
                        ),
                        "event_overlaps_reference_gap": (
                            boundary[
                                "overlaps_reference_gap"
                            ]
                        ),
                        **boundary_summary,
                    }
                )


# ============================================================
# 9. Convert to DataFrames
# ============================================================

blocks_df = pd.DataFrame(
    block_rows
)

boundaries_df = pd.DataFrame(
    boundary_rows
)

events = pd.DataFrame(
    event_rows
)

if events.empty:
    raise RuntimeError(
        "No CO or NCO events were identified."
    )


# ============================================================
# 10. Detect exact, near and broad shared events
# ============================================================

shared_information = {
    index: {
        "exact": set(),
        "near": set(),
        "broad": set(),
    }
    for index in events.index
}

shared_pair_rows = []


for _, group in events.groupby(
    [
        "family_id",
        "chromosome",
        "event_type",
    ],
    sort=False
):
    indices = group.index.tolist()

    for index_a, index_b in combinations(
        indices,
        2
    ):
        row_a = events.loc[index_a]
        row_b = events.loc[index_b]

        if (
            row_a["sample_id"]
            == row_b["sample_id"]
        ):
            continue

        exact = (
            row_a["boundary_signature"]
            == row_b["boundary_signature"]
        )

        near = (
            abs(
                row_a["event_midpoint"]
                - row_b["event_midpoint"]
            )
            <= NEAR_SHARED_DISTANCE_BP
        )

        overlap = (
            row_a["event_interval_start"]
            <= row_b["event_interval_end"]
            and
            row_b["event_interval_start"]
            <= row_a["event_interval_end"]
        )

        if exact:
            relation = "exact"

        elif near:
            relation = "near"

        elif overlap:
            relation = "broad"

        else:
            continue

        shared_information[
            index_a
        ][relation].add(
            row_b["sample_id"]
        )

        shared_information[
            index_b
        ][relation].add(
            row_a["sample_id"]
        )

        shared_pair_rows.append(
            {
                "family_id": (
                    row_a["family_id"]
                ),
                "chromosome": (
                    row_a["chromosome"]
                ),
                "event_type": (
                    row_a["event_type"]
                ),
                "event_id_1": (
                    row_a["event_id"]
                ),
                "sample_id_1": (
                    row_a["sample_id"]
                ),
                "event_id_2": (
                    row_b["event_id"]
                ),
                "sample_id_2": (
                    row_b["sample_id"]
                ),
                "shared_relation": relation,
                "midpoint_distance_bp": abs(
                    row_a["event_midpoint"]
                    - row_b["event_midpoint"]
                ),
                "event_1_start": (
                    row_a[
                        "event_interval_start"
                    ]
                ),
                "event_1_end": (
                    row_a[
                        "event_interval_end"
                    ]
                ),
                "event_2_start": (
                    row_b[
                        "event_interval_start"
                    ]
                ),
                "event_2_end": (
                    row_b[
                        "event_interval_end"
                    ]
                ),
            }
        )


def join_samples(values):
    return ",".join(
        sorted(values)
    )


events[
    "exact_shared_partner_samples"
] = [
    join_samples(
        shared_information[index]["exact"]
    )
    for index in events.index
]

events[
    "near_shared_partner_samples"
] = [
    join_samples(
        shared_information[index]["near"]
    )
    for index in events.index
]

events[
    "broad_overlap_partner_samples"
] = [
    join_samples(
        shared_information[index]["broad"]
    )
    for index in events.index
]

events["flag_exact_shared"] = [
    len(
        shared_information[index]["exact"]
    ) > 0
    for index in events.index
]

events["flag_near_shared"] = [
    len(
        shared_information[index]["near"]
    ) > 0
    for index in events.index
]

events["flag_broad_overlap"] = [
    len(
        shared_information[index]["broad"]
    ) > 0
    for index in events.index
]

events["shared_partner_samples"] = [
    join_samples(
        shared_information[index]["exact"]
        |
        shared_information[index]["near"]
        |
        shared_information[index]["broad"]
    )
    for index in events.index
]

events["shared_n_other_offspring"] = [
    len(
        shared_information[index]["exact"]
        |
        shared_information[index]["near"]
        |
        shared_information[index]["broad"]
    )
    for index in events.index
]


def classify_shared_type(row):
    if row["flag_exact_shared"]:
        return "exact"

    if row["flag_near_shared"]:
        return "near"

    if row["flag_broad_overlap"]:
        return "broad"

    return "none"


events["shared_event_type"] = (
    events.apply(
        classify_shared_type,
        axis=1
    )
)


# ============================================================
# 11. Raw event-quality flags
# ============================================================

events["flag_weak_phase"] = (
    pd.to_numeric(
        events["minimum_phase_support"],
        errors="coerce"
    ).fillna(0)
    < MIN_STRONG_PHASE_SUPPORT
)

events["flag_reference_gap"] = (
    events[
        "event_overlaps_reference_gap"
    ].astype(bool)
)

events["flag_low_resolution"] = (
    pd.to_numeric(
        events[
            "maximum_breakpoint_interval_bp"
        ],
        errors="coerce"
    ).fillna(np.inf)
    > LOW_RESOLUTION_BP
)

events[
    "flag_wide_NCO_maximum_tract"
] = (
    (
        events["event_type"] == "NCO"
    )
    &
    (
        events["event_interval_bp"]
        > LOW_RESOLUTION_BP
    )
)

events[
    "flag_single_marker_core"
] = (
    (
        events["event_type"] == "NCO"
    )
    &
    (
        pd.to_numeric(
            events["n_core_markers"],
            errors="coerce"
        ).fillna(0)
        < MIN_CORE_MARKERS
    )
) | (
    (
        events["event_context"]
        == "CO_ASSOCIATED_COMPLEX_TRACT"
    )
    &
    (
        pd.to_numeric(
            events[
                "minimum_component_core_markers"
            ],
            errors="coerce"
        ).fillna(0)
        < MIN_CORE_MARKERS
    )
)

events[
    "flag_insufficient_flank_support"
] = (
    pd.to_numeric(
        events["minimum_flank_markers"],
        errors="coerce"
    ).fillna(0)
    < MIN_FLANK_MARKERS
)

events["flag_paired_CO"] = (
    events["event_context"]
    == "PAIRED_CO"
)

events["flag_complex_CO_pattern"] = (
    events["event_context"]
    == "CO_ASSOCIATED_COMPLEX_TRACT"
)

events[
    "flag_suspicious_paired_CO"
] = (
    events["flag_paired_CO"]
    &
    (
        pd.to_numeric(
            events[
                "paired_co_tract_span_bp"
            ],
            errors="coerce"
        ).fillna(np.inf)
        <= SUSPICIOUS_PAIRED_CO_MAX_SPAN_BP
    )
    &
    (
        events["flag_exact_shared"]
        |
        events["flag_near_shared"]
    )
    &
    (
        events["flag_weak_phase"]
        |
        events["flag_reference_gap"]
    )
)


# ============================================================
# 12. Automatically assess whether each flag truly
#     lowers event confidence
# ============================================================

def assess_event(row):
    assessments = {}

    flank_adequate = not bool(
        row[
            "flag_insufficient_flank_support"
        ]
    )

    strong_phase = not bool(
        row["flag_weak_phase"]
    )

    simple_co = (
        row["event_type"] == "CO"
        and
        row["event_context"] == "SIMPLE_CO"
    )

    # --------------------------------------------------------
    # Weak phase
    # --------------------------------------------------------

    if not row["flag_weak_phase"]:
        assessments["weak_phase"] = "none"

    elif (
        simple_co
        and
        flank_adequate
        and
        not row["flag_reference_gap"]
        and
        not row["flag_exact_shared"]
        and
        not row[
            "flag_suspicious_paired_CO"
        ]
    ):
        # A persistent simple CO supported by long stable
        # blocks may still be real under a 3:2 phase vote.
        assessments["weak_phase"] = "note"

    else:
        assessments["weak_phase"] = "low"

    # --------------------------------------------------------
    # Reference gap
    # --------------------------------------------------------

    if not row["flag_reference_gap"]:
        assessments["reference_gap"] = "none"

    elif (
        simple_co
        and
        strong_phase
        and
        flank_adequate
        and
        not row["flag_exact_shared"]
    ):
        # The CO count can be reliable even though its
        # exact breakpoint inside the gap is unresolved.
        assessments["reference_gap"] = "note"

    else:
        assessments["reference_gap"] = "low"

    # --------------------------------------------------------
    # Low breakpoint resolution
    # --------------------------------------------------------

    if not row["flag_low_resolution"]:
        assessments["low_resolution"] = "none"

    elif (
        simple_co
        and
        flank_adequate
        and
        not row[
            "flag_suspicious_paired_CO"
        ]
    ):
        assessments["low_resolution"] = "note"

    else:
        assessments["low_resolution"] = "low"

    # --------------------------------------------------------
    # NCO maximum possible tract is too wide
    # --------------------------------------------------------

    if not row[
        "flag_wide_NCO_maximum_tract"
    ]:
        assessments[
            "wide_NCO_maximum_tract"
        ] = "none"

    else:
        assessments[
            "wide_NCO_maximum_tract"
        ] = "low"

    # --------------------------------------------------------
    # Core marker support
    # --------------------------------------------------------

    if not row[
        "flag_single_marker_core"
    ]:
        assessments[
            "single_marker_core"
        ] = "none"

    else:
        assessments[
            "single_marker_core"
        ] = "low"

    # --------------------------------------------------------
    # Flanking block support
    # --------------------------------------------------------

    if not row[
        "flag_insufficient_flank_support"
    ]:
        assessments[
            "insufficient_flank_support"
        ] = "none"

    else:
        assessments[
            "insufficient_flank_support"
        ] = "low"

    # --------------------------------------------------------
    # Exact shared event
    # --------------------------------------------------------

    if not row["flag_exact_shared"]:
        assessments["exact_shared"] = "none"

    elif (
        simple_co
        and
        strong_phase
        and
        flank_adequate
        and
        not row["flag_reference_gap"]
    ):
        # Could be a real hotspot or limited marker
        # resolution rather than an artifact.
        assessments["exact_shared"] = "note"

    else:
        assessments["exact_shared"] = "low"

    # --------------------------------------------------------
    # Near-shared event
    # --------------------------------------------------------

    if not row["flag_near_shared"]:
        assessments["near_shared"] = "none"

    elif (
        simple_co
        and
        flank_adequate
        and
        not row[
            "flag_suspicious_paired_CO"
        ]
    ):
        assessments["near_shared"] = "note"

    else:
        assessments["near_shared"] = "low"

    # --------------------------------------------------------
    # Broad overlap
    # --------------------------------------------------------

    if not row["flag_broad_overlap"]:
        assessments["broad_overlap"] = "none"

    else:
        # Broad interval overlap alone does not make an
        # event unreliable.
        assessments["broad_overlap"] = "note"

    # --------------------------------------------------------
    # Paired CO
    # --------------------------------------------------------

    if not row["flag_paired_CO"]:
        assessments["paired_CO"] = "none"

    else:
        # A paired CO is biologically possible and is not
        # automatically low confidence.
        assessments["paired_CO"] = "note"

    # --------------------------------------------------------
    # Suspicious paired CO
    # --------------------------------------------------------

    if not row[
        "flag_suspicious_paired_CO"
    ]:
        assessments[
            "suspicious_paired_CO"
        ] = "none"

    else:
        assessments[
            "suspicious_paired_CO"
        ] = "low"

    # --------------------------------------------------------
    # Complex CO-associated short-tract pattern
    # --------------------------------------------------------

    if not row[
        "flag_complex_CO_pattern"
    ]:
        assessments[
            "complex_CO_pattern"
        ] = "none"

    elif (
        strong_phase
        and
        flank_adequate
        and
        not row["flag_reference_gap"]
        and
        not row[
            "flag_single_marker_core"
        ]
        and
        not row["flag_exact_shared"]
    ):
        assessments[
            "complex_CO_pattern"
        ] = "note"

    else:
        assessments[
            "complex_CO_pattern"
        ] = "low"

    return assessments


assessment_rows = []

reason_codes = [
    "weak_phase",
    "reference_gap",
    "low_resolution",
    "wide_NCO_maximum_tract",
    "single_marker_core",
    "insufficient_flank_support",
    "exact_shared",
    "near_shared",
    "broad_overlap",
    "paired_CO",
    "suspicious_paired_CO",
    "complex_CO_pattern",
]


for _, row in events.iterrows():
    assessments = assess_event(
        row
    )

    triggered = [
        reason
        for reason in reason_codes
        if assessments[reason] != "none"
    ]

    effective_low = [
        reason
        for reason in reason_codes
        if assessments[reason] == "low"
    ]

    notes = [
        reason
        for reason in reason_codes
        if assessments[reason] == "note"
    ]

    initial_confidence = (
        "low"
        if triggered
        else "high"
    )

    final_confidence = (
        "low"
        if effective_low
        else "high"
    )

    if effective_low:
        review_status = (
            "LOW_MANUAL_REVIEW"
        )

    elif triggered:
        review_status = (
            "HIGH_AFTER_AUTO_REVIEW"
        )

    else:
        review_status = "HIGH"

    assessment_row = {
        "initial_confidence": (
            initial_confidence
        ),
        "raw_flag_reason_codes": (
            ";".join(triggered)
            if triggered
            else "."
        ),
        "effective_low_reason_codes": (
            ";".join(effective_low)
            if effective_low
            else "."
        ),
        "non_disqualifying_note_codes": (
            ";".join(notes)
            if notes
            else "."
        ),
        "final_confidence": (
            final_confidence
        ),
        "review_status": review_status,
        "confidence_changed_by_review": (
            initial_confidence
            != final_confidence
        ),
    }

    for reason in reason_codes:
        assessment_row[
            f"assessment_{reason}"
        ] = assessments[reason]

    assessment_rows.append(
        assessment_row
    )


assessment_df = pd.DataFrame(
    assessment_rows,
    index=events.index
)

events = pd.concat(
    [
        events,
        assessment_df,
    ],
    axis=1
)


# ============================================================
# 13. Manual-review priority
# ============================================================

def assign_review_priority(row):
    if (
        row["review_status"]
        == "LOW_MANUAL_REVIEW"
    ):
        if (
            row["event_type"] == "NCO"
            or
            row[
                "flag_suspicious_paired_CO"
            ]
            or
            row["flag_exact_shared"]
            or
            row[
                "flag_complex_CO_pattern"
            ]
        ):
            return 1

        return 2

    if (
        row["review_status"]
        == "HIGH_AFTER_AUTO_REVIEW"
    ):
        return 3

    return 4


events["manual_review_priority"] = (
    events.apply(
        assign_review_priority,
        axis=1
    )
)

events[
    "manual_BAM_review_recommended"
] = (
    events["final_confidence"]
    == "low"
) | (
    events["event_type"] == "NCO"
) | (
    events[
        "flag_complex_CO_pattern"
    ]
)

# Empty columns for later manual review.
events["manual_review_decision"] = ""
events["manual_review_reason"] = ""
events["manual_review_notes"] = ""


# ============================================================
# 14. Sort the main table
# ============================================================

events["_chromosome_rank"] = (
    events["chromosome"]
    .map(chromosome_rank)
)

events = (
    events.sort_values(
        [
            "manual_review_priority",
            "family_id",
            "sample_id",
            "_chromosome_rank",
            "event_interval_start",
        ]
    )
    .drop(
        columns="_chromosome_rank"
    )
    .reset_index(drop=True)
)


# ============================================================
# 15. Save one main event table
# ============================================================

MAIN_OUTPUT = (
    OUTPUT_ROOT
    / "CO_NCO_events.detailed_review.tsv.gz"
)

events.to_csv(
    MAIN_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip"
)


# ============================================================
# 16. Manual-review queue
# ============================================================

review_queue = events.loc[
    (
        events["review_status"]
        != "HIGH"
    )
    |
    (
        events["event_type"] == "NCO"
    )
].copy()

review_queue.to_csv(
    REVIEW_DIR
    / "CO_NCO_manual_review_queue.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)


# ============================================================
# 17. BED file for IGV/manual inspection
# ============================================================

review_bed = review_queue.copy()

review_bed["bed_start"] = (
    review_bed[
        "event_interval_start"
    ].astype(int)
    - 1
).clip(lower=0)

review_bed["bed_end"] = (
    review_bed[
        "event_interval_end"
    ].astype(int)
)

review_bed["bed_name"] = (
    review_bed["event_id"]
    + "|"
    + review_bed["event_type"]
    + "|"
    + review_bed["review_status"]
    + "|"
    + review_bed[
        "effective_low_reason_codes"
    ]
)

review_bed[
    [
        "chromosome",
        "bed_start",
        "bed_end",
        "bed_name",
    ]
].to_csv(
    REVIEW_DIR
    / "CO_NCO_manual_review_regions.bed",
    sep="\t",
    header=False,
    index=False
)


# ============================================================
# 18. Marker context for manual inspection
# ============================================================

marker_context_rows = []

for _, event in review_queue.iterrows():
    key = (
        event["family_id"],
        event["chromosome"],
    )

    chromosome_data = (
        phased_lookup[key]
    )

    positions = chromosome_data[
        "position"
    ].to_numpy()

    left_index = int(
        np.searchsorted(
            positions,
            event["event_interval_start"],
            side="left"
        )
    )

    right_index = int(
        np.searchsorted(
            positions,
            event["event_interval_end"],
            side="right"
        )
    )

    context_start = max(
        0,
        left_index
        - CONTEXT_FLANK_MARKERS
    )

    context_end = min(
        len(chromosome_data),
        right_index
        + CONTEXT_FLANK_MARKERS
    )

    context = chromosome_data.iloc[
        context_start:context_end
    ].copy()

    child = event["sample_id"]

    haplotype_column = (
        f"{child}_haplotype"
    )

    for _, marker in context.iterrows():
        position = int(
            marker["position"]
        )

        core_start = event["core_start"]
        core_end = event["core_end"]

        if (
            pd.notna(core_start)
            and
            pd.notna(core_end)
            and
            core_start
            <= position
            <= core_end
        ):
            context_role = "core"

        elif (
            position
            < event["event_interval_start"]
        ):
            context_role = "left_flank"

        elif (
            position
            > event["event_interval_end"]
        ):
            context_role = "right_flank"

        else:
            context_role = (
                "event_interval"
            )

        marker_context_rows.append(
            {
                "event_id": (
                    event["event_id"]
                ),
                "family_id": (
                    event["family_id"]
                ),
                "sample_id": child,
                "chromosome": (
                    event["chromosome"]
                ),
                "event_type": (
                    event["event_type"]
                ),
                "event_context": (
                    event["event_context"]
                ),
                "review_status": (
                    event["review_status"]
                ),
                "effective_low_reason_codes": (
                    event[
                        "effective_low_reason_codes"
                    ]
                ),
                "position": position,
                "context_role": context_role,
                "reference_allele": (
                    marker.get(
                        "reference_allele",
                        np.nan
                    )
                ),
                "alternative_allele": (
                    marker.get(
                        "alternative_allele",
                        np.nan
                    )
                ),
                "offspring_marker_state": (
                    marker.get(
                        child,
                        np.nan
                    )
                ),
                "offspring_haplotype": (
                    marker.get(
                        haplotype_column,
                        np.nan
                    )
                ),
                "phase_bit": (
                    marker.get(
                        "phase_bit",
                        np.nan
                    )
                ),
                "phase_support_from_previous": (
                    marker.get(
                        "phase_support_from_previous",
                        np.nan
                    )
                ),
                "phase_category_from_previous": (
                    marker.get(
                        "phase_category_from_previous",
                        ""
                    )
                ),
            }
        )


marker_context = pd.DataFrame(
    marker_context_rows
)

marker_context.to_csv(
    REVIEW_DIR
    / "CO_NCO_event_marker_context.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)


# ============================================================
# 19. Save intermediate/supporting files
# ============================================================

blocks_df.to_csv(
    OUTPUT_ROOT
    / "haplotype_blocks.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

boundaries_df.to_csv(
    OUTPUT_ROOT
    / "haplotype_switch_boundaries.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip"
)

shared_pairs = pd.DataFrame(
    shared_pair_rows
)

shared_pairs.to_csv(
    REVIEW_DIR
    / "shared_event_pairs.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 20. Count summary by offspring
# ============================================================

count_summary = (
    events.groupby(
        [
            "family_id",
            "sample_id",
            "event_type",
            "final_confidence",
            "review_status",
        ],
        as_index=False
    )
    .size()
    .rename(
        columns={
            "size": "n_events"
        }
    )
)

count_summary.to_csv(
    SUMMARY_DIR
    / "event_counts_by_offspring_and_status.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 21. Overall review-status summary
# ============================================================

review_summary = (
    events.groupby(
        [
            "event_type",
            "review_status",
            "final_confidence",
        ],
        as_index=False
    )
    .size()
    .rename(
        columns={
            "size": "n_events"
        }
    )
)

review_summary.to_csv(
    SUMMARY_DIR
    / "review_status_summary.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 22. Reason trigger/effect summary
# ============================================================

reason_summary_rows = []

for reason in reason_codes:
    column = (
        f"assessment_{reason}"
    )

    for assessment, n_events in (
        events[column]
        .value_counts()
        .items()
    ):
        reason_summary_rows.append(
            {
                "reason_code": reason,
                "assessment": assessment,
                "n_events": n_events,
            }
        )

reason_summary = pd.DataFrame(
    reason_summary_rows
)

reason_summary.to_csv(
    SUMMARY_DIR
    / "reason_trigger_and_effect_summary.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 23. Save parameter table
# ============================================================

parameters = pd.DataFrame(
    [
        {
            "parameter": (
                "NCO_MAX_TRACT_BP"
            ),
            "value": NCO_MAX_TRACT_BP,
        },
        {
            "parameter": (
                "LOW_RESOLUTION_BP"
            ),
            "value": LOW_RESOLUTION_BP,
        },
        {
            "parameter": (
                "NEAR_SHARED_DISTANCE_BP"
            ),
            "value": (
                NEAR_SHARED_DISTANCE_BP
            ),
        },
        {
            "parameter": (
                "MIN_FLANK_MARKERS"
            ),
            "value": MIN_FLANK_MARKERS,
        },
        {
            "parameter": (
                "MIN_CORE_MARKERS"
            ),
            "value": MIN_CORE_MARKERS,
        },
        {
            "parameter": (
                "MIN_STRONG_PHASE_SUPPORT"
            ),
            "value": (
                MIN_STRONG_PHASE_SUPPORT
            ),
        },
        {
            "parameter": (
                "SUSPICIOUS_PAIRED_CO_"
                "MAX_SPAN_BP"
            ),
            "value": (
                SUSPICIOUS_PAIRED_CO_MAX_SPAN_BP
            ),
        },
    ]
)

parameters.to_csv(
    SUMMARY_DIR
    / "analysis_parameters.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 24. Print
# ============================================================

print(
    "Detailed CO/NCO identification completed."
)

print(
    f"\nMain table:\n{MAIN_OUTPUT}"
)

print("\nReview summary:")

print(
    review_summary.to_string(
        index=False
    )
)

print(
    "\nConfidence changes:"
)

print(
    events[
        [
            "initial_confidence",
            "final_confidence",
            "review_status",
        ]
    ]
    .value_counts()
    .rename("n_events")
    .reset_index()
    .to_string(index=False)
)

print(
    f"\nManual-review files:\n"
    f"{REVIEW_DIR}"
)
