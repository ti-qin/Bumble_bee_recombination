from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path("/localdata/qinti/Project/Bee")

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

MARKER_ROOT = (
    PROJECT
    / "05_markers/family_markers"
)

OUT_ROOT = (
    PROJECT
    / "06_phasing/families"
)

SUMMARY_DIR = (
    PROJECT
    / "06_phasing/summary"
)

OUT_ROOT.mkdir(
    parents=True,
    exist_ok=True
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. Read metadata
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str
)

families = sorted(
    metadata["family_id"].unique()
)

family_information = {}

for family in families:
    subset = metadata.loc[
        metadata["family_id"] == family
    ]

    queens = subset.loc[
        subset["role"] == "queen",
        "sample_id"
    ].tolist()

    offspring = (
        subset.loc[
            subset["role"] == "offspring",
            "sample_id"
        ]
        .sort_values()
        .tolist()
    )

    if len(queens) != 1:
        raise ValueError(
            f"{family}: expected one queen, "
            f"found {len(queens)}"
        )

    if len(offspring) != 5:
        raise ValueError(
            f"{family}: expected five offspring, "
            f"found {len(offspring)}"
        )

    family_information[family] = {
        "queen": queens[0],
        "offspring": offspring,
    }


# ============================================================
# 3. Chromosome order
# ============================================================

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

chromosome_rank = {
    chromosome: index
    for index, chromosome
    in enumerate(chromosome_order)
}

chromosome_length = dict(
    zip(
        primary_bed["chromosome"],
        primary_bed["end"]
        - primary_bed["start"]
    )
)


# ============================================================
# 4. Helper functions
# ============================================================

def classify_phase_support(
    n_comparable: int,
    n_majority: int,
) -> str:
    """
    Classify the support for one adjacent-marker
    phase relationship.
    """

    if n_comparable == 0:
        return "unresolved"

    if 2 * n_majority == n_comparable:
        return "tie"

    if n_majority == n_comparable:
        return "unanimous"

    support = n_majority / n_comparable

    if support >= 0.80:
        return "strong_majority"

    return "weak_majority"


def phase_one_chromosome(
    family: str,
    chromosome: str,
    chromosome_data: pd.DataFrame,
    offspring: list[str],
):
    """
    Reconstruct the two maternal haplotypes using
    adjacent-marker majority linkage.

    phase_bit = 0:
        H0 allele = REF
        H1 allele = ALT

    phase_bit = 1:
        H0 allele = ALT
        H1 allele = REF
    """

    chromosome_data = (
        chromosome_data
        .sort_values("position")
        .reset_index(drop=True)
        .copy()
    )

    n_markers = len(chromosome_data)

    genotype_matrix = (
        chromosome_data[offspring]
        .apply(
            pd.to_numeric,
            errors="coerce"
        )
        .to_numpy(dtype=float)
    )

    finite_values = genotype_matrix[
        np.isfinite(genotype_matrix)
    ]

    invalid_values = finite_values[
        ~np.isin(
            finite_values,
            [0, 1]
        )
    ]

    if len(invalid_values) > 0:
        raise ValueError(
            f"{family} {chromosome}: "
            "offspring marker states must be "
            "0, 1 or NA."
        )

    phase_bits = np.zeros(
        n_markers,
        dtype=np.int8
    )

    support_from_previous = np.full(
        n_markers,
        np.nan,
        dtype=float
    )

    category_from_previous = np.full(
        n_markers,
        "first_marker",
        dtype=object
    )

    relation_from_previous = np.full(
        n_markers,
        "first_marker",
        dtype=object
    )

    link_rows = []

    for marker_index in range(
        1,
        n_markers
    ):
        previous_calls = genotype_matrix[
            marker_index - 1
        ]

        current_calls = genotype_matrix[
            marker_index
        ]

        comparable = (
            np.isfinite(previous_calls)
            & np.isfinite(current_calls)
        )

        n_comparable = int(
            comparable.sum()
        )

        previous_position = int(
            chromosome_data.loc[
                marker_index - 1,
                "position"
            ]
        )

        current_position = int(
            chromosome_data.loc[
                marker_index,
                "position"
            ]
        )

        if n_comparable == 0:
            # No linkage information.
            # Carry the previous orientation forward,
            # but mark this relationship unresolved.
            relation_bit = 0
            n_same = 0
            n_opposite = 0
            n_majority = 0
            support = np.nan
            category = "unresolved"
            switching_offspring = []

        else:
            previous_binary = (
                previous_calls[comparable]
                .astype(np.int8)
            )

            current_binary = (
                current_calls[comparable]
                .astype(np.int8)
            )

            same_mask = (
                previous_binary
                == current_binary
            )

            n_same = int(
                same_mask.sum()
            )

            n_opposite = (
                n_comparable - n_same
            )

            if n_same > n_opposite:
                # The REF/ALT orientation remains
                # the same between adjacent markers.
                relation_bit = 0
                n_majority = n_same

            elif n_opposite > n_same:
                # Flip REF/ALT orientation at the
                # current marker.
                relation_bit = 1
                n_majority = n_opposite

            else:
                # A tie is unresolved. Carry the
                # previous orientation forward but
                # retain an explicit warning.
                relation_bit = 0
                n_majority = n_same

            support = (
                n_majority / n_comparable
            )

            category = classify_phase_support(
                n_comparable,
                n_majority
            )

        phase_bits[marker_index] = (
            phase_bits[marker_index - 1]
            ^ relation_bit
        )

        support_from_previous[
            marker_index
        ] = support

        category_from_previous[
            marker_index
        ] = category

        relation_from_previous[
            marker_index
        ] = (
            "same_orientation"
            if relation_bit == 0
            else "flipped_orientation"
        )

        # Determine which offspring change maternal
        # haplotype after applying the inferred phase.
        switching_offspring = []

        if n_comparable > 0:
            for offspring_index, child in enumerate(
                offspring
            ):
                previous_value = previous_calls[
                    offspring_index
                ]

                current_value = current_calls[
                    offspring_index
                ]

                if not (
                    np.isfinite(previous_value)
                    and np.isfinite(current_value)
                ):
                    continue

                previous_haplotype = (
                    int(previous_value)
                    ^ int(
                        phase_bits[
                            marker_index - 1
                        ]
                    )
                )

                current_haplotype = (
                    int(current_value)
                    ^ int(
                        phase_bits[
                            marker_index
                        ]
                    )
                )

                if (
                    previous_haplotype
                    != current_haplotype
                ):
                    switching_offspring.append(
                        child
                    )

        vote_margin = abs(
            n_same - n_opposite
        )

        link_rows.append(
            {
                "family_id": family,
                "chromosome": chromosome,
                "previous_marker_index": (
                    marker_index - 1
                ),
                "current_marker_index": (
                    marker_index
                ),
                "previous_position": (
                    previous_position
                ),
                "current_position": (
                    current_position
                ),
                "interval_bp": (
                    current_position
                    - previous_position
                ),
                "n_comparable_offspring": (
                    n_comparable
                ),
                "n_same_allele_encoding": (
                    n_same
                ),
                "n_opposite_allele_encoding": (
                    n_opposite
                ),
                "phase_relation_bit": (
                    relation_bit
                ),
                "chosen_phase_relation": (
                    relation_from_previous[
                        marker_index
                    ]
                ),
                "phase_support": support,
                "vote_margin": vote_margin,
                "phase_support_category": (
                    category
                ),
                "n_offspring_switching_haplotype": (
                    len(switching_offspring)
                ),
                "offspring_switching_haplotype": (
                    ",".join(
                        switching_offspring
                    )
                ),
                "weak_or_unresolved_link": (
                    category
                    in {
                        "weak_majority",
                        "tie",
                        "unresolved",
                    }
                ),
            }
        )

    # Convert offspring REF/ALT states into
    # maternal H0/H1 states.
    haplotype_matrix = np.full(
        genotype_matrix.shape,
        np.nan,
        dtype=float
    )

    expanded_phase = np.broadcast_to(
        phase_bits[:, np.newaxis],
        genotype_matrix.shape
    )

    valid_genotypes = np.isfinite(
        genotype_matrix
    )

    haplotype_matrix[
        valid_genotypes
    ] = (
        genotype_matrix[
            valid_genotypes
        ].astype(np.int8)
        ^ expanded_phase[
            valid_genotypes
        ].astype(np.int8)
    )

    chromosome_data[
        "phase_bit"
    ] = phase_bits

    chromosome_data[
        "maternal_haplotype_0_allele"
    ] = np.where(
        phase_bits == 0,
        chromosome_data[
            "reference_allele"
        ],
        chromosome_data[
            "alternative_allele"
        ],
    )

    chromosome_data[
        "maternal_haplotype_1_allele"
    ] = np.where(
        phase_bits == 0,
        chromosome_data[
            "alternative_allele"
        ],
        chromosome_data[
            "reference_allele"
        ],
    )

    chromosome_data[
        "phase_support_from_previous"
    ] = support_from_previous

    chromosome_data[
        "phase_category_from_previous"
    ] = category_from_previous

    chromosome_data[
        "phase_relation_from_previous"
    ] = relation_from_previous

    for offspring_index, child in enumerate(
        offspring
    ):
        chromosome_data[
            f"{child}_haplotype"
        ] = haplotype_matrix[
            :,
            offspring_index
        ]

    links = pd.DataFrame(
        link_rows
    )

    # Raw switch count per offspring.
    offspring_rows = []

    for offspring_index, child in enumerate(
        offspring
    ):
        child_states = haplotype_matrix[
            :,
            offspring_index
        ]

        previous_states = child_states[:-1]
        current_states = child_states[1:]

        comparable = (
            np.isfinite(previous_states)
            & np.isfinite(current_states)
        )

        n_comparable_links = int(
            comparable.sum()
        )

        n_raw_switches = int(
            (
                previous_states[comparable]
                != current_states[comparable]
            ).sum()
        )

        offspring_rows.append(
            {
                "family_id": family,
                "chromosome": chromosome,
                "sample_id": child,
                "n_markers": n_markers,
                "n_called_markers": int(
                    np.isfinite(
                        child_states
                    ).sum()
                ),
                "n_comparable_marker_links": (
                    n_comparable_links
                ),
                "n_raw_haplotype_switches": (
                    n_raw_switches
                ),
                "raw_switch_fraction": (
                    n_raw_switches
                    / n_comparable_links
                    if n_comparable_links > 0
                    else np.nan
                ),
            }
        )

    offspring_summary = pd.DataFrame(
        offspring_rows
    )

    return (
        chromosome_data,
        links,
        offspring_summary,
    )


# ============================================================
# 5. Run phasing family by family
# ============================================================

all_link_tables = []
all_offspring_summaries = []
chromosome_summary_rows = []

for family in families:
    print("=" * 60)
    print(f"Phasing family: {family}")

    offspring = family_information[
        family
    ]["offspring"]

    marker_file = (
        MARKER_ROOT
        / family
        / f"{family}.high_confidence_markers.tsv.gz"
    )

    if not marker_file.exists():
        raise FileNotFoundError(
            f"Marker file not found: "
            f"{marker_file}"
        )

    markers = pd.read_csv(
        marker_file,
        sep="\t"
    )

    required_columns = {
        "chromosome",
        "position",
        "reference_allele",
        "alternative_allele",
        "queen_id",
    } | set(offspring)

    missing_columns = (
        required_columns
        - set(markers.columns)
    )

    if missing_columns:
        raise ValueError(
            f"{family}: missing columns: "
            f"{sorted(missing_columns)}"
        )

    markers["chromosome"] = (
        markers["chromosome"]
        .astype(str)
    )

    markers["position"] = pd.to_numeric(
        markers["position"],
        errors="raise"
    ).astype(int)

    non_primary = sorted(
        set(markers["chromosome"])
        - set(chromosome_order)
    )

    if non_primary:
        raise ValueError(
            f"{family}: markers found on "
            f"non-primary contigs: {non_primary}"
        )

    duplicated = markers.duplicated(
        subset=[
            "chromosome",
            "position",
        ],
        keep=False
    )

    if duplicated.any():
        duplicate_table = markers.loc[
            duplicated,
            [
                "chromosome",
                "position",
            ]
        ]

        raise ValueError(
            f"{family}: duplicated marker "
            "positions were detected:\n"
            f"{duplicate_table.head(20)}"
        )

    markers["_chromosome_rank"] = (
        markers["chromosome"]
        .map(chromosome_rank)
    )

    markers = (
        markers
        .sort_values(
            [
                "_chromosome_rank",
                "position",
            ]
        )
        .drop(
            columns="_chromosome_rank"
        )
        .reset_index(drop=True)
    )

    phased_chromosomes = []
    family_links = []
    family_offspring_summary = []

    for chromosome in chromosome_order:
        chromosome_markers = markers.loc[
            markers["chromosome"]
            == chromosome
        ].copy()

        if chromosome_markers.empty:
            continue

        (
            phased_chromosome,
            chromosome_links,
            offspring_summary,
        ) = phase_one_chromosome(
            family=family,
            chromosome=chromosome,
            chromosome_data=(
                chromosome_markers
            ),
            offspring=offspring,
        )

        phased_chromosomes.append(
            phased_chromosome
        )

        family_links.append(
            chromosome_links
        )

        family_offspring_summary.append(
            offspring_summary
        )

        category_counts = Counter(
            chromosome_links[
                "phase_support_category"
            ]
            if not chromosome_links.empty
            else []
        )

        chromosome_summary_rows.append(
            {
                "family_id": family,
                "chromosome": chromosome,
                "chromosome_length": (
                    chromosome_length[
                        chromosome
                    ]
                ),
                "n_markers": len(
                    phased_chromosome
                ),
                "n_marker_links": len(
                    chromosome_links
                ),
                "n_unanimous_links": (
                    category_counts.get(
                        "unanimous",
                        0
                    )
                ),
                "n_strong_majority_links": (
                    category_counts.get(
                        "strong_majority",
                        0
                    )
                ),
                "n_weak_majority_links": (
                    category_counts.get(
                        "weak_majority",
                        0
                    )
                ),
                "n_tie_links": (
                    category_counts.get(
                        "tie",
                        0
                    )
                ),
                "n_unresolved_links": (
                    category_counts.get(
                        "unresolved",
                        0
                    )
                ),
                "fraction_weak_or_unresolved": (
                    chromosome_links[
                        "weak_or_unresolved_link"
                    ].mean()
                    if len(
                        chromosome_links
                    ) > 0
                    else np.nan
                ),
                "median_phase_support": (
                    chromosome_links[
                        "phase_support"
                    ].median()
                    if len(
                        chromosome_links
                    ) > 0
                    else np.nan
                ),
                "median_marker_interval_bp": (
                    chromosome_links[
                        "interval_bp"
                    ].median()
                    if len(
                        chromosome_links
                    ) > 0
                    else np.nan
                ),
                "maximum_marker_interval_bp": (
                    chromosome_links[
                        "interval_bp"
                    ].max()
                    if len(
                        chromosome_links
                    ) > 0
                    else np.nan
                ),
            }
        )

    phased_family = pd.concat(
        phased_chromosomes,
        ignore_index=True
    )

    links_family = pd.concat(
        family_links,
        ignore_index=True
    )

    offspring_family = pd.concat(
        family_offspring_summary,
        ignore_index=True
    )

    family_out = OUT_ROOT / family

    family_out.mkdir(
        parents=True,
        exist_ok=True
    )

    phased_family.to_csv(
        family_out
        / f"{family}.phased_markers.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    links_family.to_csv(
        family_out
        / f"{family}.adjacent_phase_links.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    offspring_family.to_csv(
        family_out
        / f"{family}.raw_switch_summary.tsv",
        sep="\t",
        index=False,
    )

    all_link_tables.append(
        links_family
    )

    all_offspring_summaries.append(
        offspring_family
    )

    n_weak = int(
        links_family[
            "weak_or_unresolved_link"
        ].sum()
    )

    print(
        f"Markers: {len(phased_family):,}"
    )

    print(
        f"Adjacent links: "
        f"{len(links_family):,}"
    )

    print(
        f"Weak/unresolved links: "
        f"{n_weak:,}"
    )


# ============================================================
# 6. Combined summaries
# ============================================================

all_links = pd.concat(
    all_link_tables,
    ignore_index=True
)

all_offspring_summary = pd.concat(
    all_offspring_summaries,
    ignore_index=True
)

chromosome_summary = pd.DataFrame(
    chromosome_summary_rows
)

all_links.to_csv(
    SUMMARY_DIR
    / "all_families_adjacent_phase_links.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

all_offspring_summary.to_csv(
    SUMMARY_DIR
    / "all_families_raw_switch_summary.tsv",
    sep="\t",
    index=False,
)

chromosome_summary.to_csv(
    SUMMARY_DIR
    / "phasing_chromosome_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 7. Family-level phasing summary
# ============================================================

family_summary_rows = []

for family in families:
    subset = all_links.loc[
        all_links["family_id"]
        == family
    ]

    category_counts = Counter(
        subset[
            "phase_support_category"
        ]
    )

    family_summary_rows.append(
        {
            "family_id": family,
            "n_adjacent_links": len(subset),
            "n_unanimous_links": (
                category_counts.get(
                    "unanimous",
                    0
                )
            ),
            "n_strong_majority_links": (
                category_counts.get(
                    "strong_majority",
                    0
                )
            ),
            "n_weak_majority_links": (
                category_counts.get(
                    "weak_majority",
                    0
                )
            ),
            "n_tie_links": (
                category_counts.get(
                    "tie",
                    0
                )
            ),
            "n_unresolved_links": (
                category_counts.get(
                    "unresolved",
                    0
                )
            ),
            "fraction_weak_or_unresolved": (
                subset[
                    "weak_or_unresolved_link"
                ].mean()
            ),
            "median_phase_support": (
                subset[
                    "phase_support"
                ].median()
            ),
        }
    )

family_summary = pd.DataFrame(
    family_summary_rows
)

family_summary.to_csv(
    SUMMARY_DIR
    / "phasing_family_summary.tsv",
    sep="\t",
    index=False,
)

print("\nFamily phasing summary:")
print(
    family_summary.to_string(
        index=False
    )
)

print(
    "\nPhasing completed."
)
