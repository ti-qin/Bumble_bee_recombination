from pathlib import Path
from collections import Counter, defaultdict
import gzip
import shutil
import subprocess

import numpy as np
import pandas as pd


# ============================================================
# 1. Paths and parameters
# ============================================================

PROJECT = Path(
    "/localdata/qinti/Project/Bee"
)

VCF_ROOT = (
    PROJECT
    / "03_variants/gatk_direct_diploid_SB/"
      "families"
)

METADATA = (
    PROJECT
    / "00_rawdata/metadata/"
      "analysis_samples.noC3.tsv"
)

MAPPING_SUMMARY = (
    PROJECT
    / "02_mapping/stats/"
      "mapping_summary.tsv"
)

PRIMARY_BED = (
    PROJECT
    / "00_rawdata/reference/"
      "primary_chromosomes.bed"
)

OUT_ROOT = (
    PROJECT
    / "05_markers/family_markers"
)

SUMMARY_DIR = (
    PROJECT
    / "05_markers/summary"
)

OUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ------------------------------------------------------------
# Marker-level genotype thresholds
# ------------------------------------------------------------

MIN_DP = 10
MIN_GQ = 30

MAX_DEPTH_MULTIPLIER = 2.5
DEFAULT_MAX_DP = 100

QUEEN_HET_MIN_AB = 0.25
QUEEN_HET_MAX_AB = 0.75
MIN_QUEEN_ALLELE_DEPTH = 5

OFFSPRING_HOM_REF_MAX_AB = 0.1
OFFSPRING_HOM_ALT_MIN_AB = 0.9

MIN_CALLED_OFFSPRING = 5

# Require more than three supporting reads from each relevant
# strand. Therefore, the minimum accepted count is four.
MIN_STRAND_READS_PER_DIRECTION = 2


# ============================================================
# 2. Locate bcftools
# ============================================================

BCFTOOLS = shutil.which(
    "bcftools"
)

if BCFTOOLS is None:
    raise FileNotFoundError(
        "bcftools was not found in PATH."
    )

print(
    f"Using bcftools: {BCFTOOLS}"
)


# ============================================================
# 3. Metadata and family structure
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
)

required_metadata_columns = {
    "sample_id",
    "family_id",
    "role",
}

missing_metadata_columns = (
    required_metadata_columns
    - set(metadata.columns)
)

if missing_metadata_columns:
    raise ValueError(
        "Metadata is missing columns: "
        f"{sorted(missing_metadata_columns)}"
    )

families = sorted(
    metadata["family_id"].unique()
)

expected_families = {
    "C1",
    "C2",
    "C4",
}

if set(families) != expected_families:
    raise ValueError(
        "Expected families C1, C2 and C4, "
        f"but found: {families}"
    )

family_samples = {}

for family in families:
    subset = metadata.loc[
        metadata["family_id"] == family
    ]

    queens = subset.loc[
        subset["role"] == "queen",
        "sample_id",
    ].tolist()

    offspring = (
        subset.loc[
            subset["role"] == "offspring",
            "sample_id",
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

    family_samples[family] = {
        "queen": queens[0],
        "offspring": offspring,
    }

print(
    "Families:",
    ", ".join(families),
)


# ============================================================
# 4. Sample-specific maximum depth
# ============================================================

sample_max_dp = {
    sample_id: DEFAULT_MAX_DP
    for sample_id
    in metadata["sample_id"]
}

mapping = pd.read_csv(
    MAPPING_SUMMARY,
    sep="\t",
)

for _, row in mapping.iterrows():
    sample_id = str(
        row["sample_id"]
    )

    if sample_id not in sample_max_dp:
        continue

    try:
        mean_depth = float(
            row["mean_depth"]
        )
    except (TypeError, ValueError):
        continue

    if (
        np.isfinite(mean_depth)
        and mean_depth > 0
    ):
        sample_max_dp[sample_id] = max(
            30,
            int(
                np.ceil(
                    mean_depth
                    * MAX_DEPTH_MULTIPLIER
                )
            ),
        )


# ============================================================
# 5. Chromosome lengths
# ============================================================

chromosomes = []
chromosome_lengths = {}

with PRIMARY_BED.open() as handle:
    for line in handle:
        if not line.strip():
            continue

        chrom, start, end = (
            line.rstrip("\n")
            .split("\t")[:3]
        )

        chromosomes.append(
            chrom
        )

        chromosome_lengths[chrom] = (
            int(end)
            - int(start)
        )


# ============================================================
# 6. Parsing helpers
# ============================================================

def get_family_vcf(
    family: str,
) -> Path:
    return (
        VCF_ROOT
        / family
        / "filtered"
        / (
            f"{family}.direct_diploid_SB."
            "biallelic.snps.PASS.vcf.gz"
        )
    )


def parse_sb(
    sb_text: str,
):
    """
    Parse FORMAT/SB:
    REF_forward, REF_reverse, ALT_forward, ALT_reverse.

    SB is used as a sample-level marker filter.

    For a heterozygous queen, REF and ALT must each have
    sufficient forward- and reverse-strand support.

    For a homozygous offspring, only the called major allele
    is required to have sufficient support on both strands.
    """

    if sb_text in {
        "",
        ".",
        None,
    }:
        return None

    parts = sb_text.split(
        ","
    )

    if len(parts) != 4:
        return None

    try:
        values = tuple(
            int(value)
            for value in parts
        )
    except ValueError:
        return None

    return {
        "ref_forward": values[0],
        "ref_reverse": values[1],
        "alt_forward": values[2],
        "alt_reverse": values[3],
    }


def queen_has_bidirectional_support(
    call: dict,
) -> bool:
    """
    The queen is heterozygous, so both REF and ALT must be
    supported by at least MIN_STRAND_READS_PER_DIRECTION
    reads on both the forward and reverse strands.
    """

    sb = call.get(
        "strand_bias_counts"
    )

    if sb is None:
        return False

    required_counts = [
        sb["ref_forward"],
        sb["ref_reverse"],
        sb["alt_forward"],
        sb["alt_reverse"],
    ]

    return (
        min(required_counts)
        >= MIN_STRAND_READS_PER_DIRECTION
    )


def offspring_has_bidirectional_support(
    call: dict,
    state: int,
) -> bool:
    """
    For a REF-homozygous offspring call (state=0), require
    REF support on both strands. For an ALT-homozygous call
    (state=1), require ALT support on both strands.

    The minor allele is not required to have strand support,
    because a clean haploid/homozygous call should contain
    little or no minor-allele evidence.
    """

    sb = call.get(
        "strand_bias_counts"
    )

    if sb is None:
        return False

    if state == 0:
        required_counts = [
            sb["ref_forward"],
            sb["ref_reverse"],
        ]
    elif state == 1:
        required_counts = [
            sb["alt_forward"],
            sb["alt_reverse"],
        ]
    else:
        raise ValueError(
            f"Unexpected offspring state: {state}"
        )

    return (
        min(required_counts)
        >= MIN_STRAND_READS_PER_DIRECTION
    )


def parse_common_fields(
    sample_id: str,
    field: str,
):
    parts = field.split(
        ":"
    )

    if len(parts) != 5:
        return None

    (
        gt,
        dp_text,
        gq_text,
        ad_text,
        sb_text,
    ) = parts

    if (
        gt in {
            ".",
            "./.",
            ".|.",
        }
        or dp_text == "."
        or gq_text == "."
        or ad_text == "."
    ):
        return None

    try:
        dp = int(
            dp_text
        )

        gq = float(
            gq_text
        )
    except ValueError:
        return None

    if dp < MIN_DP:
        return None

    if gq < MIN_GQ:
        return None

    if dp > sample_max_dp.get(
        sample_id,
        DEFAULT_MAX_DP,
    ):
        return None

    ad_parts = ad_text.split(
        ","
    )

    if len(ad_parts) != 2:
        return None

    try:
        ref_depth = int(
            ad_parts[0]
        )

        alt_depth = int(
            ad_parts[1]
        )
    except ValueError:
        return None

    total_ad = (
        ref_depth
        + alt_depth
    )

    if total_ad < MIN_DP:
        return None

    allele_balance = (
        alt_depth
        / total_ad
    )

    return {
        "gt": gt.replace(
            "|",
            "/",
        ),
        "dp": dp,
        "gq": gq,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "total_ad": total_ad,
        "allele_balance": allele_balance,
        "strand_bias_counts": parse_sb(
            sb_text
        ),
        "sb_raw": sb_text,
    }


def parse_queen_heterozygote(
    sample_id: str,
    field: str,
):
    call = parse_common_fields(
        sample_id,
        field,
    )

    if call is None:
        return None

    if call["gt"] not in {
        "0/1",
        "1/0",
    }:
        return None

    if min(
        call["ref_depth"],
        call["alt_depth"],
    ) < MIN_QUEEN_ALLELE_DEPTH:
        return None

    if not (
        QUEEN_HET_MIN_AB
        <= call["allele_balance"]
        <= QUEEN_HET_MAX_AB
    ):
        return None

    if not queen_has_bidirectional_support(
        call
    ):
        return None

    return call


def parse_offspring_allele(
    sample_id: str,
    field: str,
):
    call = parse_common_fields(
        sample_id,
        field,
    )

    if call is None:
        return None

    gt = call["gt"]
    ab = call["allele_balance"]

    if gt in {
        "0",
        "0/0",
    }:
        if (
            ab
            <= OFFSPRING_HOM_REF_MAX_AB
            and offspring_has_bidirectional_support(
                call,
                state=0,
            )
        ):
            return {
                "state": 0,
                "call": call,
            }

    if gt in {
        "1",
        "1/1",
    }:
        if (
            ab
            >= OFFSPRING_HOM_ALT_MIN_AB
            and offspring_has_bidirectional_support(
                call,
                state=1,
            )
        ):
            return {
                "state": 1,
                "call": call,
            }

    return None


# ============================================================
# 7. Process one family VCF at a time
# ============================================================

marker_counts = defaultdict(
    int
)

marker_gaps = defaultdict(
    list
)

filter_summary_rows = []
offspring_qc_rows = []

for family in families:
    vcf = get_family_vcf(
        family
    )

    if not vcf.exists():
        raise FileNotFoundError(
            f"{family} PASS VCF not found: {vcf}"
        )

    queen = family_samples[
        family
    ]["queen"]

    offspring = family_samples[
        family
    ]["offspring"]

    expected_samples = [
        queen,
        *offspring,
    ]

    sample_result = subprocess.run(
        [
            BCFTOOLS,
            "query",
            "-l",
            str(vcf),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if sample_result.returncode != 0:
        raise RuntimeError(
            sample_result.stderr
        )

    vcf_samples = [
        line.strip()
        for line
        in sample_result.stdout.splitlines()
        if line.strip()
    ]

    if set(vcf_samples) != set(
        expected_samples
    ):
        raise ValueError(
            f"{family}: VCF samples differ from "
            "family metadata.\n"
            f"Expected: {sorted(expected_samples)}\n"
            f"Observed: {sorted(vcf_samples)}"
        )

    sample_index = {
        sample_id: index
        for index, sample_id
        in enumerate(vcf_samples)
    }

    family_dir = (
        OUT_ROOT
        / family
    )

    family_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    marker_file = (
        family_dir
        / (
            f"{family}."
            "high_confidence_markers.tsv.gz"
        )
    )

    bed_file = (
        family_dir
        / (
            f"{family}."
            "high_confidence_markers.bed"
        )
    )

    last_position = {}

    family_counters = Counter()

    child_callable_at_queen_sites = Counter()
    child_failed_at_queen_sites = Counter()

    header = [
        "chromosome",
        "position",
        "reference_allele",
        "alternative_allele",
        "queen_id",
        "queen_gt",
        "queen_dp",
        "queen_gq",
        "queen_ref_depth",
        "queen_alt_depth",
        "queen_min_allele_depth",
        "queen_alt_fraction",
        "n_called_offspring",
        "offspring_call_rate",
    ] + offspring

    query_format = (
        "%CHROM\\t%POS\\t%REF\\t%ALT"
        "[\\t%GT:%DP:%GQ:%AD:%SB]"
        "\\n"
    )

    process = subprocess.Popen(
        [
            BCFTOOLS,
            "query",
            "-f",
            query_format,
            str(vcf),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    if process.stdout is None:
        raise RuntimeError(
            f"Failed to query {vcf}"
        )

    with (
        gzip.open(
            marker_file,
            "wt",
        ) as output_handle,
        bed_file.open(
            "w",
        ) as bed_handle,
    ):
        output_handle.write(
            "\t".join(header)
            + "\n"
        )

        for line in process.stdout:
            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if (
                len(fields)
                != 4 + len(vcf_samples)
            ):
                family_counters[
                    "malformed_query_rows"
                ] += 1

                continue

            chromosome = fields[0]
            position = int(
                fields[1]
            )
            ref = fields[2]
            alt = fields[3]

            genotype_fields = fields[
                4:
            ]

            family_counters[
                "pass_biallelic_snp_sites_examined"
            ] += 1

            queen_field = genotype_fields[
                sample_index[queen]
            ]

            queen_call = (
                parse_queen_heterozygote(
                    queen,
                    queen_field,
                )
            )

            if queen_call is None:
                family_counters[
                    "queen_filter_failed"
                ] += 1

                continue

            family_counters[
                "queen_filter_passed"
            ] += 1

            child_results = []

            for child in offspring:
                child_field = genotype_fields[
                    sample_index[child]
                ]

                child_result = (
                    parse_offspring_allele(
                        child,
                        child_field,
                    )
                )

                child_results.append(
                    child_result
                )

                if child_result is None:
                    child_failed_at_queen_sites[
                        child
                    ] += 1
                else:
                    child_callable_at_queen_sites[
                        child
                    ] += 1

            n_called = sum(
                result is not None
                for result in child_results
            )

            if (
                n_called
                != MIN_CALLED_OFFSPRING
            ):
                family_counters[
                    "offspring_5of5_filter_failed"
                ] += 1

                continue

            family_counters[
                "offspring_5of5_filter_passed"
            ] += 1

            child_states = [
                result["state"]
                for result
                in child_results
            ]

            call_rate = (
                n_called
                / len(offspring)
            )

            row = [
                chromosome,
                str(position),
                ref,
                alt,
                queen,
                queen_call["gt"],
                str(
                    queen_call["dp"]
                ),
                (
                    f"{queen_call['gq']:.2f}"
                ),
                str(
                    queen_call[
                        "ref_depth"
                    ]
                ),
                str(
                    queen_call[
                        "alt_depth"
                    ]
                ),
                str(
                    min(
                        queen_call[
                            "ref_depth"
                        ],
                        queen_call[
                            "alt_depth"
                        ],
                    )
                ),
                (
                    f"{queen_call['allele_balance']:.4f}"
                ),
                str(n_called),
                f"{call_rate:.4f}",
            ]

            row.extend(
                str(state)
                for state
                in child_states
            )

            output_handle.write(
                "\t".join(row)
                + "\n"
            )

            # BED: 0-based, half-open.
            bed_handle.write(
                f"{chromosome}\t"
                f"{position - 1}\t"
                f"{position}\n"
            )

            marker_counts[
                (
                    family,
                    chromosome,
                )
            ] += 1

            family_counters[
                "high_confidence_markers"
            ] += 1

            previous = last_position.get(
                chromosome
            )

            if previous is not None:
                marker_gaps[
                    (
                        family,
                        chromosome,
                    )
                ].append(
                    position
                    - previous
                )

            last_position[
                chromosome
            ] = position

    return_code = process.wait()

    stderr = ""

    if process.stderr is not None:
        stderr = (
            process.stderr.read()
        )

    if return_code != 0:
        raise RuntimeError(
            f"bcftools query failed for "
            f"{family}:\n{stderr}"
        )

    for key, value in sorted(
        family_counters.items()
    ):
        filter_summary_rows.append(
            {
                "family_id": family,
                "stage": key,
                "n_sites": value,
            }
        )

    for child in offspring:
        n_callable = (
            child_callable_at_queen_sites[
                child
            ]
        )

        n_failed = (
            child_failed_at_queen_sites[
                child
            ]
        )

        total = (
            n_callable
            + n_failed
        )

        offspring_qc_rows.append(
            {
                "family_id": family,
                "sample_id": child,
                "n_queen_pass_sites": total,
                "n_callable_sites": n_callable,
                "n_failed_sites": n_failed,
                "call_rate_at_queen_pass_sites": (
                    n_callable / total
                    if total > 0
                    else np.nan
                ),
                "n_final_family_markers": (
                    family_counters[
                        "high_confidence_markers"
                    ]
                ),
            }
        )

        filter_summary_rows.append(
            {
                "family_id": family,
                "stage": (
                    f"{child}:"
                    "callable_at_queen_pass_sites"
                ),
                "n_sites": n_callable,
            }
        )

        filter_summary_rows.append(
            {
                "family_id": family,
                "stage": (
                    f"{child}:"
                    "failed_at_queen_pass_sites"
                ),
                "n_sites": n_failed,
            }
        )

    print(
        f"{family}: "
        f"{family_counters['high_confidence_markers']:,} "
        "high-confidence markers"
    )


# ============================================================
# 8. Chromosome-level marker summary
# ============================================================

summary_rows = []

for family in families:
    for chromosome in chromosomes:
        n_markers = marker_counts[
            (
                family,
                chromosome,
            )
        ]

        chromosome_length = (
            chromosome_lengths[
                chromosome
            ]
        )

        gaps = marker_gaps.get(
            (
                family,
                chromosome,
            ),
            [],
        )

        if gaps:
            median_gap = float(
                np.median(gaps)
            )

            mean_gap = float(
                np.mean(gaps)
            )

            p95_gap = float(
                np.percentile(
                    gaps,
                    95,
                )
            )

            max_gap = int(
                np.max(gaps)
            )
        else:
            median_gap = np.nan
            mean_gap = np.nan
            p95_gap = np.nan
            max_gap = np.nan

        summary_rows.append(
            {
                "family_id": family,
                "chromosome": chromosome,
                "chromosome_length": (
                    chromosome_length
                ),
                "n_markers": n_markers,
                "markers_per_Mb": (
                    n_markers
                    / chromosome_length
                    * 1_000_000
                ),
                "median_marker_gap_bp": (
                    median_gap
                ),
                "mean_marker_gap_bp": (
                    mean_gap
                ),
                "p95_marker_gap_bp": (
                    p95_gap
                ),
                "max_marker_gap_bp": (
                    max_gap
                ),
            }
        )

marker_summary = pd.DataFrame(
    summary_rows
)

marker_summary.to_csv(
    SUMMARY_DIR
    / "family_marker_chromosome_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 9. Family summary
# ============================================================

family_rows = []

for family in families:
    n_markers = sum(
        marker_counts[
            (
                family,
                chromosome,
            )
        ]
        for chromosome
        in chromosomes
    )

    family_rows.append(
        {
            "family_id": family,
            "queen_id": (
                family_samples[
                    family
                ]["queen"]
            ),
            "n_offspring": len(
                family_samples[
                    family
                ]["offspring"]
            ),
            "n_high_confidence_markers": (
                n_markers
            ),
        }
    )

family_summary = pd.DataFrame(
    family_rows
)

family_summary.to_csv(
    SUMMARY_DIR
    / "family_marker_summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 10. Marker-filter summary and parameters
# ============================================================

filter_summary = pd.DataFrame(
    filter_summary_rows
)

filter_summary.to_csv(
    SUMMARY_DIR
    / "family_marker_filter_summary.tsv",
    sep="\t",
    index=False,
)

offspring_qc = pd.DataFrame(
    offspring_qc_rows
)

offspring_qc.to_csv(
    SUMMARY_DIR
    / "offspring_marker_missingness.tsv",
    sep="\t",
    index=False,
)

parameter_rows = [
    {
        "parameter": "MIN_DP",
        "value": MIN_DP,
    },
    {
        "parameter": "MIN_GQ",
        "value": MIN_GQ,
    },
    {
        "parameter": "MAX_DEPTH_MULTIPLIER",
        "value": MAX_DEPTH_MULTIPLIER,
    },
    {
        "parameter": "QUEEN_HET_MIN_AB",
        "value": QUEEN_HET_MIN_AB,
    },
    {
        "parameter": "QUEEN_HET_MAX_AB",
        "value": QUEEN_HET_MAX_AB,
    },
    {
        "parameter": "MIN_QUEEN_ALLELE_DEPTH",
        "value": MIN_QUEEN_ALLELE_DEPTH,
    },
    {
        "parameter": "OFFSPRING_HOM_REF_MAX_AB",
        "value": OFFSPRING_HOM_REF_MAX_AB,
    },
    {
        "parameter": "OFFSPRING_HOM_ALT_MIN_AB",
        "value": OFFSPRING_HOM_ALT_MIN_AB,
    },
    {
        "parameter": "MIN_CALLED_OFFSPRING",
        "value": MIN_CALLED_OFFSPRING,
    },
    {
        "parameter": "MIN_STRAND_READS_PER_DIRECTION",
        "value": MIN_STRAND_READS_PER_DIRECTION,
    },
    {
        "parameter": "FORMAT_SB",
        "value": (
            "hard_filter_queen_both_alleles_"
            "offspring_called_allele"
        ),
    },
]

pd.DataFrame(
    parameter_rows
).to_csv(
    SUMMARY_DIR
    / "family_marker_parameters.tsv",
    sep="\t",
    index=False,
)

print(
    "\nFamily marker summary:"
)

print(
    family_summary.to_string(
        index=False
    )
)

print(
    "\nSaved marker and QC tables under:"
)

print(
    OUT_ROOT
)

print(
    SUMMARY_DIR
)
