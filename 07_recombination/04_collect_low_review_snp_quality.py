from pathlib import Path
from collections import defaultdict
import shutil
import subprocess

import numpy as np
import pandas as pd


# ============================================================
# 1. Paths
# ============================================================

PROJECT = Path(
    "/localdata/qinti/Project/Bee"
)

VCF = (
    PROJECT
    / "03_variants/gatk_final/filtered/"
      "C1_C2_C4.biallelic.snps.PASS.vcf.gz"
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

CONTEXT_FILE = (
    EVENT_ROOT
    / "manual_review/"
      "CO_NCO_event_marker_context.tsv.gz"
)

OUT_DIR = (
    EVENT_ROOT
    / "manual_review/"
      "LOW_MANUAL_REVIEW_SNP_QC"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. Current marker thresholds
# ============================================================

MIN_DP = 10
MIN_GQ = 30

MAX_DEPTH_MULTIPLIER = 2.5
DEFAULT_MAX_DP = 100

QUEEN_MIN_AB = 0.25
QUEEN_MAX_AB = 0.75
MIN_QUEEN_ALLELE_DEPTH = 5

OFFSPRING_REF_MAX_AB = 0.10
OFFSPRING_ALT_MIN_AB = 0.90


# ============================================================
# 3. Locate bcftools
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
# 4. Metadata
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
)

family_samples = {}

for family, subset in metadata.groupby(
    "family_id",
    sort=False,
):
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
            f"{family}: expected one queen."
        )

    if len(offspring) != 5:
        raise ValueError(
            f"{family}: expected five offspring."
        )

    family_samples[family] = {
        "queen": queens[0],
        "offspring": offspring,
    }


# ============================================================
# 5. Sample-specific maximum depth
# ============================================================

sample_max_dp = {
    sample_id: DEFAULT_MAX_DP
    for sample_id in metadata[
        "sample_id"
    ]
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
# 6. Chromosome order
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
    dtype={
        "chromosome": str,
    },
)

chromosome_rank = {
    chromosome: index
    for index, chromosome
    in enumerate(
        primary["chromosome"]
    )
}


# ============================================================
# 7. Select LOW_MANUAL_REVIEW events
# ============================================================

events = pd.read_csv(
    EVENT_FILE,
    sep="\t",
)

low_events = events.loc[
    events["review_status"]
    == "LOW_MANUAL_REVIEW"
].copy()

if low_events.empty:
    raise RuntimeError(
        "No LOW_MANUAL_REVIEW events found."
    )

print(
    "LOW_MANUAL_REVIEW events:",
    len(low_events),
)

event_columns = [
    "event_id",
    "family_id",
    "sample_id",
    "chromosome",
    "event_type",
    "event_context",
    "event_interval_start",
    "event_interval_end",
    "core_start",
    "core_end",
    "minimum_phase_support",
    "raw_flag_reason_codes",
    "effective_low_reason_codes",
    "review_status",
]

event_columns = [
    column
    for column in event_columns
    if column in low_events.columns
]

low_event_info = low_events[
    event_columns
].copy()


# ============================================================
# 8. Obtain marker context
# ============================================================

context = pd.read_csv(
    CONTEXT_FILE,
    sep="\t",
)

low_context = context.loc[
    context["event_id"].isin(
        low_events["event_id"]
    )
].copy()

if low_context.empty:
    raise RuntimeError(
        "No marker context was found for "
        "LOW_MANUAL_REVIEW events."
    )

# Keep one record per event-marker.
low_context = (
    low_context.sort_values(
        [
            "event_id",
            "position",
        ]
    )
    .drop_duplicates(
        [
            "event_id",
            "chromosome",
            "position",
        ]
    )
    .reset_index(drop=True)
)

# Add event-level information.
extra_event_columns = [
    column
    for column in low_event_info.columns
    if column not in {
        "event_id",
        "family_id",
        "sample_id",
        "chromosome",
        "event_type",
        "event_context",
        "review_status",
        "effective_low_reason_codes",
    }
]

low_context = low_context.merge(
    low_event_info[
        ["event_id"] + extra_event_columns
    ],
    on="event_id",
    how="left",
    suffixes=("", "_event"),
)


# ============================================================
# 9. Make a unique marker BED
# ============================================================

unique_positions = (
    low_context[
        [
            "chromosome",
            "position",
        ]
    ]
    .drop_duplicates()
    .copy()
)

unique_positions[
    "_chromosome_rank"
] = (
    unique_positions["chromosome"]
    .map(chromosome_rank)
    .fillna(
        len(chromosome_rank)
    )
)

unique_positions = (
    unique_positions.sort_values(
        [
            "_chromosome_rank",
            "position",
        ]
    )
    .drop(
        columns="_chromosome_rank"
    )
)

region_bed = (
    OUT_DIR
    / "LOW_MANUAL_REVIEW_marker_positions.bed"
)

bed = pd.DataFrame(
    {
        "chromosome": (
            unique_positions[
                "chromosome"
            ]
        ),
        "start": (
            unique_positions[
                "position"
            ].astype(int) - 1
        ),
        "end": (
            unique_positions[
                "position"
            ].astype(int)
        ),
    }
)

bed.to_csv(
    region_bed,
    sep="\t",
    header=False,
    index=False,
)

print(
    "Unique marker positions:",
    len(unique_positions),
)


# ============================================================
# 10. VCF sample order
# ============================================================

sample_result = subprocess.run(
    [
        BCFTOOLS,
        "query",
        "-l",
        str(VCF),
    ],
    capture_output=True,
    text=True,
    check=True,
)

vcf_samples = [
    line.strip()
    for line
    in sample_result.stdout.splitlines()
    if line.strip()
]

analysis_samples = (
    metadata["sample_id"]
    .tolist()
)

missing_samples = (
    set(analysis_samples)
    - set(vcf_samples)
)

if missing_samples:
    raise ValueError(
        "Samples missing from VCF: "
        f"{sorted(missing_samples)}"
    )

# Query samples in a defined order.
query_samples = [
    sample_id
    for sample_id in vcf_samples
    if sample_id in set(
        analysis_samples
    )
]


# ============================================================
# 11. Query site and genotype metrics
# ============================================================

query_format = (
    "%CHROM\\t"
    "%POS\\t"
    "%REF\\t"
    "%ALT\\t"
    "%QUAL\\t"
    "%FILTER\\t"
    "%INFO/QD\\t"
    "%INFO/FS\\t"
    "%INFO/SOR\\t"
    "%INFO/MQ\\t"
    "%INFO/MQRankSum\\t"
    "%INFO/ReadPosRankSum"
    "[\\t%GT:%DP:%GQ:%AD]"
    "\\n"
)

query_command = [
    BCFTOOLS,
    "query",
    "-R",
    str(region_bed),
    "-s",
    ",".join(query_samples),
    "-f",
    query_format,
    str(VCF),
]

process = subprocess.run(
    query_command,
    capture_output=True,
    text=True,
    check=False,
)

if process.returncode != 0:
    raise RuntimeError(
        process.stderr
    )


# ============================================================
# 12. Parsing helpers
# ============================================================

def parse_float(value):
    if value in {
        None,
        "",
        ".",
        "nan",
    }:
        return np.nan

    try:
        return float(value)
    except ValueError:
        return np.nan


def parse_int(value):
    if value in {
        None,
        "",
        ".",
    }:
        return np.nan

    try:
        return int(value)
    except ValueError:
        return np.nan


def parse_sample_field(field):
    parts = str(field).split(":")

    if len(parts) != 4:
        return {
            "gt": ".",
            "dp": np.nan,
            "gq": np.nan,
            "ref_depth": np.nan,
            "alt_depth": np.nan,
            "total_ad": np.nan,
            "alt_fraction": np.nan,
        }

    gt, dp_text, gq_text, ad_text = (
        parts
    )

    gt = gt.replace("|", "/")

    dp = parse_int(dp_text)
    gq = parse_float(gq_text)

    ref_depth = np.nan
    alt_depth = np.nan

    if ad_text not in {
        ".",
        "",
    }:
        ad_parts = ad_text.split(",")

        if len(ad_parts) >= 2:
            ref_depth = parse_int(
                ad_parts[0]
            )

            alt_depth = parse_int(
                ad_parts[1]
            )

    if (
        pd.notna(ref_depth)
        and pd.notna(alt_depth)
    ):
        total_ad = (
            ref_depth
            + alt_depth
        )
    else:
        total_ad = np.nan

    if (
        pd.notna(total_ad)
        and total_ad > 0
    ):
        alt_fraction = (
            alt_depth / total_ad
        )
    else:
        alt_fraction = np.nan

    return {
        "gt": gt,
        "dp": dp,
        "gq": gq,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "total_ad": total_ad,
        "alt_fraction": alt_fraction,
    }


# ============================================================
# 13. Store queried VCF records
# ============================================================

site_records = {}

for line in process.stdout.splitlines():
    fields = line.rstrip("\n").split(
        "\t"
    )

    expected_fields = (
        12 + len(query_samples)
    )

    if len(fields) != expected_fields:
        continue

    chromosome = fields[0]
    position = int(fields[1])

    site_records[
        (
            chromosome,
            position,
        )
    ] = {
        "chromosome": chromosome,
        "position": position,
        "reference_allele_vcf": (
            fields[2]
        ),
        "alternative_allele_vcf": (
            fields[3]
        ),
        "QUAL": parse_float(
            fields[4]
        ),
        "FILTER": fields[5],
        "QD": parse_float(
            fields[6]
        ),
        "FS": parse_float(
            fields[7]
        ),
        "SOR": parse_float(
            fields[8]
        ),
        "MQ": parse_float(
            fields[9]
        ),
        "MQRankSum": parse_float(
            fields[10]
        ),
        "ReadPosRankSum": parse_float(
            fields[11]
        ),
        "sample_fields": {
            sample_id: field
            for sample_id, field
            in zip(
                query_samples,
                fields[12:],
            )
        },
    }


# ============================================================
# 14. Genotype QC assessment
# ============================================================

def assess_sample_qc(
    sample_id,
    role,
    call,
):
    fail_reasons = []
    warning_reasons = []

    gt = call["gt"]
    dp = call["dp"]
    gq = call["gq"]
    ref_depth = call["ref_depth"]
    alt_depth = call["alt_depth"]
    total_ad = call["total_ad"]
    ab = call["alt_fraction"]

    if gt in {
        ".",
        "./.",
        "",
    }:
        fail_reasons.append(
            "missing_GT"
        )

    if pd.isna(dp):
        fail_reasons.append(
            "missing_DP"
        )
    elif dp < MIN_DP:
        fail_reasons.append(
            "low_DP"
        )
    elif dp > sample_max_dp.get(
        sample_id,
        DEFAULT_MAX_DP,
    ):
        fail_reasons.append(
            "high_DP"
        )
    elif dp <= 12:
        warning_reasons.append(
            "DP_near_minimum"
        )

    if pd.isna(gq):
        fail_reasons.append(
            "missing_GQ"
        )
    elif gq < MIN_GQ:
        fail_reasons.append(
            "low_GQ"
        )
    elif gq < 40:
        warning_reasons.append(
            "GQ_30_to_39"
        )

    if pd.isna(total_ad):
        fail_reasons.append(
            "missing_AD"
        )
    elif total_ad < MIN_DP:
        fail_reasons.append(
            "low_total_AD"
        )

    inferred_state = np.nan
    major_allele_fraction = np.nan

    if role == "queen":
        if gt not in {
            "0/1",
            "1/0",
        }:
            fail_reasons.append(
                "queen_not_heterozygous"
            )

        if (
            pd.notna(ref_depth)
            and pd.notna(alt_depth)
        ):
            minimum_allele_depth = min(
                ref_depth,
                alt_depth,
            )

            if (
                minimum_allele_depth
                < MIN_QUEEN_ALLELE_DEPTH
            ):
                fail_reasons.append(
                    "queen_min_AD_below_5"
                )
            elif minimum_allele_depth <= 6:
                warning_reasons.append(
                    "queen_min_AD_5_to_6"
                )
        else:
            minimum_allele_depth = (
                np.nan
            )

        if pd.isna(ab):
            fail_reasons.append(
                "queen_missing_AB"
            )
        elif not (
            QUEEN_MIN_AB
            <= ab
            <= QUEEN_MAX_AB
        ):
            fail_reasons.append(
                "queen_AB_outside_range"
            )
        elif (
            ab <= 0.30
            or ab >= 0.70
        ):
            warning_reasons.append(
                "queen_AB_near_edge"
            )

    else:
        minimum_allele_depth = (
            np.nan
        )

        if (
            gt in {"0", "0/0"}
            and pd.notna(ab)
            and ab <= OFFSPRING_REF_MAX_AB
        ):
            inferred_state = 0
            major_allele_fraction = (
                1 - ab
            )

        elif (
            gt in {"1", "1/1"}
            and pd.notna(ab)
            and ab >= OFFSPRING_ALT_MIN_AB
        ):
            inferred_state = 1
            major_allele_fraction = ab

        else:
            fail_reasons.append(
                "offspring_ambiguous_allele"
            )

        if (
            pd.notna(
                major_allele_fraction
            )
            and
            major_allele_fraction < 0.95
        ):
            warning_reasons.append(
                "offspring_AF_0.90_to_0.95"
            )

    qc_pass = (
        len(fail_reasons) == 0
    )

    return {
        "sample_qc_pass": qc_pass,
        "sample_qc_fail_reasons": (
            ";".join(fail_reasons)
            if fail_reasons
            else "."
        ),
        "sample_qc_warning_reasons": (
            ";".join(warning_reasons)
            if warning_reasons
            else "."
        ),
        "inferred_offspring_state": (
            inferred_state
        ),
        "major_allele_fraction": (
            major_allele_fraction
        ),
        "minimum_queen_allele_depth": (
            minimum_allele_depth
        ),
    }


# ============================================================
# 15. Build event-marker-sample long table
# ============================================================

long_rows = []

for _, marker in low_context.iterrows():
    family = marker["family_id"]

    if family not in family_samples:
        raise ValueError(
            f"Unknown family: {family}"
        )

    queen = family_samples[
        family
    ]["queen"]

    offspring = family_samples[
        family
    ]["offspring"]

    samples = [
        queen,
        *offspring,
    ]

    key = (
        str(marker["chromosome"]),
        int(marker["position"]),
    )

    site = site_records.get(key)

    if site is None:
        for sample_id in samples:
            role = (
                "queen"
                if sample_id == queen
                else "offspring"
            )

            long_rows.append(
                {
                    **marker.to_dict(),
                    "sample_id_qc": (
                        sample_id
                    ),
                    "sample_role": role,
                    "site_found_in_vcf": (
                        False
                    ),
                    "FILTER": ".",
                    "QUAL": np.nan,
                    "QD": np.nan,
                    "FS": np.nan,
                    "SOR": np.nan,
                    "MQ": np.nan,
                    "MQRankSum": np.nan,
                    "ReadPosRankSum": (
                        np.nan
                    ),
                    "GT": ".",
                    "DP": np.nan,
                    "GQ": np.nan,
                    "AD_REF": np.nan,
                    "AD_ALT": np.nan,
                    "AD_TOTAL": np.nan,
                    "ALT_fraction": np.nan,
                    "major_allele_fraction": (
                        np.nan
                    ),
                    "minimum_queen_allele_depth": (
                        np.nan
                    ),
                    "sample_qc_pass": False,
                    "sample_qc_fail_reasons": (
                        "site_not_found_in_vcf"
                    ),
                    "sample_qc_warning_reasons": (
                        "."
                    ),
                }
            )

        continue

    for sample_id in samples:
        role = (
            "queen"
            if sample_id == queen
            else "offspring"
        )

        sample_field = site[
            "sample_fields"
        ].get(
            sample_id,
            ".:.:.:.",
        )

        call = parse_sample_field(
            sample_field
        )

        assessment = assess_sample_qc(
            sample_id,
            role,
            call,
        )

        long_rows.append(
            {
                **marker.to_dict(),
                "sample_id_qc": sample_id,
                "sample_role": role,
                "site_found_in_vcf": True,
                "reference_allele_vcf": (
                    site[
                        "reference_allele_vcf"
                    ]
                ),
                "alternative_allele_vcf": (
                    site[
                        "alternative_allele_vcf"
                    ]
                ),
                "FILTER": site["FILTER"],
                "QUAL": site["QUAL"],
                "QD": site["QD"],
                "FS": site["FS"],
                "SOR": site["SOR"],
                "MQ": site["MQ"],
                "MQRankSum": (
                    site["MQRankSum"]
                ),
                "ReadPosRankSum": (
                    site[
                        "ReadPosRankSum"
                    ]
                ),
                "GT": call["gt"],
                "DP": call["dp"],
                "GQ": call["gq"],
                "AD_REF": (
                    call["ref_depth"]
                ),
                "AD_ALT": (
                    call["alt_depth"]
                ),
                "AD_TOTAL": (
                    call["total_ad"]
                ),
                "ALT_fraction": (
                    call["alt_fraction"]
                ),
                **assessment,
            }
        )


long_qc = pd.DataFrame(
    long_rows
)

long_qc.to_csv(
    OUT_DIR
    / "LOW_MANUAL_REVIEW_marker_genotype_QC.long.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


# ============================================================
# 16. One-row-per-event-marker summary
# ============================================================

marker_summary_rows = []

group_columns = [
    "event_id",
    "family_id",
    "chromosome",
    "position",
]

for keys, group in long_qc.groupby(
    group_columns,
    sort=False,
):
    (
        event_id,
        family,
        chromosome,
        position,
    ) = keys

    first = group.iloc[0]

    failed = group.loc[
        ~group["sample_qc_pass"]
    ]

    warned = group.loc[
        group[
            "sample_qc_warning_reasons"
        ] != "."
    ]

    queen_rows = group.loc[
        group["sample_role"] == "queen"
    ]

    offspring_rows = group.loc[
        group[
            "sample_role"
        ] == "offspring"
    ]

    queen_row = (
        queen_rows.iloc[0]
        if not queen_rows.empty
        else None
    )

    marker_summary_rows.append(
        {
            "event_id": event_id,
            "family_id": family,
            "chromosome": chromosome,
            "position": int(position),
            "event_type": first.get(
                "event_type",
                "",
            ),
            "event_context": first.get(
                "event_context",
                "",
            ),
            "context_role": first.get(
                "context_role",
                "",
            ),
            "effective_low_reason_codes": (
                first.get(
                    "effective_low_reason_codes",
                    "",
                )
            ),
            "phase_support_from_previous": (
                first.get(
                    "phase_support_from_previous",
                    np.nan,
                )
            ),
            "phase_category_from_previous": (
                first.get(
                    "phase_category_from_previous",
                    "",
                )
            ),
            "FILTER": first["FILTER"],
            "QUAL": first["QUAL"],
            "QD": first["QD"],
            "FS": first["FS"],
            "SOR": first["SOR"],
            "MQ": first["MQ"],
            "MQRankSum": (
                first["MQRankSum"]
            ),
            "ReadPosRankSum": (
                first["ReadPosRankSum"]
            ),
            "queen_GT": (
                queen_row["GT"]
                if queen_row is not None
                else "."
            ),
            "queen_DP": (
                queen_row["DP"]
                if queen_row is not None
                else np.nan
            ),
            "queen_GQ": (
                queen_row["GQ"]
                if queen_row is not None
                else np.nan
            ),
            "queen_AD_REF": (
                queen_row["AD_REF"]
                if queen_row is not None
                else np.nan
            ),
            "queen_AD_ALT": (
                queen_row["AD_ALT"]
                if queen_row is not None
                else np.nan
            ),
            "queen_AB": (
                queen_row["ALT_fraction"]
                if queen_row is not None
                else np.nan
            ),
            "queen_min_allele_depth": (
                queen_row[
                    "minimum_queen_allele_depth"
                ]
                if queen_row is not None
                else np.nan
            ),
            "minimum_offspring_DP": (
                offspring_rows["DP"].min()
            ),
            "minimum_offspring_GQ": (
                offspring_rows["GQ"].min()
            ),
            "minimum_offspring_major_AF": (
                offspring_rows[
                    "major_allele_fraction"
                ].min()
            ),
            "n_family_samples_failed": (
                len(failed)
            ),
            "failed_samples": ",".join(
                failed["sample_id_qc"]
                .astype(str)
                .tolist()
            ),
            "sample_failure_reasons": (
                ";".join(
                    sorted(
                        set(
                            failed[
                                "sample_qc_fail_reasons"
                            ]
                            .astype(str)
                            .tolist()
                        )
                    )
                )
                if len(failed) > 0
                else "."
            ),
            "n_family_samples_with_warnings": (
                len(warned)
            ),
            "warning_samples": ",".join(
                warned["sample_id_qc"]
                .astype(str)
                .tolist()
            ),
            "sample_warning_reasons": (
                ";".join(
                    sorted(
                        set(
                            warned[
                                "sample_qc_warning_reasons"
                            ]
                            .astype(str)
                            .tolist()
                        )
                    )
                )
                if len(warned) > 0
                else "."
            ),
            "all_family_genotypes_pass": (
                len(failed) == 0
            ),
            "marker_requires_BAM_check": (
                len(failed) > 0
                or len(warned) > 0
            ),
        }
    )


marker_summary = pd.DataFrame(
    marker_summary_rows
)

marker_summary.to_csv(
    OUT_DIR
    / "LOW_MANUAL_REVIEW_marker_QC.summary.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)


# ============================================================
# 17. Event-level QC summary
# ============================================================

event_summary_rows = []

for event_id, group in marker_summary.groupby(
    "event_id",
    sort=False,
):
    first = group.iloc[0]

    event_summary_rows.append(
        {
            "event_id": event_id,
            "family_id": (
                first["family_id"]
            ),
            "event_type": (
                first["event_type"]
            ),
            "event_context": (
                first["event_context"]
            ),
            "effective_low_reason_codes": (
                first[
                    "effective_low_reason_codes"
                ]
            ),
            "n_context_markers": len(
                group
            ),
            "n_core_markers_in_context": int(
                (
                    group["context_role"]
                    == "core"
                ).sum()
            ),
            "n_markers_with_sample_failure": int(
                (
                    group[
                        "n_family_samples_failed"
                    ] > 0
                ).sum()
            ),
            "n_markers_with_warnings": int(
                (
                    group[
                        "n_family_samples_with_warnings"
                    ] > 0
                ).sum()
            ),
            "minimum_QUAL": (
                group["QUAL"].min()
            ),
            "minimum_QD": (
                group["QD"].min()
            ),
            "maximum_FS": (
                group["FS"].max()
            ),
            "maximum_SOR": (
                group["SOR"].max()
            ),
            "minimum_MQ": (
                group["MQ"].min()
            ),
            "minimum_MQRankSum": (
                group[
                    "MQRankSum"
                ].min()
            ),
            "minimum_ReadPosRankSum": (
                group[
                    "ReadPosRankSum"
                ].min()
            ),
            "minimum_queen_GQ": (
                group["queen_GQ"].min()
            ),
            "minimum_queen_min_AD": (
                group[
                    "queen_min_allele_depth"
                ].min()
            ),
            "minimum_offspring_GQ": (
                group[
                    "minimum_offspring_GQ"
                ].min()
            ),
            "minimum_offspring_major_AF": (
                group[
                    "minimum_offspring_major_AF"
                ].min()
            ),
            "all_context_markers_pass": (
                (
                    group[
                        "n_family_samples_failed"
                    ] == 0
                ).all()
            ),
            "snp_quality_review_result": (
                "SNP_QC_FAILURE_DETECTED"
                if (
                    group[
                        "n_family_samples_failed"
                    ] > 0
                ).any()
                else (
                    "PASS_WITH_BORDERLINE_MARKERS"
                    if (
                        group[
                            "n_family_samples_with_warnings"
                        ] > 0
                    ).any()
                    else
                    "ALL_CONTEXT_MARKERS_PASS"
                )
            ),
        }
    )


event_summary = pd.DataFrame(
    event_summary_rows
)

event_summary.to_csv(
    OUT_DIR
    / "LOW_MANUAL_REVIEW_event_SNP_QC.summary.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 18. Markers prioritized for BAM review
# ============================================================

bam_review_markers = marker_summary.loc[
    marker_summary[
        "marker_requires_BAM_check"
    ]
].copy()

bam_review_markers.to_csv(
    OUT_DIR
    / "LOW_MANUAL_REVIEW_markers_for_BAM_review.tsv.gz",
    sep="\t",
    index=False,
    compression="gzip",
)

bam_bed = pd.DataFrame(
    {
        "chromosome": (
            bam_review_markers[
                "chromosome"
            ]
        ),
        "start": (
            bam_review_markers[
                "position"
            ].astype(int) - 1
        ),
        "end": (
            bam_review_markers[
                "position"
            ].astype(int)
        ),
        "name": (
            bam_review_markers[
                "event_id"
            ].astype(str)
            + "|"
            + bam_review_markers[
                "context_role"
            ].astype(str)
        ),
    }
)

bam_bed.to_csv(
    OUT_DIR
    / "LOW_MANUAL_REVIEW_markers_for_BAM_review.bed",
    sep="\t",
    header=False,
    index=False,
)


# ============================================================
# 19. Print summary
# ============================================================

print("\nEvent SNP-QC results:")

print(
    event_summary[
        "snp_quality_review_result"
    ]
    .value_counts()
    .to_string()
)

print(
    "\nMarker-level sample failures:"
)

print(
    marker_summary[
        "n_family_samples_failed"
    ]
    .value_counts()
    .sort_index()
    .to_string()
)

print(
    f"\nOutputs saved under:\n"
    f"{OUT_DIR}"
)
