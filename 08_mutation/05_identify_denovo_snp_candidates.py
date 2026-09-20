from pathlib import Path
from collections import Counter
import gzip
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

VCF_ROOT = (
    PROJECT
    / "03_variants/gatk_direct_diploid_SB/families"
)

METADATA = (
    PROJECT
    / "00_rawdata/metadata/analysis_samples.noC3.tsv"
)

MAPPING_SUMMARY = (
    PROJECT
    / "02_mapping/stats/mapping_summary.tsv"
)

OUT_ROOT = (
    PROJECT
    / "06_mutations/snp_candidates"
)

SUMMARY_DIR = (
    PROJECT
    / "06_mutations/summary"
)

OUT_ROOT.mkdir(
    parents=True,
    exist_ok=True,
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 2. Inherited genotype-level thresholds
# ============================================================

# These values intentionally match the previous marker script.

# Minimum per-sample FORMAT/DP.
MIN_DP = 10

# Minimum per-sample genotype quality.
MIN_GQ = 30

# Reject abnormally high-depth calls, which are enriched in
# repeats, collapsed regions and CNVs.
MAX_DEPTH_MULTIPLIER = 2.5
DEFAULT_MAX_DP = 100

# A homozygous REF call must have ALT fraction <= 0.10.
HOM_REF_MAX_ALT_FRACTION = 0.05

# A homozygous ALT call must have ALT fraction >= 0.90.
HOM_ALT_MIN_ALT_FRACTION = 0.95

# This is the value actually used in the previous uploaded
# marker script. It requires at least two supporting reads on
# each relevant strand. To require "more than three reads per
# direction", change this value from 2 to 4.
MIN_STRAND_READS_PER_DIRECTION = 1

# Family inheritance pattern:
# queen is homozygous; all five offspring are callable;
# exactly four match the queen and one is discordant.
N_EXPECTED_OFFSPRING = 5
N_QUEEN_MATCHING_OFFSPRING = 4
N_DISCORDANT_OFFSPRING = 1


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
# 4. General helpers
# ============================================================

def capture_command(
    command,
):
    result = subprocess.run(
        [
            str(value)
            for value in command
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            "Command failed:\n"
            + " ".join(
                str(value)
                for value in command
            )
            + "\n\n"
            + result.stderr
        )

    return result.stdout


def run_command(
    command,
):
    print(
        "Running:",
        " ".join(
            str(value)
            for value in command
        ),
    )

    subprocess.run(
        [
            str(value)
            for value in command
        ],
        check=True,
    )


def remove_vcf_and_indexes(
    vcf_path: Path,
):
    for path in [
        vcf_path,
        Path(
            str(vcf_path) + ".tbi"
        ),
        Path(
            str(vcf_path) + ".csi"
        ),
    ]:
        if path.exists():
            path.unlink()


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


def parse_number(
    text,
):
    if text in {
        "",
        ".",
        None,
    }:
        return np.nan

    try:
        return float(
            text
        )
    except ValueError:
        return np.nan


def parse_sb(
    sb_text: str,
):
    """
    FORMAT/SB order:
    REF_forward, REF_reverse, ALT_forward, ALT_reverse.
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


def called_allele_has_bidirectional_support(
    sb: dict,
    allele_index: int,
) -> bool:
    if allele_index == 0:
        required_counts = [
            sb["ref_forward"],
            sb["ref_reverse"],
        ]
    elif allele_index == 1:
        required_counts = [
            sb["alt_forward"],
            sb["alt_reverse"],
        ]
    else:
        raise ValueError(
            f"Unexpected allele index: {allele_index}"
        )

    return (
        min(required_counts)
        >= MIN_STRAND_READS_PER_DIRECTION
    )


def parse_homozygous_call(
    sample_id: str,
    field: str,
    sample_max_dp: dict,
    require_diploid_gt: bool,
):
    """
    Return a high-quality homozygous/haploid allele call.

    Queen:
      require_diploid_gt=True, so only 0/0 or 1/1 is accepted.

    Offspring:
      accept either haploid GT (0 or 1) or diagnostic
      diploid-homozygous GT (0/0 or 1/1).

    Heterozygous GT (0/1) is always rejected.
    """

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
        or sb_text == "."
    ):
        return None

    gt = gt.replace(
        "|",
        "/",
    )

    if require_diploid_gt:
        if gt == "0/0":
            allele_index = 0
        elif gt == "1/1":
            allele_index = 1
        else:
            return None
    else:
        if gt in {
            "0",
            "0/0",
        }:
            allele_index = 0
        elif gt in {
            "1",
            "1/1",
        }:
            allele_index = 1
        else:
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

    # Retain the previous rule that effective allele depth
    # must independently satisfy the minimum depth threshold.
    if total_ad < MIN_DP:
        return None

    alt_fraction = (
        alt_depth
        / total_ad
    )

    if allele_index == 0:
        if (
            alt_fraction
            > HOM_REF_MAX_ALT_FRACTION
        ):
            return None

        called_depth = ref_depth
        called_fraction = (
            1.0
            - alt_fraction
        )
    else:
        if (
            alt_fraction
            < HOM_ALT_MIN_ALT_FRACTION
        ):
            return None

        called_depth = alt_depth
        called_fraction = alt_fraction

    sb = parse_sb(
        sb_text
    )

    if sb is None:
        return None

    if not called_allele_has_bidirectional_support(
        sb,
        allele_index,
    ):
        return None

    return {
        "gt": gt,
        "allele_index": allele_index,
        "dp": dp,
        "gq": gq,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "total_ad": total_ad,
        "alt_fraction": alt_fraction,
        "called_depth": called_depth,
        "called_fraction": called_fraction,
        "sb_raw": sb_text,
        "strand_bias_counts": sb,
    }


def substitution_class(
    ancestral_base: str,
    derived_base: str,
) -> str:
    pair = frozenset(
        [
            ancestral_base.upper(),
            derived_base.upper(),
        ]
    )

    if pair in {
        frozenset(
            [
                "A",
                "G",
            ]
        ),
        frozenset(
            [
                "C",
                "T",
            ]
        ),
    }:
        return "transition"

    return "transversion"


# ============================================================
# 5. Metadata and family structure
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

    if (
        len(offspring)
        != N_EXPECTED_OFFSPRING
    ):
        raise ValueError(
            f"{family}: expected "
            f"{N_EXPECTED_OFFSPRING} offspring, "
            f"found {len(offspring)}"
        )

    family_samples[family] = {
        "queen": queens[0],
        "offspring": offspring,
    }


# ============================================================
# 6. Sample-specific maximum depth
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
# 7. Scan family-specific PASS SNP VCFs
# ============================================================

candidate_rows = []
genotype_rows = []
filter_summary_rows = []

query_format = (
    "%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER"
    "\\t%INFO/QD\\t%INFO/FS\\t%INFO/SOR\\t%INFO/MQ"
    "\\t%INFO/MQRankSum\\t%INFO/ReadPosRankSum"
    "[\\t%GT:%DP:%GQ:%AD:%SB]"
    "\\n"
)

N_STATIC_FIELDS = 12

for family in families:
    vcf = get_family_vcf(
        family
    )

    if not vcf.exists():
        raise FileNotFoundError(
            f"{family} PASS SNP VCF not found: {vcf}"
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

    vcf_samples = [
        line.strip()
        for line in capture_command(
            [
                BCFTOOLS,
                "query",
                "-l",
                str(vcf),
            ]
        ).splitlines()
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

    family_counters = Counter()

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

    for line in process.stdout:
        fields = (
            line.rstrip("\n")
            .split("\t")
        )

        if (
            len(fields)
            != N_STATIC_FIELDS
            + len(vcf_samples)
        ):
            family_counters[
                "malformed_query_rows"
            ] += 1

            continue

        (
            chromosome,
            position_text,
            ref,
            alt,
            qual_text,
            filter_value,
            qd_text,
            fs_text,
            sor_text,
            mq_text,
            mq_rank_sum_text,
            read_pos_rank_sum_text,
        ) = fields[
            :N_STATIC_FIELDS
        ]

        position = int(
            position_text
        )

        genotype_fields = fields[
            N_STATIC_FIELDS:
        ]

        family_counters[
            "pass_biallelic_snp_sites_examined"
        ] += 1

        queen_call = parse_homozygous_call(
            queen,
            genotype_fields[
                sample_index[queen]
            ],
            sample_max_dp,
            require_diploid_gt=True,
        )

        if queen_call is None:
            family_counters[
                "queen_homozygous_filter_failed"
            ] += 1

            continue

        family_counters[
            "queen_homozygous_filter_passed"
        ] += 1

        offspring_calls = {}
        offspring_call_failed = False

        for child in offspring:
            child_call = parse_homozygous_call(
                child,
                genotype_fields[
                    sample_index[child]
                ],
                sample_max_dp,
                require_diploid_gt=False,
            )

            if child_call is None:
                offspring_call_failed = True
                break

            offspring_calls[
                child
            ] = child_call

        if offspring_call_failed:
            family_counters[
                "one_or_more_offspring_not_callable"
            ] += 1

            continue

        family_counters[
            "all_five_offspring_callable"
        ] += 1

        queen_allele_index = (
            queen_call[
                "allele_index"
            ]
        )

        matching_offspring = [
            child
            for child in offspring
            if (
                offspring_calls[
                    child
                ]["allele_index"]
                == queen_allele_index
            )
        ]

        discordant_offspring = [
            child
            for child in offspring
            if (
                offspring_calls[
                    child
                ]["allele_index"]
                != queen_allele_index
            )
        ]

        if not (
            len(matching_offspring)
            == N_QUEEN_MATCHING_OFFSPRING
            and len(discordant_offspring)
            == N_DISCORDANT_OFFSPRING
        ):
            family_counters[
                "offspring_4plus1_pattern_failed"
            ] += 1

            continue

        mutant_sample = (
            discordant_offspring[0]
        )

        mutant_call = offspring_calls[
            mutant_sample
        ]

        mutant_allele_index = (
            mutant_call[
                "allele_index"
            ]
        )

        queen_base = (
            ref
            if queen_allele_index == 0
            else alt
        )

        mutant_base = (
            ref
            if mutant_allele_index == 0
            else alt
        )

        candidate_id = (
            f"{family}_DNM_SNP_"
            f"{chromosome}_{position}_"
            f"{queen_base}_to_{mutant_base}_"
            f"{mutant_sample}"
        )

        candidate_rows.append(
            {
                "candidate_id": candidate_id,
                "family_id": family,
                "chromosome": chromosome,
                "position": position,
                "reference_allele": ref,
                "alternative_allele": alt,
                "queen_id": queen,
                "queen_allele_index": (
                    queen_allele_index
                ),
                "queen_allele": queen_base,
                "mutant_offspring_id": (
                    mutant_sample
                ),
                "mutant_allele_index": (
                    mutant_allele_index
                ),
                "mutant_allele": mutant_base,
                "mutation": (
                    f"{queen_base}>{mutant_base}"
                ),
                "substitution_class": (
                    substitution_class(
                        queen_base,
                        mutant_base,
                    )
                ),
                "n_queen_matching_offspring": (
                    len(matching_offspring)
                ),
                "matching_offspring_ids": (
                    ",".join(
                        matching_offspring
                    )
                ),
                "site_QUAL": parse_number(
                    qual_text
                ),
                "site_FILTER": filter_value,
                "site_QD": parse_number(
                    qd_text
                ),
                "site_FS": parse_number(
                    fs_text
                ),
                "site_SOR": parse_number(
                    sor_text
                ),
                "site_MQ": parse_number(
                    mq_text
                ),
                "site_MQRankSum": parse_number(
                    mq_rank_sum_text
                ),
                "site_ReadPosRankSum": (
                    parse_number(
                        read_pos_rank_sum_text
                    )
                ),
                "queen_GT": queen_call["gt"],
                "queen_DP": queen_call["dp"],
                "queen_GQ": queen_call["gq"],
                "queen_AD_REF": (
                    queen_call[
                        "ref_depth"
                    ]
                ),
                "queen_AD_ALT": (
                    queen_call[
                        "alt_depth"
                    ]
                ),
                "queen_called_allele_fraction": (
                    queen_call[
                        "called_fraction"
                    ]
                ),
                "queen_SB": (
                    queen_call[
                        "sb_raw"
                    ]
                ),
                "mutant_GT": mutant_call["gt"],
                "mutant_DP": mutant_call["dp"],
                "mutant_GQ": mutant_call["gq"],
                "mutant_AD_REF": (
                    mutant_call[
                        "ref_depth"
                    ]
                ),
                "mutant_AD_ALT": (
                    mutant_call[
                        "alt_depth"
                    ]
                ),
                "mutant_called_allele_depth": (
                    mutant_call[
                        "called_depth"
                    ]
                ),
                "mutant_called_allele_fraction": (
                    mutant_call[
                        "called_fraction"
                    ]
                ),
                "mutant_SB": (
                    mutant_call[
                        "sb_raw"
                    ]
                ),
            }
        )

        all_calls = {
            queen: queen_call,
            **offspring_calls,
        }

        for sample_id in expected_samples:
            call = all_calls[
                sample_id
            ]

            if sample_id == queen:
                sample_role = "queen"
                inheritance_class = (
                    "queen_homozygous"
                )
            elif sample_id == mutant_sample:
                sample_role = "offspring"
                inheritance_class = (
                    "discordant_offspring"
                )
            else:
                sample_role = "offspring"
                inheritance_class = (
                    "queen_matching_offspring"
                )

            genotype_rows.append(
                {
                    "candidate_id": candidate_id,
                    "family_id": family,
                    "chromosome": chromosome,
                    "position": position,
                    "sample_id": sample_id,
                    "sample_role": sample_role,
                    "inheritance_class": (
                        inheritance_class
                    ),
                    "GT": call["gt"],
                    "allele_index": (
                        call[
                            "allele_index"
                        ]
                    ),
                    "DP": call["dp"],
                    "max_DP_threshold": (
                        sample_max_dp[
                            sample_id
                        ]
                    ),
                    "GQ": call["gq"],
                    "AD_REF": (
                        call[
                            "ref_depth"
                        ]
                    ),
                    "AD_ALT": (
                        call[
                            "alt_depth"
                        ]
                    ),
                    "total_AD": (
                        call[
                            "total_ad"
                        ]
                    ),
                    "ALT_fraction": (
                        call[
                            "alt_fraction"
                        ]
                    ),
                    "called_allele_fraction": (
                        call[
                            "called_fraction"
                        ]
                    ),
                    "SB": call["sb_raw"],
                }
            )

        family_counters[
            "preliminary_denovo_snp_candidates"
        ] += 1

    return_code = process.wait()

    stderr = ""

    if process.stderr is not None:
        stderr = process.stderr.read()

    if return_code != 0:
        raise RuntimeError(
            f"bcftools query failed for "
            f"{family}:\n{stderr}"
        )

    for stage, n_sites in sorted(
        family_counters.items()
    ):
        filter_summary_rows.append(
            {
                "family_id": family,
                "stage": stage,
                "n_sites": n_sites,
            }
        )


# ============================================================
# 8. Candidate recurrence flags
# ============================================================

candidate_columns = [
    "candidate_id",
    "family_id",
    "chromosome",
    "position",
    "reference_allele",
    "alternative_allele",
    "queen_id",
    "queen_allele_index",
    "queen_allele",
    "mutant_offspring_id",
    "mutant_allele_index",
    "mutant_allele",
    "mutation",
    "substitution_class",
    "n_queen_matching_offspring",
    "matching_offspring_ids",
    "site_QUAL",
    "site_FILTER",
    "site_QD",
    "site_FS",
    "site_SOR",
    "site_MQ",
    "site_MQRankSum",
    "site_ReadPosRankSum",
    "queen_GT",
    "queen_DP",
    "queen_GQ",
    "queen_AD_REF",
    "queen_AD_ALT",
    "queen_called_allele_fraction",
    "queen_SB",
    "mutant_GT",
    "mutant_DP",
    "mutant_GQ",
    "mutant_AD_REF",
    "mutant_AD_ALT",
    "mutant_called_allele_depth",
    "mutant_called_allele_fraction",
    "mutant_SB",
]

candidate_df = pd.DataFrame(
    candidate_rows,
    columns=candidate_columns,
)

if not candidate_df.empty:
    site_group_columns = [
        "chromosome",
        "position",
        "reference_allele",
        "alternative_allele",
    ]

    candidate_df[
        "candidate_site_n_families"
    ] = (
        candidate_df.groupby(
            site_group_columns
        )["family_id"]
        .transform(
            "nunique"
        )
    )

    candidate_df[
        "candidate_site_n_records"
    ] = (
        candidate_df.groupby(
            site_group_columns
        )["candidate_id"]
        .transform(
            "size"
        )
    )

    candidate_df[
        "recurrent_candidate_site"
    ] = (
        candidate_df[
            "candidate_site_n_records"
        ]
        > 1
    )

    candidate_df[
        "review_status"
    ] = np.where(
        candidate_df[
            "recurrent_candidate_site"
        ],
        "MANUAL_REVIEW_RECURRENT_SITE",
        "PRELIMINARY_CANDIDATE",
    )

    candidate_df = candidate_df.sort_values(
        [
            "family_id",
            "chromosome",
            "position",
            "mutant_offspring_id",
        ]
    )
else:
    candidate_df[
        "candidate_site_n_families"
    ] = pd.Series(
        dtype=int
    )

    candidate_df[
        "candidate_site_n_records"
    ] = pd.Series(
        dtype=int
    )

    candidate_df[
        "recurrent_candidate_site"
    ] = pd.Series(
        dtype=bool
    )

    candidate_df[
        "review_status"
    ] = pd.Series(
        dtype=str
    )


# ============================================================
# 9. Write combined and family-level outputs
# ============================================================

combined_candidate_file = (
    OUT_ROOT
    / "all_families.preliminary_denovo_SNP_candidates.tsv.gz"
)

candidate_df.to_csv(
    combined_candidate_file,
    sep="\t",
    index=False,
    compression="gzip",
)

genotype_df = pd.DataFrame(
    genotype_rows
)

combined_genotype_file = (
    OUT_ROOT
    / "all_families.preliminary_denovo_SNP_genotypes.long.tsv.gz"
)

genotype_df.to_csv(
    combined_genotype_file,
    sep="\t",
    index=False,
    compression="gzip",
)

for family in families:
    family_dir = (
        OUT_ROOT
        / family
    )

    family_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    family_candidates = candidate_df.loc[
        candidate_df["family_id"]
        == family
    ].copy()

    family_tsv = (
        family_dir
        / (
            f"{family}."
            "preliminary_denovo_SNP_candidates.tsv.gz"
        )
    )

    family_bed = (
        family_dir
        / (
            f"{family}."
            "preliminary_denovo_SNP_candidates.bed"
        )
    )

    family_vcf = (
        family_dir
        / (
            f"{family}."
            "preliminary_denovo_SNP_candidates.vcf.gz"
        )
    )

    family_candidates.to_csv(
        family_tsv,
        sep="\t",
        index=False,
        compression="gzip",
    )

    with family_bed.open(
        "w"
    ) as bed_handle:
        for _, row in family_candidates.iterrows():
            position = int(
                row["position"]
            )

            bed_handle.write(
                f"{row['chromosome']}\t"
                f"{position - 1}\t"
                f"{position}\n"
            )

    remove_vcf_and_indexes(
        family_vcf
    )

    source_vcf = get_family_vcf(
        family
    )

    if family_candidates.empty:
        # Produce a valid header-only VCF instead of leaving a
        # stale candidate VCF from an earlier run.
        run_command(
            [
                BCFTOOLS,
                "view",
                "--header-only",
                "--output-type",
                "z",
                "--output-file",
                str(family_vcf),
                str(source_vcf),
            ]
        )
    else:
        run_command(
            [
                BCFTOOLS,
                "view",
                "--regions-file",
                str(family_bed),
                "--output-type",
                "z",
                "--output-file",
                str(family_vcf),
                str(source_vcf),
            ]
        )

    run_command(
        [
            BCFTOOLS,
            "index",
            "--force",
            "--tbi",
            str(family_vcf),
        ]
    )


# ============================================================
# 10. Summaries and recorded parameters
# ============================================================

filter_summary = pd.DataFrame(
    filter_summary_rows
)

filter_summary.to_csv(
    SUMMARY_DIR
    / "denovo_SNP_filter_summary.tsv",
    sep="\t",
    index=False,
)

family_summary_rows = []

for family in families:
    family_candidates = candidate_df.loc[
        candidate_df["family_id"]
        == family
    ]

    family_summary_rows.append(
        {
            "family_id": family,
            "queen_id": (
                family_samples[
                    family
                ]["queen"]
            ),
            "n_offspring": (
                len(
                    family_samples[
                        family
                    ]["offspring"]
                )
            ),
            "n_preliminary_denovo_SNP_candidates": (
                len(
                    family_candidates
                )
            ),
            "n_recurrent_candidate_records": (
                int(
                    family_candidates[
                        "recurrent_candidate_site"
                    ].sum()
                )
                if not family_candidates.empty
                else 0
            ),
        }
    )

family_summary = pd.DataFrame(
    family_summary_rows
)

family_summary.to_csv(
    SUMMARY_DIR
    / "denovo_SNP_candidate_summary.tsv",
    sep="\t",
    index=False,
)

parameter_rows = [
    {
        "parameter": "INPUT_VARIANTS",
        "value": (
            "family-specific PASS biallelic SNPs"
        ),
        "comment": (
            "Site-level GATK hard filters were "
            "applied upstream."
        ),
    },
    {
        "parameter": "MIN_DP",
        "value": MIN_DP,
        "comment": (
            "Inherited from the recombination "
            "marker workflow."
        ),
    },
    {
        "parameter": "MIN_GQ",
        "value": MIN_GQ,
        "comment": (
            "Inherited from the recombination "
            "marker workflow."
        ),
    },
    {
        "parameter": "MAX_DEPTH_MULTIPLIER",
        "value": MAX_DEPTH_MULTIPLIER,
        "comment": (
            "Sample-specific maximum DP equals "
            "ceil(mean depth times this value)."
        ),
    },
    {
        "parameter": "DEFAULT_MAX_DP",
        "value": DEFAULT_MAX_DP,
        "comment": (
            "Used only when sample mean depth is "
            "unavailable."
        ),
    },
    {
        "parameter": "HOM_REF_MAX_ALT_FRACTION",
        "value": HOM_REF_MAX_ALT_FRACTION,
        "comment": (
            "Equivalent to REF major-allele "
            "fraction >= 0.90."
        ),
    },
    {
        "parameter": "HOM_ALT_MIN_ALT_FRACTION",
        "value": HOM_ALT_MIN_ALT_FRACTION,
        "comment": (
            "Equivalent to ALT major-allele "
            "fraction >= 0.90."
        ),
    },
    {
        "parameter": "MIN_STRAND_READS_PER_DIRECTION",
        "value": (
            MIN_STRAND_READS_PER_DIRECTION
        ),
        "comment": (
            "Inherited actual value. Change to 4 "
            "to require more than three reads on "
            "each relevant strand."
        ),
    },
    {
        "parameter": "N_EXPECTED_OFFSPRING",
        "value": N_EXPECTED_OFFSPRING,
        "comment": (
            "All five offspring must be callable."
        ),
    },
    {
        "parameter": "N_QUEEN_MATCHING_OFFSPRING",
        "value": (
            N_QUEEN_MATCHING_OFFSPRING
        ),
        "comment": (
            "Exactly four offspring must match "
            "the homozygous queen."
        ),
    },
    {
        "parameter": "N_DISCORDANT_OFFSPRING",
        "value": N_DISCORDANT_OFFSPRING,
        "comment": (
            "Exactly one offspring must carry "
            "the opposite allele."
        ),
    },
]

pd.DataFrame(
    parameter_rows
).to_csv(
    SUMMARY_DIR
    / "denovo_SNP_calling_parameters.tsv",
    sep="\t",
    index=False,
)


# ============================================================
# 11. Final report
# ============================================================

print(
    "\nPreliminary de novo SNP candidate summary:"
)

print(
    family_summary.to_string(
        index=False
    )
)

print(
    "\nCombined candidate table:"
)

print(
    combined_candidate_file
)

print(
    "\nCombined per-sample genotype evidence:"
)

print(
    combined_genotype_file
)

print(
    "\nThese are preliminary candidates. "
    "Repeat/CNV/mappability masks and BAM-level "
    "review have not yet been applied."
)

