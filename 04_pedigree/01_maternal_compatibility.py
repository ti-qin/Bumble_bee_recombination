from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 1. Paths and parameters
# ============================================================

PROJECT = Path("/localdata/qinti/Project/Bee")

VCF = (
    PROJECT
    / "03_variants/filtered/"
      "all_samples.diploid_for_ploidy.biallelic.snp.qual30.vcf.gz"
)

METADATA = PROJECT / "00_rawdata/metadata/clean_samples.tsv"

MAPPING_SUMMARY = (
    PROJECT
    / "02_mapping/stats/mapping_summary.tsv"
)

OUT_DIR = PROJECT / "04_pedigree/compatibility"
FIG_DIR = PROJECT / "04_pedigree/figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

MIN_DP = 10
DEFAULT_MAX_DP = 100

# Strict allele-balance criteria for homozygous calls
HOM_REF_MAX_AB = 0.10
HOM_ALT_MIN_AB = 0.90


# ============================================================
# 2. Read metadata
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str
)

required_columns = {
    "sample_id",
    "family_id",
    "role",
}

missing_columns = required_columns - set(metadata.columns)

if missing_columns:
    raise ValueError(
        f"Missing metadata columns: {sorted(missing_columns)}"
    )

queens = (
    metadata.loc[
        metadata["role"] == "queen",
        "sample_id"
    ]
    .sort_values()
    .tolist()
)

offspring = (
    metadata.loc[
        metadata["role"] == "offspring",
        ["sample_id", "family_id"]
    ]
    .sort_values(["family_id", "sample_id"])
    .reset_index(drop=True)
)

offspring_ids = offspring["sample_id"].tolist()

if len(queens) == 0 or len(offspring_ids) == 0:
    raise ValueError(
        "No queen or offspring samples were found."
    )

print(f"Queens: {len(queens)}")
print(f"Offspring: {len(offspring_ids)}")
print("Queen samples:", ", ".join(queens))


# ============================================================
# 3. Sample-specific depth limits
# ============================================================

sample_max_dp = {
    sample_id: DEFAULT_MAX_DP
    for sample_id in metadata["sample_id"]
}

if MAPPING_SUMMARY.exists():
    mapping = pd.read_csv(
        MAPPING_SUMMARY,
        sep="\t"
    )

    if {
        "sample_id",
        "mean_depth"
    }.issubset(mapping.columns):

        for _, row in mapping.iterrows():
            sample_id = str(row["sample_id"])

            try:
                mean_depth = float(row["mean_depth"])
            except (TypeError, ValueError):
                continue

            if np.isfinite(mean_depth) and mean_depth > 0:
                sample_max_dp[sample_id] = max(
                    30,
                    int(np.ceil(mean_depth * 2.5))
                )


# ============================================================
# 4. Check VCF sample order
# ============================================================

sample_cmd = [
    "bcftools",
    "query",
    "-l",
    str(VCF),
]

vcf_samples = subprocess.check_output(
    sample_cmd,
    text=True
).strip().splitlines()

metadata_samples = set(metadata["sample_id"])

if set(vcf_samples) != metadata_samples:
    missing_in_vcf = sorted(
        metadata_samples - set(vcf_samples)
    )

    extra_in_vcf = sorted(
        set(vcf_samples) - metadata_samples
    )

    raise ValueError(
        "VCF and metadata samples differ.\n"
        f"Missing in VCF: {missing_in_vcf}\n"
        f"Extra in VCF: {extra_in_vcf}"
    )

sample_index = {
    sample: index
    for index, sample in enumerate(vcf_samples)
}


# ============================================================
# 5. Parse a high-confidence homozygous allele
# ============================================================

def parse_homozygous_allele(
    sample_id: str,
    field: str
):
    """
    Return:
        0 for high-confidence reference homozygote
        1 for high-confidence alternative homozygote
        None for missing, heterozygous or ambiguous calls
    """

    parts = field.split(":")

    if len(parts) < 3:
        return None

    gt, dp_text, ad_text = parts[:3]

    if (
        gt in {".", "./.", ".|."}
        or dp_text == "."
        or ad_text == "."
    ):
        return None

    try:
        dp = int(dp_text)
    except ValueError:
        return None

    if dp < MIN_DP:
        return None

    if dp > sample_max_dp.get(
        sample_id,
        DEFAULT_MAX_DP
    ):
        return None

    ad_parts = ad_text.split(",")

    if len(ad_parts) != 2:
        return None

    try:
        ref_depth = int(ad_parts[0])
        alt_depth = int(ad_parts[1])
    except ValueError:
        return None

    total_ad = ref_depth + alt_depth

    if total_ad < MIN_DP:
        return None

    allele_balance = alt_depth / total_ad

    normalized_gt = gt.replace("|", "/")

    if normalized_gt in {"0", "0/0"}:
        if allele_balance <= HOM_REF_MAX_AB:
            return 0

    if normalized_gt in {"1", "1/1"}:
        if allele_balance >= HOM_ALT_MIN_AB:
            return 1

    return None


# ============================================================
# 6. Initialize counters
# ============================================================

counts = {}

for offspring_id in offspring_ids:
    for queen_id in queens:
        counts[(offspring_id, queen_id)] = {
            "informative_sites": 0,
            "mismatches": 0,
            "queen_hom_ref_sites": 0,
            "queen_hom_alt_sites": 0,
        }


# ============================================================
# 7. Stream genotypes from VCF
# ============================================================

query_format = (
    "%CHROM\\t%POS"
    "[\\t%GT:%DP:%AD]"
    "\\n"
)

query_cmd = [
    "bcftools",
    "query",
    "-f",
    query_format,
    str(VCF),
]

process = subprocess.Popen(
    query_cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    bufsize=1,
)

if process.stdout is None:
    raise RuntimeError(
        "Failed to open bcftools query output."
    )

n_sites = 0

for line in process.stdout:
    fields = line.rstrip("\n").split("\t")

    if len(fields) != 2 + len(vcf_samples):
        continue

    n_sites += 1

    genotype_fields = fields[2:]

    homozygous_calls = {}

    for sample_id in queens + offspring_ids:
        index = sample_index[sample_id]

        homozygous_calls[sample_id] = (
            parse_homozygous_allele(
                sample_id,
                genotype_fields[index]
            )
        )

    for queen_id in queens:
        queen_allele = homozygous_calls[queen_id]

        # A heterozygous or ambiguous queen site is
        # not informative for candidate-mother assignment.
        if queen_allele is None:
            continue

        for offspring_id in offspring_ids:
            offspring_allele = homozygous_calls[
                offspring_id
            ]

            if offspring_allele is None:
                continue

            counter = counts[
                (offspring_id, queen_id)
            ]

            counter["informative_sites"] += 1

            if queen_allele == 0:
                counter["queen_hom_ref_sites"] += 1
            else:
                counter["queen_hom_alt_sites"] += 1

            if offspring_allele != queen_allele:
                counter["mismatches"] += 1

return_code = process.wait()

stderr = ""

if process.stderr is not None:
    stderr = process.stderr.read()

if return_code != 0:
    print(stderr, file=sys.stderr)

    raise RuntimeError(
        "bcftools query failed."
    )

print(f"Processed variant sites: {n_sites:,}")


# ============================================================
# 8. Construct long-format result
# ============================================================

rows = []

family_lookup = (
    metadata
    .set_index("sample_id")["family_id"]
    .to_dict()
)

for offspring_id in offspring_ids:
    expected_queen = (
        f"{family_lookup[offspring_id]}_Q"
    )

    for queen_id in queens:
        counter = counts[
            (offspring_id, queen_id)
        ]

        informative = counter[
            "informative_sites"
        ]

        mismatches = counter[
            "mismatches"
        ]

        if informative > 0:
            mismatch_rate = (
                mismatches / informative
            )

            compatibility = (
                1.0 - mismatch_rate
            )
        else:
            mismatch_rate = np.nan
            compatibility = np.nan

        rows.append(
            {
                "offspring_id": offspring_id,
                "family_id": family_lookup[
                    offspring_id
                ],
                "expected_queen": expected_queen,
                "candidate_queen": queen_id,
                "informative_sites": informative,
                "queen_hom_ref_sites": counter[
                    "queen_hom_ref_sites"
                ],
                "queen_hom_alt_sites": counter[
                    "queen_hom_alt_sites"
                ],
                "mismatches": mismatches,
                "mismatch_rate": mismatch_rate,
                "compatibility": compatibility,
            }
        )

result = pd.DataFrame(rows)

result.to_csv(
    OUT_DIR / "maternal_compatibility_long.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 9. Compatibility and mismatch matrices
# ============================================================

compatibility_matrix = (
    result.pivot(
        index="offspring_id",
        columns="candidate_queen",
        values="compatibility"
    )
    .reindex(
        index=offspring_ids,
        columns=queens
    )
)

mismatch_matrix = (
    result.pivot(
        index="offspring_id",
        columns="candidate_queen",
        values="mismatch_rate"
    )
    .reindex(
        index=offspring_ids,
        columns=queens
    )
)

informative_matrix = (
    result.pivot(
        index="offspring_id",
        columns="candidate_queen",
        values="informative_sites"
    )
    .reindex(
        index=offspring_ids,
        columns=queens
    )
)

compatibility_matrix.to_csv(
    OUT_DIR / "maternal_compatibility_matrix.tsv",
    sep="\t"
)

mismatch_matrix.to_csv(
    OUT_DIR / "maternal_mismatch_rate_matrix.tsv",
    sep="\t"
)

informative_matrix.to_csv(
    OUT_DIR / "maternal_informative_sites_matrix.tsv",
    sep="\t"
)


# ============================================================
# 10. Assign the most compatible queen
# ============================================================

assignment_rows = []

for offspring_id in offspring_ids:
    values = (
        compatibility_matrix
        .loc[offspring_id]
        .dropna()
        .sort_values(
            ascending=False
        )
    )

    expected_queen = (
        f"{family_lookup[offspring_id]}_Q"
    )

    if len(values) == 0:
        best_queen = None
        best_value = np.nan
        second_value = np.nan
        margin = np.nan
    else:
        best_queen = values.index[0]
        best_value = values.iloc[0]

        if len(values) > 1:
            second_value = values.iloc[1]
            margin = (
                best_value - second_value
            )
        else:
            second_value = np.nan
            margin = np.nan

    assignment_rows.append(
        {
            "offspring_id": offspring_id,
            "family_id": family_lookup[
                offspring_id
            ],
            "expected_queen": expected_queen,
            "assigned_queen": best_queen,
            "assignment_matches_expected": (
                best_queen == expected_queen
            ),
            "best_compatibility": best_value,
            "second_best_compatibility": second_value,
            "compatibility_margin": margin,
        }
    )

assignment = pd.DataFrame(
    assignment_rows
)

assignment.to_csv(
    OUT_DIR / "maternal_assignment_summary.tsv",
    sep="\t",
    index=False
)

print("\nMaternal assignment summary:")
print(assignment.to_string(index=False))


# ============================================================
# 11. Plot compatibility heatmap
# ============================================================

matrix_values = compatibility_matrix.to_numpy(
    dtype=float
)

fig, ax = plt.subplots(
    figsize=(4.2, 6.0),
    dpi=150
)

image = ax.imshow(
    matrix_values,
    aspect="auto"
)

ax.set_xticks(
    np.arange(len(queens))
)

ax.set_xticklabels(
    queens,
    rotation=45,
    ha="right"
)

ax.set_yticks(
    np.arange(len(offspring_ids))
)

ax.set_yticklabels(
    offspring_ids,
    fontsize=6
)

ax.set_xlabel("Candidate queen")
ax.set_ylabel("Offspring")

for row_index in range(
    len(offspring_ids)
):
    for column_index in range(
        len(queens)
    ):
        value = matrix_values[
            row_index,
            column_index
        ]

        if np.isfinite(value):
            ax.text(
                column_index,
                row_index,
                f"{value:.4f}",
                ha="center",
                va="center",
                fontsize=4.5
            )

colorbar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.035,
    pad=0.03
)

colorbar.set_label(
    "Maternal genotype compatibility"
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "maternal_compatibility_heatmap.pdf"
)

plt.close(fig)


# ============================================================
# 12. Plot mismatch-rate heatmap
# ============================================================

mismatch_values = mismatch_matrix.to_numpy(
    dtype=float
)

fig, ax = plt.subplots(
    figsize=(4.2, 6.0),
    dpi=150
)

image = ax.imshow(
    mismatch_values,
    aspect="auto"
)

ax.set_xticks(
    np.arange(len(queens))
)

ax.set_xticklabels(
    queens,
    rotation=45,
    ha="right"
)

ax.set_yticks(
    np.arange(len(offspring_ids))
)

ax.set_yticklabels(
    offspring_ids,
    fontsize=6
)

ax.set_xlabel("Candidate queen")
ax.set_ylabel("Offspring")

for row_index in range(
    len(offspring_ids)
):
    for column_index in range(
        len(queens)
    ):
        value = mismatch_values[
            row_index,
            column_index
        ]

        if np.isfinite(value):
            ax.text(
                column_index,
                row_index,
                f"{value:.4e}",
                ha="center",
                va="center",
                fontsize=4.2
            )

colorbar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.035,
    pad=0.03
)

colorbar.set_label(
    "Opposite-homozygote mismatch rate"
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "maternal_mismatch_rate_heatmap.pdf"
)

plt.close(fig)

print(
    "\nCompatibility analysis completed."
)
