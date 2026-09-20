from pathlib import Path
import copy
import subprocess
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import (
        DistanceMatrix,
        DistanceTreeConstructor,
    )
except ImportError as error:
    raise ImportError(
        "Biopython is required. Install with:\n"
        "mamba install -c conda-forge biopython"
    ) from error


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

PCA_DIR = PROJECT / "04_pedigree/pca"
DIST_DIR = PROJECT / "04_pedigree/distance"
TREE_DIR = PROJECT / "04_pedigree/tree"
FIG_DIR = PROJECT / "04_pedigree/figures"
INTERMEDIATE_DIR = PROJECT / "04_pedigree/intermediate"

for directory in [
    PCA_DIR,
    DIST_DIR,
    TREE_DIR,
    FIG_DIR,
    INTERMEDIATE_DIR,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True
    )

MIN_DP = 10
DEFAULT_MAX_DP = 100

HOM_REF_MAX_AB = 0.10
HOM_ALT_MIN_AB = 0.90

# Queen heterozygous calls are accepted only when
# allele balance is reasonably centered.
QUEEN_HET_MIN_AB = 0.25
QUEEN_HET_MAX_AB = 0.75

MIN_CALL_RATE = 0.90
MIN_MAF = 0.05

# Physical thinning to reduce strong local LD
MIN_SNP_DISTANCE_BP = 20_000

MAX_RETAINED_SNPS = 50_000


# ============================================================
# 2. Metadata and sample-specific depth
# ============================================================

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str
)

metadata = metadata.sort_values(
    ["family_id", "role", "sample_id"]
).reset_index(drop=True)

role_lookup = (
    metadata
    .set_index("sample_id")["role"]
    .to_dict()
)

family_lookup = (
    metadata
    .set_index("sample_id")["family_id"]
    .to_dict()
)

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
# 3. VCF sample order
# ============================================================

vcf_samples = subprocess.check_output(
    [
        "bcftools",
        "query",
        "-l",
        str(VCF),
    ],
    text=True
).strip().splitlines()

if set(vcf_samples) != set(
    metadata["sample_id"]
):
    raise ValueError(
        "VCF and metadata sample names differ."
    )

sample_roles = [
    role_lookup[sample_id]
    for sample_id in vcf_samples
]


# ============================================================
# 4. Parse normalized dosage
# ============================================================

def parse_normalized_dosage(
    sample_id: str,
    role: str,
    field: str
):
    """
    Normalized alternate-allele dosage:

    Diploid queen:
        0/0 -> 0.0
        0/1 -> 0.5
        1/1 -> 1.0

    Haploid offspring represented in diploid-called VCF:
        0/0 -> 0.0
        1/1 -> 1.0
        0/1 -> missing
    """

    parts = field.split(":")

    if len(parts) < 3:
        return np.nan

    gt, dp_text, ad_text = parts[:3]

    if (
        gt in {".", "./.", ".|."}
        or dp_text == "."
        or ad_text == "."
    ):
        return np.nan

    try:
        dp = int(dp_text)
    except ValueError:
        return np.nan

    if dp < MIN_DP:
        return np.nan

    if dp > sample_max_dp.get(
        sample_id,
        DEFAULT_MAX_DP
    ):
        return np.nan

    ad_parts = ad_text.split(",")

    if len(ad_parts) != 2:
        return np.nan

    try:
        ref_depth = int(ad_parts[0])
        alt_depth = int(ad_parts[1])
    except ValueError:
        return np.nan

    total_ad = ref_depth + alt_depth

    if total_ad < MIN_DP:
        return np.nan

    allele_balance = alt_depth / total_ad
    normalized_gt = gt.replace("|", "/")

    if normalized_gt in {"0", "0/0"}:
        if allele_balance <= HOM_REF_MAX_AB:
            return 0.0

        return np.nan

    if normalized_gt in {"1", "1/1"}:
        if allele_balance >= HOM_ALT_MIN_AB:
            return 1.0

        return np.nan

    if normalized_gt in {"0/1", "1/0"}:
        if role == "queen":
            if (
                QUEEN_HET_MIN_AB
                <= allele_balance
                <= QUEEN_HET_MAX_AB
            ):
                return 0.5

        # Offspring heterozygous calls are treated
        # as errors or ambiguous regions.
        return np.nan

    return np.nan


# ============================================================
# 5. Stream and thin SNPs
# ============================================================

query_format = (
    "%CHROM\\t%POS"
    "[\\t%GT:%DP:%AD]"
    "\\n"
)

process = subprocess.Popen(
    [
        "bcftools",
        "query",
        "-f",
        query_format,
        str(VCF),
    ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    bufsize=1,
)

if process.stdout is None:
    raise RuntimeError(
        "Failed to open bcftools output."
    )

last_retained_position = {}

retained_vectors = []
retained_sites = []

n_input_sites = 0

for line in process.stdout:
    fields = line.rstrip("\n").split("\t")

    if len(fields) != 2 + len(vcf_samples):
        continue

    chromosome = fields[0]

    try:
        position = int(fields[1])
    except ValueError:
        continue

    n_input_sites += 1

    previous_position = last_retained_position.get(
        chromosome
    )

    if (
        previous_position is not None
        and position - previous_position
        < MIN_SNP_DISTANCE_BP
    ):
        continue

    genotype_fields = fields[2:]

    dosages = np.array(
        [
            parse_normalized_dosage(
                sample_id,
                role,
                genotype_field
            )
            for (
                sample_id,
                role,
                genotype_field
            ) in zip(
                vcf_samples,
                sample_roles,
                genotype_fields
            )
        ],
        dtype=float
    )

    called = np.isfinite(dosages)
    call_rate = called.mean()

    if call_rate < MIN_CALL_RATE:
        continue

    allele_frequency = np.nanmean(
        dosages
    )

    maf = min(
        allele_frequency,
        1.0 - allele_frequency
    )

    if maf < MIN_MAF:
        continue

    if np.nanstd(dosages) == 0:
        continue

    retained_vectors.append(
        dosages
    )

    retained_sites.append(
        {
            "chromosome": chromosome,
            "position": position,
            "call_rate": call_rate,
            "allele_frequency": allele_frequency,
            "minor_allele_frequency": maf,
        }
    )

    last_retained_position[
        chromosome
    ] = position

    if len(retained_vectors) >= MAX_RETAINED_SNPS:
        break

return_code = process.wait()

stderr = ""

if process.stderr is not None:
    stderr = process.stderr.read()

if return_code != 0:
    print(stderr, file=sys.stderr)

    raise RuntimeError(
        "bcftools query failed."
    )

if len(retained_vectors) < 100:
    raise RuntimeError(
        "Too few SNPs were retained for "
        "relationship analysis."
    )

print(f"Input SNPs examined: {n_input_sites:,}")
print(f"Retained SNPs: {len(retained_vectors):,}")


# ============================================================
# 6. Construct sample × SNP dosage matrix
# ============================================================

# Current shape:
# retained SNP × sample
site_by_sample = np.vstack(
    retained_vectors
)

# Desired shape:
# sample × retained SNP
dosage_matrix = site_by_sample.T

site_table = pd.DataFrame(
    retained_sites
)

site_table.to_csv(
    INTERMEDIATE_DIR
    / "relationship_snps.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 7. PCA with mean imputation and SNP standardization
# ============================================================

site_means = np.nanmean(
    dosage_matrix,
    axis=0
)

missing_rows, missing_columns = np.where(
    ~np.isfinite(dosage_matrix)
)

dosage_imputed = dosage_matrix.copy()

dosage_imputed[
    missing_rows,
    missing_columns
] = site_means[
    missing_columns
]

column_means = dosage_imputed.mean(
    axis=0
)

column_sd = dosage_imputed.std(
    axis=0,
    ddof=1
)

valid_columns = (
    np.isfinite(column_sd)
    & (column_sd > 0)
)

dosage_imputed = dosage_imputed[
    :,
    valid_columns
]

column_means = column_means[
    valid_columns
]

column_sd = column_sd[
    valid_columns
]

standardized = (
    dosage_imputed - column_means
) / column_sd

U, singular_values, Vt = np.linalg.svd(
    standardized,
    full_matrices=False
)

scores = U * singular_values

variance = singular_values ** 2

explained_variance_ratio = (
    variance / variance.sum()
)

n_components = min(
    10,
    scores.shape[1]
)

pca_columns = [
    f"PC{index + 1}"
    for index in range(n_components)
]

pca = pd.DataFrame(
    scores[:, :n_components],
    columns=pca_columns
)

pca.insert(
    0,
    "sample_id",
    vcf_samples
)

pca["family_id"] = pca[
    "sample_id"
].map(family_lookup)

pca["role"] = pca[
    "sample_id"
].map(role_lookup)

pca.to_csv(
    PCA_DIR
    / "mixed_ploidy_pca_scores.tsv",
    sep="\t",
    index=False
)

variance_table = pd.DataFrame(
    {
        "PC": [
            f"PC{index + 1}"
            for index in range(
                n_components
            )
        ],
        "explained_variance_ratio": (
            explained_variance_ratio[
                :n_components
            ]
        ),
    }
)

variance_table.to_csv(
    PCA_DIR
    / "mixed_ploidy_pca_variance.tsv",
    sep="\t",
    index=False
)


# ============================================================
# 8. PCA plot
# ============================================================

families = sorted(
    pca["family_id"].dropna().unique()
)

default_colors = (
    plt.rcParams[
        "axes.prop_cycle"
    ]
    .by_key()["color"]
)

family_colors = {
    family: default_colors[
        index % len(default_colors)
    ]
    for index, family in enumerate(
        families
    )
}

fig, ax = plt.subplots(
    figsize=(5.2, 4.2),
    dpi=150
)

for family in families:
    for role, marker in [
        ("queen", "s"),
        ("offspring", "o"),
    ]:
        subset = pca.loc[
            (pca["family_id"] == family)
            & (pca["role"] == role)
        ]

        if subset.empty:
            continue

        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            marker=marker,
            s=36 if role == "queen" else 24,
            label=f"{family} {role}",
            color=family_colors[family]
        )

        for _, row in subset.iterrows():
            ax.text(
                row["PC1"],
                row["PC2"],
                row["sample_id"],
                fontsize=5,
                ha="left",
                va="bottom"
            )

pc1_percent = (
    explained_variance_ratio[0] * 100
)

pc2_percent = (
    explained_variance_ratio[1] * 100
)

ax.set_xlabel(
    f"PC1 ({pc1_percent:.2f}%)"
)

ax.set_ylabel(
    f"PC2 ({pc2_percent:.2f}%)"
)

ax.legend(
    frameon=False,
    fontsize=6,
    ncol=2
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "mixed_ploidy_pca.pdf"
)

plt.close(fig)


# ============================================================
# 9. Pairwise mixed-ploidy genetic distances
# ============================================================

n_samples = len(vcf_samples)

distance_matrix = np.zeros(
    (n_samples, n_samples),
    dtype=float
)

for sample_i in range(n_samples):
    for sample_j in range(
        sample_i + 1,
        n_samples
    ):
        values_i = dosage_matrix[
            sample_i
        ]

        values_j = dosage_matrix[
            sample_j
        ]

        valid = (
            np.isfinite(values_i)
            & np.isfinite(values_j)
        )

        if valid.sum() == 0:
            distance = np.nan
        else:
            # Mean normalized allele-dosage difference
            # ranges approximately from 0 to 1.
            distance = np.mean(
                np.abs(
                    values_i[valid]
                    - values_j[valid]
                )
            )

        distance_matrix[
            sample_i,
            sample_j
        ] = distance

        distance_matrix[
            sample_j,
            sample_i
        ] = distance

distance_df = pd.DataFrame(
    distance_matrix,
    index=vcf_samples,
    columns=vcf_samples
)

distance_df.to_csv(
    DIST_DIR
    / "mixed_ploidy_genetic_distance.tsv",
    sep="\t"
)


# ============================================================
# 10. Genetic-distance heatmap
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.2, 6.4),
    dpi=150
)

image = ax.imshow(
    distance_matrix,
    aspect="equal"
)

ax.set_xticks(
    np.arange(n_samples)
)

ax.set_xticklabels(
    vcf_samples,
    rotation=90,
    fontsize=5
)

ax.set_yticks(
    np.arange(n_samples)
)

ax.set_yticklabels(
    vcf_samples,
    fontsize=5
)

colorbar = fig.colorbar(
    image,
    ax=ax,
    fraction=0.035,
    pad=0.03
)

colorbar.set_label(
    "Mean normalized genotype difference"
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "mixed_ploidy_genetic_distance_heatmap.pdf"
)

plt.close(fig)


# ============================================================
# 11. Neighbor-Joining tree
# ============================================================

lower_triangle = [
    distance_matrix[
        row_index,
        :row_index + 1
    ].tolist()
    for row_index in range(
        n_samples
    )
]

biopython_distance = DistanceMatrix(
    names=vcf_samples,
    matrix=lower_triangle
)

constructor = DistanceTreeConstructor()

nj_tree = constructor.nj(
    biopython_distance
)

Phylo.write(
    nj_tree,
    TREE_DIR
    / "mixed_ploidy_NJ_tree.nwk",
    "newick"
)

# Midpoint rooting is only used for visualization.
# It does not imply a biological outgroup.
plot_tree = copy.deepcopy(
    nj_tree
)

try:
    plot_tree.root_at_midpoint()
except Exception:
    pass

fig = plt.figure(
    figsize=(7.0, 7.5),
    dpi=150
)

ax = fig.add_subplot(111)

Phylo.draw(
    plot_tree,
    axes=ax,
    do_show=False
)

ax.set_xlabel(
    "Mean normalized genotype difference"
)

fig.tight_layout()

fig.savefig(
    FIG_DIR
    / "mixed_ploidy_NJ_tree.pdf"
)

plt.close(fig)

print(
    "Mixed-ploidy PCA, distance analysis "
    "and NJ tree completed."
)
