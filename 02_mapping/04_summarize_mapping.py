from pathlib import Path
import re
import pandas as pd

project = Path("/localdata/qinti/Project/Bee")

manifest = project / "00_rawdata/metadata/clean_samples.tsv"
stats_dir = project / "02_mapping/stats"
depth_dir = project / "02_mapping/depth"
out_file = project / "02_mapping/stats/mapping_summary.tsv"

samples = pd.read_csv(manifest, sep="\t")

rows = []


def parse_flagstat(flagstat_file):
    total_reads = None
    mapped_reads = None
    mapped_rate = None
    properly_paired_rate = None

    if not flagstat_file.exists():
        return total_reads, mapped_reads, mapped_rate, properly_paired_rate

    text = flagstat_file.read_text()

    m_total = re.search(
        r"^(\d+) \+ \d+ in total",
        text,
        re.M
    )

    m_mapped = re.search(
        r"^(\d+) \+ \d+ mapped \(([\d.]+)%",
        text,
        re.M
    )

    m_pair = re.search(
        r"^(\d+) \+ \d+ properly paired \(([\d.]+)%",
        text,
        re.M
    )

    if m_total:
        total_reads = int(m_total.group(1))

    if m_mapped:
        mapped_reads = int(m_mapped.group(1))
        mapped_rate = float(m_mapped.group(2))

    if m_pair:
        properly_paired_rate = float(m_pair.group(2))

    return total_reads, mapped_reads, mapped_rate, properly_paired_rate


def parse_samtools_stats(stats_file):
    """
    Parse samtools stats SN lines.

    Example:
    SN    raw total sequences:    12345678    # excluding supplementary and secondary reads

    The numeric value is the third tab-separated field, not the last field.
    """
    raw_total = None
    reads_duplicated = None

    if not stats_file.exists():
        return None

    with stats_file.open() as handle:
        for line in handle:
            if not line.startswith("SN\t"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            label = parts[1].rstrip(":")
            value = parts[2]

            if label == "raw total sequences":
                raw_total = int(value)

            elif label == "reads duplicated":
                reads_duplicated = int(value)

    if raw_total is not None and reads_duplicated is not None and raw_total > 0:
        duplicate_rate = reads_duplicated / raw_total * 100
    else:
        duplicate_rate = None

    return duplicate_rate


def parse_mosdepth_summary(mosdepth_file):
    mean_depth = None

    if not mosdepth_file.exists():
        return mean_depth

    md = pd.read_csv(
        mosdepth_file,
        sep="\t"
    )

    if "chrom" not in md.columns or "mean" not in md.columns:
        return mean_depth

    total_row = md.loc[md["chrom"] == "total"]

    if len(total_row) == 1:
        mean_depth = float(total_row["mean"].iloc[0])

    return mean_depth


for _, row in samples.iterrows():
    sample_id = row["sample_id"]

    flagstat_file = stats_dir / f"{sample_id}.flagstat.txt"
    stats_file = stats_dir / f"{sample_id}.samtools.stats.txt"
    mosdepth_file = depth_dir / f"{sample_id}.mosdepth.summary.txt"

    (
        total_reads,
        mapped_reads,
        mapped_rate,
        properly_paired_rate
    ) = parse_flagstat(flagstat_file)

    duplicate_rate = parse_samtools_stats(stats_file)

    mean_depth = parse_mosdepth_summary(mosdepth_file)

    rows.append(
        {
            "sample_id": sample_id,
            "family_id": row["family_id"],
            "role": row["role"],
            "total_reads": total_reads,
            "mapped_reads": mapped_reads,
            "mapped_rate_percent": mapped_rate,
            "properly_paired_percent": properly_paired_rate,
            "duplicate_rate_percent": duplicate_rate,
            "mean_depth": mean_depth,
        }
    )

summary = pd.DataFrame(rows)

summary.to_csv(
    out_file,
    sep="\t",
    index=False
)

print(summary.to_string(index=False))
print(f"\nSaved to: {out_file}")
