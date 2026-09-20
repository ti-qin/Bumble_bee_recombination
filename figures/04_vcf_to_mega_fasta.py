#!/usr/bin/env python3

import argparse
from pathlib import Path
import pysam


def gt_to_base(gt, ref, alt):
    if gt is None or any(x is None for x in gt):
        raise ValueError("Missing genotype detected")

    alleles = set(gt)

    if alleles == {0}:
        return ref
    if alleles == {1}:
        return alt

    raise ValueError(f"Non-homozygous genotype detected: {gt}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--thin-bp",
        type=int,
        default=0,
        help="Keep at most one SNP per physical interval; 0 keeps all SNPs",
    )
    args = parser.parse_args()

    vcf = pysam.VariantFile(str(args.vcf))
    samples = list(vcf.header.samples)
    sequences = {sample: [] for sample in samples}

    selected_records = []
    best_by_bin = {}

    for record in vcf.fetch():
        if record.alts is None or len(record.alts) != 1:
            continue

        ref = record.ref.upper()
        alt = record.alts[0].upper()

        if len(ref) != 1 or len(alt) != 1:
            continue

        if args.thin_bp <= 0:
            selected_records.append(record)
        else:
            key = (record.contig, (record.pos - 1) // args.thin_bp)
            old = best_by_bin.get(key)

            if old is None or (record.qual or 0) > (old.qual or 0):
                best_by_bin[key] = record

    if args.thin_bp > 0:
        selected_records = sorted(
            best_by_bin.values(),
            key=lambda x: (x.contig, x.pos),
        )

    for record in selected_records:
        ref = record.ref.upper()
        alt = record.alts[0].upper()

        for sample in samples:
            gt = record.samples[sample]["GT"]
            sequences[sample].append(gt_to_base(gt, ref, alt))

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w") as handle:
        for sample in samples:
            sequence = "".join(sequences[sample])
            handle.write(f">{sample}\n")

            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    print(f"Samples:       {len(samples)}")
    print(f"Selected SNPs: {len(selected_records):,}")
    print(f"Output FASTA:  {args.output}")


if __name__ == "__main__":
    main()
