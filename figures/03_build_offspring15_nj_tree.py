#!/usr/bin/env python3

import argparse
import copy
import csv
from pathlib import Path

import numpy as np
import pysam
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def gt_to_binary(gt):
    if gt is None or len(gt) == 0:
        return None
    if any(x is None for x in gt):
        return None
    if any(x not in (0, 1) for x in gt):
        return None
    if len(set(gt)) != 1:
        return None
    return int(gt[0])


def read_vcf(vcf_file, thin_bp):
    vcf = pysam.VariantFile(str(vcf_file))
    samples = list(vcf.header.samples)

    if len(samples) != 15:
        raise ValueError(f"Expected 15 samples, found {len(samples)}")

    sites = []
    best_by_bin = {}

    for record in vcf.fetch():
        if record.alts is None or len(record.alts) != 1:
            continue
        if len(record.ref) != 1 or len(record.alts[0]) != 1:
            continue

        genotypes = []
        valid = True

        for sample in samples:
            gt = gt_to_binary(record.samples[sample]["GT"])
            if gt is None:
                valid = False
                break
            genotypes.append(gt)

        if not valid or len(set(genotypes)) < 2:
            continue

        site = {
            "chrom": record.contig,
            "pos": record.pos,
            "ref": record.ref,
            "alt": record.alts[0],
            "qual": float(record.qual or 0),
            "gt": np.asarray(genotypes, dtype=np.int8),
        }

        if thin_bp <= 0:
            sites.append(site)
        else:
            key = (record.contig, (record.pos - 1) // thin_bp)
            old = best_by_bin.get(key)

            if old is None or site["qual"] > old["qual"]:
                best_by_bin[key] = site

    vcf.close()

    if thin_bp > 0:
        sites = list(best_by_bin.values())

    if not sites:
        raise ValueError("No usable SNPs remain.")

    return samples, sites


def calculate_distance(sites, n_samples):
    genotype_matrix = np.column_stack([site["gt"] for site in sites])
    distance = np.zeros((n_samples, n_samples), dtype=float)

    for i in range(n_samples):
        for j in range(i):
            value = np.mean(genotype_matrix[i] != genotype_matrix[j])
            distance[i, j] = value
            distance[j, i] = value

    return distance


def write_distance(path, samples, distance):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample"] + samples)

        for sample, row in zip(samples, distance):
            writer.writerow([sample] + [f"{x:.10f}" for x in row])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--thin-bp", type=int, default=5000)
    args = parser.parse_args()

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    samples, sites = read_vcf(args.vcf, args.thin_bp)
    distance = calculate_distance(sites, len(samples))

    distance_file = Path(f"{args.out_prefix}.distance_matrix.tsv")
    snp_file = Path(f"{args.out_prefix}.selected_snps.tsv")
    unrooted_file = Path(f"{args.out_prefix}.NJ.unrooted.nwk")
    midpoint_file = Path(f"{args.out_prefix}.NJ.midpoint_rooted.nwk")

    write_distance(distance_file, samples, distance)

    with open(snp_file, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["chromosome", "position", "ref", "alt", "qual"])

        for site in sites:
            writer.writerow([
                site["chrom"],
                site["pos"],
                site["ref"],
                site["alt"],
                site["qual"],
            ])

    lower_triangle = [
        [float(distance[i, j]) for j in range(i + 1)]
        for i in range(len(samples))
    ]

    dm = DistanceMatrix(names=samples, matrix=lower_triangle)
    tree = DistanceTreeConstructor().nj(dm)

    for clade in tree.get_nonterminals():
        clade.name = None

    tree.rooted = False
    Phylo.write(tree, unrooted_file, "newick")

    midpoint_tree = copy.deepcopy(tree)
    midpoint_tree.root_at_midpoint()

    for clade in midpoint_tree.get_nonterminals():
        clade.name = None

    Phylo.write(midpoint_tree, midpoint_file, "newick")

    print(f"Samples: {len(samples)}")
    print(f"Selected SNPs: {len(sites):,}")
    print(f"Distance matrix: {distance_file}")
    print(f"Unrooted tree:   {unrooted_file}")
    print(f"Midpoint tree:   {midpoint_file}")
    print()
    Phylo.draw_ascii(tree)


if __name__ == "__main__":
    main()
