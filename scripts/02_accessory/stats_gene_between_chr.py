#!/usr/bin/env python3

import argparse
import pandas as pd
from itertools import combinations


def parse_gff3_attributes(attr_str):
    """Parse GFF3 attribute string into a dictionary."""
    attributes = {}
    for attr in attr_str.strip().split(";"):
        if "=" in attr:
            key, value = attr.split("=", 1)
            attributes[key] = value
    return attributes


def main():
    parser = argparse.ArgumentParser(description="Find shared genes between chromosomes")
    parser.add_argument("-c", "--cluster_file", required=True,
                        help="MMseqs easycluster output (TSV: cluster_rep, gene_id)")
    parser.add_argument("-g", "--gff3_file", required=True,
                        help="GFF3 annotation file")
    parser.add_argument("-s", "--selected_contig", default=None,
                        help="File with list of contigs to include (TSV, single column)")
    parser.add_argument("-o", "--output_file", default="gene_between_chr.tsv",
                        help="Output file for gene counts between chromosomes")

    args = parser.parse_args()

    print("Reading GFF3 file...")
    gene_to_chrom = {}
    with open(args.gff3_file, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom = parts[0]
            feature = parts[2]
            attributes = parse_gff3_attributes(parts[8])

            if feature == "mRNA" and "ID" in attributes:
                gene_id = attributes["ID"]
                gene_to_chrom[gene_id] = chrom

    print(f"Gene-to-chromosome mapping: {len(gene_to_chrom)} genes")

    print("Reading MMseqs cluster file...")
    cluster_to_genes = {}
    with open(args.cluster_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            cluster_rep = parts[0]
            gene_full_id = parts[1]

            gene_id = gene_full_id.split("|")[-1]

            if cluster_rep not in cluster_to_genes:
                cluster_to_genes[cluster_rep] = set()
            cluster_to_genes[cluster_rep].add(gene_id)

    print(f"Clusters: {len(cluster_to_genes)}")

    selected_contigs = set()
    if args.selected_contig:
        selected_contig_df = pd.read_csv(args.selected_contig, sep="\t", header=None)
        selected_contigs = set(selected_contig_df[0])
        print(f"Selected contigs: {len(selected_contigs)}")

    print("Finding genes shared between different chromosomes...")
    chr_pair_counts = {}

    for cluster_rep, genes in cluster_to_genes.items():
        chroms = set()
        for gene_id in genes:
            if gene_id in gene_to_chrom:
                chroms.add(gene_to_chrom[gene_id])

        if len(chroms) >= 2:
            chroms_sorted = sorted(chroms)
            for chrom1, chrom2 in combinations(chroms_sorted, 2):
                if not selected_contigs or (chrom1 in selected_contigs or chrom2 in selected_contigs):
                    pair_key = (chrom1, chrom2)
                    if pair_key not in chr_pair_counts:
                        chr_pair_counts[pair_key] = 0
                    chr_pair_counts[pair_key] += 1

    output_data = []
    for (chrom1, chrom2), count in sorted(chr_pair_counts.items()):
        output_data.append({"Contig1": chrom1, "Contig2": chrom2, "Gene_count": count})

    output_df = pd.DataFrame(output_data)
    output_df.to_csv(args.output_file, sep="\t", index=False)

    print(f"Gene count saved to: {args.output_file}")
    print(f"Total chromosome pairs with shared genes: {len(output_df)}")
    print("Done!")


if __name__ == "__main__":
    main()
