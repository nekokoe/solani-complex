#!/usr/bin/env python3

import argparse
import os
import warnings
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Get HOG list on specified chromosomes")
    parser.add_argument("-g", "--gene_hog_file", required=True,
                        help="Gene-HOG mapping file (TSV, header: GeneID, SampleID, HOG_ID, OG_ID)")
    parser.add_argument("-m", "--mapping_file", required=True,
                        help="Gene-chromosome mapping file (TSV, header: TranscriptID, GeneID, ContigID)")
    parser.add_argument("-c", "--contig_list", required=True,
                        help="Chromosome ID list file (TSV, single column, no header)")
    parser.add_argument("-o", "--output_file", default="chr_hog_list.tsv",
                        help="Output HOG list file (TSV, single column, no header)")

    args = parser.parse_args()

    # Read gene-HOG mapping
    print("Reading gene-HOG mapping...")
    gene_hog_df = pd.read_csv(args.gene_hog_file, sep="\t", header=0)
    print(f"Total gene-HOG records: {len(gene_hog_df)}")

    # Read gene-chromosome mapping
    print("Reading gene-chromosome mapping...")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t", header=0)
    # Build gene -> contig mapping (using TranscriptID from file2 = GeneID from file1)
    gene_to_contig = dict(zip(mapping_df["TranscriptID"], mapping_df["ContigID"]))
    print(f"Gene-to-chromosome mapping: {len(gene_to_contig)} genes")

    # Read selected contigs
    print("Reading contig list...")
    if os.path.exists(args.contig_list):
        contig_df = pd.read_csv(args.contig_list, sep="\t", header=None)
        selected_contigs = set(contig_df[0])
    else:
        warnings.warn(
            f"File '{args.contig_list}' not found. Treating it as a single contig ID.",
            UserWarning
        )
        selected_contigs = {args.contig_list}
    print(f"Selected contigs: {len(selected_contigs)}")

    # Find genes on selected chromosomes
    genes_on_selected_chr = set()
    for gene_id, contig_id in gene_to_contig.items():
        if contig_id in selected_contigs:
            genes_on_selected_chr.add(gene_id)

    print(f"Genes on selected chromosomes: {len(genes_on_selected_chr)}")

    # Get HOGs for these genes
    filtered_hog_df = gene_hog_df[gene_hog_df["GeneID"].isin(genes_on_selected_chr)]
    hog_list = filtered_hog_df["HOG_ID"].unique()

    # Write output
    print(f"Unique HOGs on selected chromosomes: {len(hog_list)}")
    with open(args.output_file, "w") as f:
        for hog_id in sorted(hog_list):
            f.write(f"{hog_id}\n")

    print(f"HOG list saved to: {args.output_file}")


if __name__ == "__main__":
    main()
