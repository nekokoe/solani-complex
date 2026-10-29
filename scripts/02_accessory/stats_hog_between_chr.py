#!/usr/bin/env python3

import argparse
import pandas as pd
from itertools import combinations

def main():
    parser = argparse.ArgumentParser(description="Find shared HOGs between chromosomes")
    parser.add_argument("-g", "--gene_chrom_file", required=True,
                        help="Gene to chromosome affiliation file (TranscriptID, GeneID, ContigID)")
    parser.add_argument("-m", "--gene_hog_file", required=True,
                        help="Gene to HOG mapping file (GeneID, SampleID, HOG_ID, OG_ID)")
    parser.add_argument("-p", "--hog_interval_file", required=True,
                        help="HOG to presence rate interval file (HOG, Presence_Rate, Presence_Rate_Interval)")
    parser.add_argument("--selected_contig", default=None,
                        help="File with list of contigs to include (tsv, single column)")
    parser.add_argument("--sample_info", default=None,
                        help="Sample info file (tsv: BMU_ID, MedEnv, Host, Species)")
    parser.add_argument("-o", "--output_prefix", default="hog_between_chr",
                        help="Output file prefix")
    parser.add_argument("--min_presence_rate", type=float, default=80,
                        help="Filter out HOGs with Presence_Rate >= this value (default: 80)")

    args = parser.parse_args()

    print("Reading input files...")

    first_line = open(args.gene_chrom_file).readline().strip()
    has_header = first_line.startswith("TranscriptID")
    if has_header:
        gene_chrom = pd.read_csv(args.gene_chrom_file, sep="\t")
    else:
        gene_chrom = pd.read_csv(args.gene_chrom_file, sep="\t", header=None,
                                  names=["TranscriptID", "GeneID", "ContigID"])

    gene_chrom.columns = ["TranscriptID", "GeneID", "ContigID"]

    first_line = open(args.gene_hog_file).readline().strip()
    has_header = first_line.startswith("GeneID")
    if has_header:
        gene_hog = pd.read_csv(args.gene_hog_file, sep="\t")
    else:
        gene_hog = pd.read_csv(args.gene_hog_file, sep="\t", header=None,
                                names=["GeneID", "SampleID", "HOG_ID", "OG_ID"])

    first_line = open(args.hog_interval_file).readline().strip()
    has_header = first_line.startswith("HOG")
    if has_header:
        hog_interval = pd.read_csv(args.hog_interval_file, sep="\t")
    else:
        hog_interval = pd.read_csv(args.hog_interval_file, sep="\t", header=None,
                                     names=["HOG", "Presence_Rate", "Presence_Rate_Interval"])

    print(f"Gene-Chrom mapping: {len(gene_chrom)} rows")
    print(f"Gene-HOG mapping: {len(gene_hog)} rows")
    print(f"HOG-Interval mapping: {len(hog_interval)} rows")

    print(f"Filtering HOGs with Presence_Rate >= {args.min_presence_rate}...")
    hog_interval_filtered = hog_interval[hog_interval["Presence_Rate"] < args.min_presence_rate]
    print(f"HOGs after filtering: {len(hog_interval_filtered)}")

    hog_set = set(hog_interval_filtered["HOG"])

    print("Merging data...")
    merged = gene_chrom[["TranscriptID", "ContigID"]].merge(
        gene_hog[["GeneID", "HOG_ID"]], left_on="TranscriptID", right_on="GeneID"
    )
    merged = merged[merged["HOG_ID"].isin(hog_set)]
    print(f"Merged data: {len(merged)} rows")

    selected_contigs = set()
    if args.selected_contig:
        selected_contig_df = pd.read_csv(args.selected_contig, sep="\t", header=None)
        selected_contigs = set(selected_contig_df[0])
        print(f"Selected contigs: {len(selected_contigs)}")

    print("Finding HOGs shared between different contigs...")
    hog_contig = merged.groupby("HOG_ID")["ContigID"].apply(set).reset_index()
    hog_contig = hog_contig[hog_contig["ContigID"].apply(len) >= 2]

    print(f"HOGs appearing on >= 2 contigs: {len(hog_contig)}")

    all_contigs = set()
    for contigs_set in hog_contig["ContigID"]:
        all_contigs.update(contigs_set)
    print(f"Total unique contigs: {len(all_contigs)}")

    contig_annotation = pd.DataFrame({"ContigID": list(all_contigs)})
    contig_annotation["Core_Accessory"] = contig_annotation["ContigID"].apply(
        lambda x: "Core" if x not in selected_contigs else "Accessory"
    )

    hog_pairs = []
    for _, row in hog_contig.iterrows():
        hog_id = row["HOG_ID"]
        contigs = sorted(row["ContigID"])
        for contig1, contig2 in combinations(contigs, 2):
            if not selected_contigs or (contig1 in selected_contigs and contig2 in selected_contigs):
                hog_pairs.append({"Contig1": contig1, "Contig2": contig2, "HOG": hog_id})

    hog_pairs_df = pd.DataFrame(hog_pairs)
    hog_pairs_df = hog_pairs_df.sort_values(["Contig1", "Contig2", "HOG"])

    if args.sample_info:
        print("Reading sample info file...")
        sample_info = pd.read_csv(args.sample_info, sep="\t")

        gene_to_contig = gene_chrom[["TranscriptID", "ContigID"]].drop_duplicates()
        gene_to_sample = gene_hog[["GeneID", "SampleID"]].drop_duplicates()

        contig_to_sample = gene_to_contig.merge(
            gene_to_sample, left_on="TranscriptID", right_on="GeneID", how="inner"
        )

        contig_to_sample = contig_to_sample[["ContigID", "SampleID"]].drop_duplicates()

        contig_to_info = contig_to_sample.merge(
            sample_info, left_on="SampleID", right_on="BMU_ID", how="left"
        )

        contig_to_info = contig_to_info.dropna(subset=["ContigID"])
        contig_to_info = contig_to_info.drop_duplicates(subset=["ContigID"])

        contig_annotation = contig_annotation.merge(
            contig_to_info[["ContigID", "SampleID", "Species", "Host", "MedEnv"]],
            on="ContigID",
            how="left"
        )

        contig_annotation = contig_annotation[["ContigID", "SampleID", "Species", "Host", "MedEnv", "Core_Accessory"]]
    else:
        contig_annotation = contig_annotation[["ContigID", "Core_Accessory"]]

    output_annotation = f"{args.output_prefix}_annotation.tsv"
    contig_annotation.to_csv(output_annotation, sep="\t", index=False)
    print(f"Contig annotation saved to: {output_annotation}")

    output_list = f"{args.output_prefix}_list.tsv"
    hog_pairs_df.to_csv(output_list, sep="\t", index=False)
    print(f"HOG list saved to: {output_list}")

    hog_count = hog_pairs_df.groupby(["Contig1", "Contig2"]).size().reset_index(name="HOG_count")
    hog_count = hog_count.sort_values(["Contig1", "Contig2"])

    output_count = f"{args.output_prefix}_count.tsv"
    hog_count.to_csv(output_count, sep="\t", index=False)
    print(f"HOG count saved to: {output_count}")

    print("Done!")

if __name__ == "__main__":
    main()
