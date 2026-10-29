#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a master effector table for all FSSC isolates.

For every isolate in the 211-sample list, reads the per-isolate EffectorP output
(per_isolates/<sample>/signalp/effectorp.tsv) and produces:

1. effector_master.gene_level.tsv  - one row per gene:
   GeneID, SampleID, Cytoplasmic_score, Apoplastic_score, NonEffector_score, Prediction, EffectorClass
2. effector_master.hog_level.tsv   - one row per HOG (only HOGs carrying >=1 gene in the panel):
   HOG_ID, n_genes, n_effector, n_cytoplasmic, n_apoplastic, n_both, EffectorHOG (>=1 effector gene), DominantClass
3. effector_master.coverage.tsv    - sample-level coverage summary (which samples lack effectorp.tsv)

Effector definition follows EffectorP v3 consensus Prediction column:
  Non-effector                       -> not an effector
  Cytoplasmic effector / Apoplastic effector / Cytoplasmic/apoplastic effector
  / Apoplastic/cytoplasmic effector  -> effector (class = both if both terms present)

Usage:
    python3 build_effector_master.py \
        --sample-list sample.solaniSC.211.list \
        --per-isolate-dir 02.annotation/per_isolates \
        --gene-vs-hog orthofinder-shorter/N0.gene_vs_HOG.tsv \
        --output-dir 08.pan_020_genes
"""

import argparse
import os
import sys

import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Build master effector table for FSSC isolates")
    parser.add_argument("--sample-list", required=True, help="One sample ID per line (211 samples)")
    parser.add_argument("--per-isolate-dir", required=True,
                        help="Directory holding per_isolates/<sample>/signalp/effectorp.tsv")
    parser.add_argument("--gene-vs-hog", required=True,
                        help="N0.gene_vs_HOG.tsv: GeneID SampleID HOG_ID OG_ID")
    parser.add_argument("--output-dir", required=True, help="Output directory for tables")
    return parser.parse_args()


def classify(prediction: str) -> str:
    """Map EffectorP consensus prediction to a class label."""
    low = prediction.lower()
    if low == "non-effector":
        return "non-effector"
    cyto = "cytoplasmic" in low
    apo = "apoplastic" in low
    if cyto and apo:
        return "both"
    if cyto:
        return "cytoplasmic"
    if apo:
        return "apoplastic"
    return "unknown"


def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.sample_list) as fh:
        samples = [line.strip() for line in fh if line.strip()]

    frames = []
    missing = []
    for sample in samples:
        path = os.path.join(args.per_isolate_dir, sample, "signalp", "effectorp.tsv")
        if not os.path.isfile(path):
            missing.append(sample)
            continue
        # effectorp.tsv header line starts with '#' and the first column holds
        # "<GeneID> <short-id>"; read positionally and take the first token as GeneID.
        df = pd.read_csv(path, sep="\t", header=None, skiprows=1, comment="#",
                         names=["Identifier", "Cytoplasmic_score", "Apoplastic_score",
                                "NonEffector_score", "Prediction"])
        df["GeneID"] = df["Identifier"].str.split().str[0]
        df["SampleID"] = sample
        df["EffectorClass"] = df["Prediction"].map(classify)
        df["IsEffector"] = df["EffectorClass"] != "non-effector"
        frames.append(df[["GeneID", "SampleID", "Cytoplasmic_score", "Apoplastic_score",
                          "NonEffector_score", "Prediction", "EffectorClass", "IsEffector"]])

    if not frames:
        sys.exit(f"ERROR: no effectorp.tsv files found under {args.per_isolate_dir}")
    gene_df = pd.concat(frames, ignore_index=True)
    gene_df = gene_df.drop_duplicates(subset=["GeneID"]).reset_index(drop=True)

    gene_out = os.path.join(args.output_dir, "effector_master.gene_level.tsv")
    gene_df.to_csv(gene_out, sep="\t", index=False)
    n_eff = int(gene_df["IsEffector"].sum())
    print(f"gene-level: {len(gene_df)} genes, {n_eff} effectors "
          f"({n_eff / len(gene_df) * 100:.2f}%) from {gene_df['SampleID'].nunique()} samples")
    print(f"  -> {gene_out}")

    # HOG-level table (only HOGs present in the gene-vs-HOG mapping)
    gvh = pd.read_csv(args.gene_vs_hog, sep="\t",
                      usecols=["GeneID", "SampleID", "HOG_ID"])
    gvh = gvh.drop_duplicates(subset=["GeneID"])
    merged = gene_df.merge(gvh[["GeneID", "HOG_ID"]], on="GeneID", how="left")

    hog_stats = (
        merged.groupby("HOG_ID")
        .agg(
            n_genes=("GeneID", "count"),
            n_effector=("IsEffector", "sum"),
            n_cytoplasmic=("EffectorClass", lambda s: s.eq("cytoplasmic").sum() + s.eq("both").sum()),
            n_apoplastic=("EffectorClass", lambda s: s.eq("apoplastic").sum() + s.eq("both").sum()),
            n_both=("EffectorClass", lambda s: s.eq("both").sum()),
        )
        .reset_index()
    )
    hog_stats["EffectorHOG"] = hog_stats["n_effector"] > 0
    hog_stats["DominantClass"] = (
        hog_stats[["n_cytoplasmic", "n_apoplastic"]]
        .idxmax(axis=1)
        .map({"n_cytoplasmic": "cytoplasmic", "n_apoplastic": "apoplastic"})
    )
    hog_stats.loc[hog_stats["n_effector"] == 0, "DominantClass"] = "non-effector"

    # drop HOGs not in the annotation/interval universe
    hog_stats = hog_stats[hog_stats["HOG_ID"].str.startswith("N0.HOG", na=False)]

    hog_out = os.path.join(args.output_dir, "effector_master.hog_level.tsv")
    hog_stats.to_csv(hog_out, sep="\t", index=False)
    n_hog_eff = int(hog_stats["EffectorHOG"].sum())
    print(f"hog-level: {len(hog_stats)} HOGs, {n_hog_eff} effector HOGs "
          f"({n_hog_eff / len(hog_stats) * 100:.2f}%)")
    print(f"  -> {hog_out}")

    # Coverage summary
    cov = pd.DataFrame({"Sample": samples})
    cov["HasEffectorP"] = ~cov["Sample"].isin(missing)
    cov_out = os.path.join(args.output_dir, "effector_master.coverage.tsv")
    cov.to_csv(cov_out, sep="\t", index=False)
    print(f"coverage: {cov['HasEffectorP'].sum()}/{len(cov)} samples have effectorp.tsv")
    if missing:
        print("  missing samples:", ", ".join(missing))
    print(f"  -> {cov_out}")


if __name__ == "__main__":
    main()
