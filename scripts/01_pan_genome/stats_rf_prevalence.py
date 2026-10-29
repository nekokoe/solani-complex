#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prevalence and functional-composition statistics for the RF-selected
clinical-associated HOGs, per species (F. solani / F. falciformis).

For each selected HOG:
  clinical_rate    = fraction of clinical isolates of the species cohort carrying it
  nonclinical_rate = fraction of non-clinical isolates carrying it
  dominance        = clinical / non-clinical / equal
Functional types per HOG (HOG-level):
  effector    - HOG contains >=1 EffectorP-3 predicted effector gene
  kinase      - HOG has a protein-kinase Pfam domain
  transporter - HOG's primary Pfam class is Transporter
  regulatory  - HOG's primary Pfam class is Transcription factor

Outputs:
  suppl_table.rf_prevalence.tsv      per-HOG rates + dominance
  suppl_table.rf_dom_functions.tsv   functional composition by species x dominance
  figure_data.rf_v2.tsv              pie + summary numbers for the figure

Usage:
    python3 stats_rf_prevalence.py \
        --rf-solani .../hog_count.solani/top_significant_features.txt \
        --rf-falciformis .../hog_count.falciformis/top_significant_features.txt \
        --hog-count 05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv \
        --group-file bmu_ref_groups.211.txt \
        --solani-cohort .../all.solani.sample.list \
        --falciformis-cohort .../all.falciformis.sample.list \
        --effector-hog 08.pan_020_genes/effector_master.hog_level.tsv \
        --annotation-file 05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv \
        --output-dir 09.rf_clinical_genes
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_020_genes import classify_pfam  # noqa: E402

N_TOTAL_HOGS = 29080


def parse_arguments():
    p = argparse.ArgumentParser(description="RF-selected HOG prevalence statistics")
    p.add_argument("--rf-solani", required=True)
    p.add_argument("--rf-falciformis", required=True)
    p.add_argument("--hog-count", required=True)
    p.add_argument("--group-file", required=True)
    p.add_argument("--solani-cohort", required=True)
    p.add_argument("--falciformis-cohort", required=True)
    p.add_argument("--effector-hog", required=True)
    p.add_argument("--annotation-file", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    groups = pd.read_csv(args.group_file, sep="\t")
    groups.columns = [c.strip() for c in groups.columns]
    med_env = dict(zip(groups[groups.columns[0]], groups["MedEnv"]))

    hc = pd.read_csv(args.hog_count, sep="\t")
    hc = hc[pd.to_numeric(hc["GeneCount"], errors="coerce").fillna(0) > 0]
    carriers = hc.groupby("HOG_ID")["SampleID"].apply(set)

    eff = pd.read_csv(args.effector_hog, sep="\t", usecols=["HOG_ID", "EffectorHOG"])
    eff_set = set(eff.loc[eff["EffectorHOG"], "HOG_ID"])
    ann = pd.read_csv(args.annotation_file, sep="\t",
                      usecols=["HOG_ID", "Pfam_Domain"])
    ann = ann.set_index("HOG_ID")
    ann = ann[~ann.index.duplicated(keep="first")]
    kin_set = set(ann.index[ann["Pfam_Domain"].str.contains(
        "pkinase|hiska", case=False, regex=True, na=False)])
    pfam_class = ann["Pfam_Domain"].dropna().apply(
        lambda v: classify_pfam(str(v).split(";")))
    trans_set = set(pfam_class[pfam_class == "Transporter"].index)
    regul_set = set(pfam_class[pfam_class == "Transcription factor"].index)

    rows, fig_rows = [], []
    for sp, feat_path, cohort_path in [("solani", args.rf_solani, args.solani_cohort),
                                       ("falciformis", args.rf_falciformis, args.falciformis_cohort)]:
        feat = pd.read_csv(feat_path, sep="\t")
        # the count screen calls the id column Feature, the PAV screen HOG_ID
        id_col = "Feature" if "Feature" in feat.columns else "HOG_ID"
        top = feat[id_col].tolist()
        cohort = [l.strip() for l in open(cohort_path) if l.strip()]
        clin = [s for s in cohort if med_env.get(s) == "med"]
        nonc = [s for s in cohort if med_env.get(s) == "env"]

        def rate(h, iso):
            return len(carriers.get(h, set()) & set(iso)) / len(iso)

        for h in top:
            cr, nr = rate(h, clin), rate(h, nonc)
            dom = "clinical" if cr > nr else ("non-clinical" if cr < nr else "equal")
            rows.append({
                "Species": sp, "HOG_ID": h,
                "Clinical_rate": cr, "NonClinical_rate": nr, "Dominance": dom,
                "Effector": h in eff_set,
                "Kinase": h in kin_set,
                "Transporter": h in trans_set,
                "Regulatory": h in regul_set,
            })
        n_top = len(top)
        fig_rows.append({"panel": "pie", "species": sp, "category": "selected",
                         "value": n_top / N_TOTAL_HOGS * 100, "n": n_top})
        fig_rows.append({"panel": "pie", "species": sp, "category": "rest",
                         "value": (N_TOTAL_HOGS - n_top) / N_TOTAL_HOGS * 100,
                         "n": N_TOTAL_HOGS - n_top})
        fig_rows.append({"panel": "pie", "species": sp, "category": "total",
                         "value": 100.0, "n": N_TOTAL_HOGS})

    prev = pd.DataFrame(rows)
    prev_out = os.path.join(args.output_dir, "suppl_table.rf_prevalence.tsv")
    prev.to_csv(prev_out, sep="\t", index=False)
    print(f"prevalence -> {prev_out}")

    # functional composition by species x dominance
    func_cols = ["Effector", "Kinase", "Transporter", "Regulatory"]
    dom_rows = []
    for sp in ["solani", "falciformis"]:
        for dom in ["clinical", "non-clinical"]:
            sub = prev[(prev["Species"] == sp) & (prev["Dominance"] == dom)]
            n = len(sub)
            for fc in func_cols:
                k = int(sub[fc].sum())
                dom_rows.append({"Species": sp, "Dominance": dom, "Function": fc,
                                 "n": n, "n_Function": k,
                                 "pct": k / n * 100 if n else np.nan})
            # Others = HOGs belonging to none of the four functional classes
            other = ~(sub["Effector"] | sub["Kinase"] | sub["Transporter"]
                      | sub["Regulatory"])
            k_other = int(other.sum())
            dom_rows.append({"Species": sp, "Dominance": dom, "Function": "Others",
                             "n": n, "n_Function": k_other,
                             "pct": k_other / n * 100 if n else np.nan})
    dom_df = pd.DataFrame(dom_rows)
    dom_out = os.path.join(args.output_dir, "suppl_table.rf_dom_functions.tsv")
    dom_df.to_csv(dom_out, sep="\t", index=False)
    print(f"dominance functions -> {dom_out}")
    print(dom_df.pivot_table(index=["Species", "Dominance"], columns="Function",
                             values="pct").round(1).to_string())

    fig_df = pd.DataFrame(fig_rows)
    fig_out = os.path.join(args.output_dir, "figure_data.rf_v2.tsv")
    fig_df.to_csv(fig_out, sep="\t", index=False)
    print(f"figure data -> {fig_out}")


if __name__ == "__main__":
    main()
