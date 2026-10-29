#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functional comparison of A-A (accessory-accessory) and C-A (core-accessory)
transfer-involved genes.

Transfer-involved genes come from the 211-genome clustering analysis
(gene_transfer_211/suppl_table.gene_transfer_211.tsv). Each gene is annotated
with:
  - Pfam primary functional class (classify_pfam) and protein-kinase status
  - KEGG BRITE level-2 category (eggNOG ko -> kegg_pathway_categories)
  - CAZyme status + major class (dbCAN)
  - Effector status/class (EffectorP, effector_master.gene_level.tsv)

For every functional category we test A-A vs C-A with a two-sided Fisher
exact test and BH-FDR correction.

Outputs (06.accessory-chrs/gene_transfer_211/):
  suppl_table.transfer_function.tsv         per-category comparison
  suppl_table.transfer_gene_functions.tsv   per-gene functional flags
  figure_transfer_function.pdf              top enriched / depleted categories
"""

import argparse
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

import stats_020_genes as s020

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]


def parse_args():
    p = argparse.ArgumentParser(description="A-A vs C-A functional comparison")
    p.add_argument("--gene-table", required=True)
    p.add_argument("--annotation-root", required=True)
    p.add_argument("--effector-master", required=True)
    p.add_argument("--kegg-categories", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_kegg_level2(path):
    """ko id -> level2 category."""
    m = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                m[parts[0]] = parts[2].strip()
    return m


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    df = pd.read_csv(args.gene_table, sep="\t")
    aa = df[(df["AC_related_transfer"]) & (df["Transfer_type"] == "A-A")]
    ca = df[(df["AC_related_transfer"]) & (df["Transfer_type"] == "C-A")]
    print(f"A-A genes: {len(aa)}; C-A genes: {len(ca)}")

    kegg2 = load_kegg_level2(args.kegg_categories)

    # ---- functional annotations per gene ----
    pfam_class, kinase = {}, {}
    cazy, cazy_class = {}, {}
    kegg_cat = {}
    for iso in ISOLATES:
        base = os.path.join(args.annotation_root, iso)
        pf = {}
        with open(os.path.join(base, f"{iso}.gene_pfam.tsv")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    pf.setdefault(parts[0], []).append(parts[2])
        for g, doms in pf.items():
            pfam_class[g] = s020.classify_pfam(doms)
            kinase[g] = any(d.lower().startswith(("pkinase", "hiska"))
                            or "kinase" in d.lower() for d in doms)
        db = defaultdict(list)
        with open(os.path.join(base, f"{iso}.gene_dbcan.tsv")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    db[parts[0]].append(parts[1])
        for g, fams in db.items():
            cazy[g] = True
            cls = set(f[:2] for f in fams if len(f) >= 2 and f[:2].isalpha())
            cazy_class[g] = ";".join(sorted(cls)) if cls else "NA"
        eg = {}
        with open(os.path.join(base, f"{iso}.gene_eggnog.tsv")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5 and parts[4].startswith("ko:"):
                    eg[parts[0]] = parts[4][3:]
        for g, ko in eg.items():
            kegg_cat[g] = kegg2.get(ko)

    eff = pd.read_csv(args.effector_master, sep="\t")
    eff_map = dict(zip(eff["GeneID"], eff["IsEffector"]))
    eff_cls = dict(zip(eff["GeneID"], eff["EffectorClass"]))

    # ---- annotate each transfer gene ----
    def annotate(grp):
        rows = []
        for _, r in grp.iterrows():
            g = r["Transcript"]
            rows.append({
                "Isolate": r["Isolate"], "Transcript": g,
                "Transfer_type": r["Transfer_type"],
                "Contig_type": r["Contig_type"],
                "Pfam_class": pfam_class.get(g, "No Pfam"),
                "Kinase": kinase.get(g, False),
                "CAZyme": cazy.get(g, False),
                "CAZy_class": cazy_class.get(g, "NA"),
                "KEGG_level2": kegg_cat.get(g, "No KEGG"),
                "Effector": bool(eff_map.get(g, False)),
                "EffectorClass": eff_cls.get(g, "no-annotation"),
            })
        return pd.DataFrame(rows)

    aa_a = annotate(aa)
    ca_a = annotate(ca)
    all_a = pd.concat([aa_a, ca_a])
    all_a.to_csv(os.path.join(args.output_dir,
                              "suppl_table.transfer_gene_functions.tsv"),
                 sep="\t", index=False)

    # ---- per-category comparison ----
    rows = []
    def test(col, value, tag):
        a = (aa_a[col] == value).sum()
        b = (ca_a[col] == value).sum()
        if a == 0 and b == 0:
            return
        or_, p = stats.fisher_exact([[a, len(aa_a) - a], [b, len(ca_a) - b]])
        rows.append({"Category": tag, "A-A_frac": a / len(aa_a),
                     "C-A_frac": b / len(ca_a), "n_A-A": a, "n_C-A": b,
                     "Odds_ratio": or_, "p": p})

    for cls in sorted({v for v in pfam_class.values()}):
        test("Pfam_class", cls, f"Pfam_{cls}")
    test("Kinase", True, "Kinase")
    test("CAZyme", True, "CAZyme")
    for c in sorted({v for v in cazy_class.values() if v != "NA"}):
        test("CAZy_class", c, f"CAZy_{c}")
    for k in sorted({v for v in kegg_cat.values() if v}):
        test("KEGG_level2", k, f"KEGG_{k}")
    test("Effector", True, "Effector")
    for c in ["cytoplasmic", "apoplastic", "both"]:
        test("EffectorClass", c, f"Effector_{c}")

    comp = pd.DataFrame(rows)
    comp["FDR"] = stats.false_discovery_control(comp["p"], method="bh")
    comp = comp.sort_values("p")
    comp.to_csv(os.path.join(args.output_dir,
                             "suppl_table.transfer_function.tsv"),
                sep="\t", index=False)
    print("\n=== A-A vs C-A functional comparison (sorted by p) ===")
    print(comp.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
