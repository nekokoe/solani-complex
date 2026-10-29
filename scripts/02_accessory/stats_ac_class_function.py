#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functional comparison across AC transfer classes (13 TGS assemblies only).

Q2: genes transferred by Type-1 (block-transfer) ACs vs genes transferred by
    Type-2 (scattered/single-gene) ACs.
Q3: genes in cross-species transfer clusters vs all AC-resident transfer genes.

Annotations: Pfam functional class (classify_pfam), protein kinase, CAZyme,
KEGG BRITE level-2, effector (EffectorP).

Outputs (06.accessory-chrs/gene_transfer_211/):
  suppl_table.ac_class_function.tsv        Q2 comparison
  suppl_table.cross_species_function.tsv   Q3 comparison
"""

import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

import stats_020_genes as s020

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]


def parse_args():
    p = argparse.ArgumentParser(description="AC class / cross-species function")
    p.add_argument("--gene-class-table", required=True)
    p.add_argument("--annotation-root", required=True)
    p.add_argument("--effector-master", required=True)
    p.add_argument("--kegg-categories", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def load_kegg_level2(path):
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
    kegg2 = load_kegg_level2(args.kegg_categories)

    # functional annotations
    pfam_class, kinase, cazy, cazy_class, kegg_cat = {}, {}, {}, {}, {}
    for iso in ISOLATES:
        base = os.path.join(args.annotation_root, iso)
        pf = defaultdict(list)
        with open(os.path.join(base, f"{iso}.gene_pfam.tsv")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    pf[parts[0]].append(parts[2])
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
            cls = {f[:2] for f in fams if len(f) >= 2 and f[:2].isalpha()}
            cazy_class[g] = ";".join(sorted(cls)) if cls else "NA"
        with open(os.path.join(base, f"{iso}.gene_eggnog.tsv")) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5 and parts[4].startswith("ko:"):
                    kegg_cat[parts[0]] = kegg2.get(parts[4][3:])

    eff = pd.read_csv(args.effector_master, sep="\t")
    eff_map = dict(zip(eff["GeneID"], eff["IsEffector"]))
    eff_cls = dict(zip(eff["GeneID"], eff["EffectorClass"]))

    g = pd.read_csv(args.gene_class_table, sep="\t")
    g["Pfam_class"] = g["Transcript"].map(lambda x: pfam_class.get(x, "No Pfam"))
    g["Kinase"] = g["Transcript"].map(lambda x: kinase.get(x, False))
    g["CAZyme"] = g["Transcript"].map(lambda x: cazy.get(x, False))
    g["CAZy_class"] = g["Transcript"].map(lambda x: cazy_class.get(x, "NA"))
    g["KEGG_level2"] = g["Transcript"].map(lambda x: kegg_cat.get(x, "No KEGG"))
    g["Effector"] = g["Transcript"].map(lambda x: bool(eff_map.get(x, False)))

    def compare(test, control, out_name):
        rows = []
        def t(col, value, tag):
            a = (test[col] == value).sum()
            b = (control[col] == value).sum()
            if a == 0 and b == 0:
                return
            or_, p = stats.fisher_exact([[a, len(test) - a],
                                         [b, len(control) - b]])
            rows.append({"Category": tag, "test_frac": a / len(test),
                         "control_frac": b / len(control),
                         "n_test": a, "n_control": b,
                         "Odds_ratio": or_, "p": p})
        for cls in sorted({v for v in pfam_class.values()}):
            t("Pfam_class", cls, f"Pfam_{cls}")
        t("Kinase", True, "Kinase")
        t("CAZyme", True, "CAZyme")
        for k in sorted({v for v in kegg_cat.values() if v}):
            t("KEGG_level2", k, f"KEGG_{k}")
        t("Effector", True, "Effector")
        comp = pd.DataFrame(rows)
        comp["FDR"] = stats.false_discovery_control(comp["p"], method="bh")
        comp = comp.sort_values("p")
        comp.to_csv(os.path.join(args.output_dir, out_name), sep="\t",
                    index=False)
        print(f"\n=== {out_name} (test n={len(test)}, control n={len(control)}) ===")
        print(comp.head(18).round(4).to_string(index=False))

    # Q2: Type-1 block vs Type-2 scattered
    t1 = g[g["AC_class"] == "Type1_block"]
    t2 = g[g["AC_class"] == "Type2_scattered"]
    compare(t1, t2, "suppl_table.ac_class_function.tsv")

    # Q3: cross-species vs all AC-resident transfer genes
    xs = g[g["Cross_species"]]
    compare(xs, g, "suppl_table.cross_species_function.tsv")


if __name__ == "__main__":
    main()
