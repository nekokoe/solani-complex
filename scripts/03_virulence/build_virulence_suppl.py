#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_virulence_suppl.py — regenerate the virulence supplementary workbook.

`00.manuscript/supplementary/suppl_table.S19_virulence_assoc.xlsx` (Supplementary Table
S19) had no generator
script in the repository; its three association sheets are reconstructions of the
mapping used before (verified row-for-row against the archived workbook before the
21-strain rerun), and the two enrichment sheets now come from the corrected
hypergeometric test written by scripts/derive_virulence_modes.py.

Sheets:
  Significant_HOGs_S19   HOGs whose best per-temperature association reaches P<0.001
  Nominal_P005_HOGs      HOGs with any concentration trait at P<0.05, with the
                         per-temperature mean |rho| and a 28C/37C/shared label
  Temp_effect_HOGs       HOGs associated with the temperature-preference trait (P<0.05)
  sp37_toolkit_enrichment  Pfam enrichment of the 37 C-biased toolkit (hog_sp37.tsv, 15 HOGs)
  sp28_toolkit_enrichment  the same for the 28 C-biased toolkit (hog_sp28.tsv, 18 HOGs)

Usage: python3 scripts/build_virulence_suppl.py [--out PATH]
"""
import argparse
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
ASSOC = "07.lethal_groups/genomewide_assoc/hog_phenotype_assoc_auc_no10e5.tsv"
ANNOT = "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv"
ENR = "07.lethal_groups/genomewide_assoc/enrich_{tag}_mode.tsv"
DEFAULT_OUT = "00.manuscript/supplementary/suppl_table.S19_virulence_assoc.xlsx"
# Phenotypes: 1 - AUC at 10^3 and 10^4 plus the two-dose means and their difference.
CONC = ["vir_28Ce3", "vir_28Ce4", "vir_37Ce3", "vir_37Ce4"]
T28 = ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"]
T37 = ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"]
TE = "vir_temp_effect2"
MODE_GAP = 0.25


def best(sub, traits):
    """rho / p / trait at the strongest |rho| within `traits`."""
    g = sub[sub.trait.isin(traits)].groupby("HOG_ID")
    return pd.DataFrame({
        "rho": g.apply(lambda x: x.loc[x["rho"].abs().idxmax(), "rho"], include_groups=False),
        "p": g.apply(lambda x: x.loc[x["rho"].abs().idxmax(), "p"], include_groups=False),
        "trait": g.apply(lambda x: x.loc[x["rho"].abs().idxmax(), "trait"], include_groups=False)})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    a = pd.read_csv(ASSOC, sep="\t")
    # every trait in this table was tested (the saturated 10^5 rates are not in it)
    ann = pd.read_csv(ANNOT, sep="\t").set_index("HOG_ID")
    b28, b37 = best(a, T28), best(a, T37)
    mean28 = a[a.trait.isin(T28)].groupby("HOG_ID").rho.apply(lambda s: s.abs().mean())
    mean37 = a[a.trait.isin(T37)].groupby("HOG_ID").rho.apply(lambda s: s.abs().mean())

    def label(rho28, rho37):
        d = abs(rho28) - abs(rho37)
        return "28C" if d >= MODE_GAP else ("37C" if d <= -MODE_GAP else "shared")

    def annotate(rr):
        out = rr.copy()
        out["gene"] = [str(ann.Gene_Name.get(h, "")) if str(ann.Gene_Name.get(h, "")) not in ("nan",) else ""
                       for h in out.index]
        out["pfam"] = [str(ann.Pfam_Domain.get(h, "")) if str(ann.Pfam_Domain.get(h, "")) not in ("nan",) else ""
                       for h in out.index]
        out["nog"] = [str(ann.NOG_Description.get(h, "")) if str(ann.NOG_Description.get(h, "")) not in ("nan",) else ""
                      for h in out.index]
        return out

    # ---- Significant_HOGs_S24 ----
    sig = pd.DataFrame({"rho28_best": b28.rho, "p28_best": b28.p, "trait28": b28.trait,
                        "rho37_best": b37.rho, "p37_best": b37.p, "trait37": b37.trait})
    sig["p_best"] = sig[["p28_best", "p37_best"]].min(axis=1)
    sig = annotate(sig[sig.p_best < 0.001])
    sig["mode"] = [label(r.rho28_best, r.rho37_best) for r in sig.itertuples()]
    sig["spec"] = sig["mode"]
    sig = sig.reset_index().rename(columns={"HOG_ID": "HOG"})
    cols = ["HOG", "mode", "spec", "rho28_best", "p28_best", "trait28",
            "rho37_best", "p37_best", "trait37", "gene", "pfam", "nog", "p_best"]
    sig = sig[cols].sort_values("p_best")

    # ---- Nominal_P005_HOGs ----
    conc = a[a.trait.isin(CONC)]
    pmin = conc.groupby("HOG_ID").p.min()
    nom = pd.DataFrame({"p_min_concentration": pmin,
                        "mean_rho_28C": mean28, "mean_rho_37C": mean37})
    nom = nom[nom.p_min_concentration < 0.05]
    nom["mode"] = [label(nom.mean_rho_28C[h], nom.mean_rho_37C[h]) for h in nom.index]
    nom = annotate(nom).reset_index().rename(columns={"HOG_ID": "HOG"})
    nom = nom[["HOG", "p_min_concentration", "mean_rho_28C", "mean_rho_37C", "mode",
               "gene", "pfam", "nog"]].sort_values("p_min_concentration")

    # ---- Temp_effect_HOGs ----
    te = a[a.trait == TE][["HOG_ID", "rho", "p"]].set_index("HOG_ID")
    te = te[te.p < 0.05]
    te["highlight"] = ""
    if len(te):
        for h in te.rho.abs().sort_values(ascending=False).head(2).index:
            te.loc[h, "highlight"] = "Fig7c-candidate"
    te = annotate(te).reset_index().rename(columns={"HOG_ID": "HOG", "rho": "rho_temp_effect",
                                                    "p": "p_temp_effect"})
    te = te[["HOG", "rho_temp_effect", "p_temp_effect", "highlight", "gene", "pfam", "nog"]]
    te = te.sort_values("p_temp_effect")

    # Enrichment of the CMH-defined toolkits (the 2 + 2 grouping), recomputed here so the
    # workbook does not depend on the retired derive_virulence_modes.py output.
    from math import comb
    from scipy import stats as _st
    _annot = pd.read_csv(ANNOT, sep="\t", dtype=str).drop_duplicates("HOG_ID").set_index("HOG_ID")
    _uni = list(a.HOG_ID.unique())

    def _doms(h):
        v = _annot.at[h, "Pfam_Domain"] if h in _annot.index else ""
        return [d for d in str(v).split(";") if d and d.lower() not in ("nan", "-")]

    _bg = pd.Series([d for h in _uni for d in set(_doms(h))]).value_counts()

    def _enrich(hog_file):
        ids = [l.split("\t")[0] for l in open(hog_file).read().splitlines()[1:] if l.strip()]
        ids = [h for h in ids if h in set(_uni)]
        K = len(ids)
        cnt = pd.Series([d for h in ids for d in set(_doms(h))]).value_counts()
        rows = [{"domain": d, "in_set": int(c), "background": int(_bg.get(d, 0)), "set_size": K,
                 "p": _st.hypergeom.sf(int(c) - 1, len(_uni), int(_bg.get(d, 0)), K)}
                for d, c in cnt.items()]
        e = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
        m = len(e)
        e["q"] = np.minimum(1, (e.p * m / (e.index + 1))[::-1].cummin()[::-1])
        return e

    e28 = _enrich("07.lethal_groups/genomewide_assoc/hog_sp37.tsv")
    e37 = _enrich("07.lethal_groups/genomewide_assoc/hog_sp28.tsv")

    out = os.path.abspath(args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with pd.ExcelWriter(out) as w:
        sig.to_excel(w, sheet_name="Significant_HOGs_S19", index=False)
        nom.to_excel(w, sheet_name="Nominal_P005_HOGs", index=False)
        te.to_excel(w, sheet_name="Temp_effect_HOGs", index=False)
        e28.to_excel(w, sheet_name="sp37_toolkit_enrichment", index=False)
        e37.to_excel(w, sheet_name="sp28_toolkit_enrichment", index=False)
    print(f"Significant_HOGs_S19 {len(sig)}   modes: {sig['mode'].value_counts().to_dict()}")
    print(f"Nominal_P005_HOGs    {len(nom)}   modes: {nom['mode'].value_counts().to_dict()}")
    print(f"Temp_effect_HOGs     {len(te)}")
    print(f"28C_enrichment       {len(e28)}   q_min {e28.q.min() if len(e28) else float('nan'):.3f}")
    print(f"37C_enrichment       {len(e37)}   q_min {e37.q.min() if len(e37) else float('nan'):.3f}")
    print(f"wrote {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
