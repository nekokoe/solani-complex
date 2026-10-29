#!/usr/bin/env python3
"""Calibrate the S7 module threshold against a permutation null, and report the sweep.

The module rule has always been "a variable orthogroup whose max |Spearman rho| against
that temperature's phenotypes reaches a threshold", with the threshold set by convention
(0.50 in the script, 0.55 in the Methods). Neither value was calibrated. This script
calibrates it: permute the trait labels, take max|rho| over the 4,484 variable orthogroups,
and use the quantiles of that null as candidate thresholds — a member then lies beyond what
chance produces for the most extreme orthogroup.

It also reports what the enrichment test does at each candidate threshold, because the
answer turns out to depend on it: a Pfam enrichment appears only inside a narrow band and
is carried by one or two orthogroups.

Writes
    07.lethal_groups/genomewide_assoc/S7_threshold_sensitivity.tsv
    and appends that table to suppl_table.S7_temperature_modules.xlsx as
    `Threshold_sensitivity`.
"""

import os
from collections import Counter

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GW = os.path.join(BASE, "07.lethal_groups/genomewide_assoc")
PA = os.path.join(BASE, "07.lethal_groups/phenotype_analysis")
XLSX = os.path.join(BASE, "00.manuscript/supplementary",
                    "suppl_table.S7_temperature_modules.xlsx")
HOGCOUNT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
ANNOT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
QUANTILES = (0.50, 0.90, 0.95, 0.99, 0.995)
N_PERM = 500


def main():
    traits_df = pd.read_csv(f"{PA}/phenotype_metrics_auc_no10e5.tsv", sep="\t", index_col=0)
    strains = list(traits_df.index)
    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = (hog[hog.SampleID.isin(strains)].pivot(index="HOG_ID", columns="SampleID",
                                                 values="GeneCount")
           .reindex(columns=strains).fillna(0))
    keep = mat[((mat > 0).sum(1) >= 3) & ((mat > 0).sum(1) <= len(strains) - 2)
               & (mat[strains].nunique(1) > 1)]
    ids, N = list(keep.index), len(keep)
    R = np.apply_along_axis(rankdata, 1, keep[strains].to_numpy().astype(float))
    Rc = R - R.mean(1, keepdims=True)

    def rho(y):
        ry = rankdata(y)
        yc = ry - ry.mean()
        return np.clip((Rc * yc).sum(1) / np.sqrt((Rc ** 2).sum(1) * (yc ** 2).sum()), -1, 1)

    TL = {"28C": ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"],
          "37C": ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"]}
    assoc = pd.read_csv(f"{GW}/hog_phenotype_assoc_auc_no10e5.tsv", sep="\t")
    mx = {t: assoc[assoc.trait.isin(tl)].groupby("HOG_ID")["rho"].apply(lambda x: x.abs().max())
          for t, tl in TL.items()}

    rng = np.random.default_rng(7)
    null = {t: np.array([np.max([np.abs(rho(rng.permutation(traits_df[c].values))) for c in tl],
                                 axis=0) for _ in range(N_PERM)]) for t, tl in TL.items()}

    annot = pd.read_csv(ANNOT, sep="\t", dtype=str).drop_duplicates("HOG_ID").set_index("HOG_ID")
    dom = {h: set(x for x in str(annot.at[h, "Pfam_Domain"] if h in annot.index else "").split(";")
                  if x and x.lower() not in ("nan", "-")) for h in ids}
    bg = Counter(x for h in ids for x in dom[h])

    def bh(p):
        p = np.asarray(p, float)
        n = len(p)
        o = np.argsort(p)
        r = p[o] * n / (np.arange(n) + 1)
        q = np.empty(n)
        q[o] = np.minimum(1, np.minimum.accumulate(r[::-1])[::-1])
        return q

    rows = []
    print(f"{'quantile':>9s} {'temp':>4s} {'thr':>6s} {'module':>7s} {'best domain':>18s} "
          f"{'k/K':>8s} {'q':>7s} {'q<0.05':>7s}")
    for qq in QUANTILES:
        for tag in ("28C", "37C"):
            thr = float(np.quantile(null[tag], qq))
            mod = sorted(set(mx[tag][mx[tag] >= thr].index) & set(ids))
            if not mod:
                continue
            sc = Counter(x for h in mod for x in dom[h])
            pr = [(d, c, bg.get(d, 0), stats.hypergeom.sf(c - 1, N, bg.get(d, 0), len(mod)))
                  for d, c in sc.items()]
            qs = bh([x[3] for x in pr])
            i = int(np.argmin(qs))
            rows.append(dict(quantile=qq, temperature=tag, rho_threshold=round(thr, 4),
                             module_size=len(mod), best_domain=pr[i][0],
                             best_k=pr[i][1], background=pr[i][2],
                             best_p=pr[i][3], best_q=qs[i],
                             n_passing_q05=int((qs < 0.05).sum())))
            r = rows[-1]
            print(f"{qq:9.3f} {tag:>4s} {thr:6.3f} {len(mod):7d} {pr[i][0]:>18s} "
                  f"{pr[i][1]:3d}/{len(mod):<4d} {qs[i]:7.3f} {int((qs < 0.05).sum()):7d}")

    out = pd.DataFrame(rows)
    path = f"{GW}/S7_threshold_sensitivity.tsv"
    out.to_csv(path, sep="\t", index=False)

    wb = openpyxl.load_workbook(XLSX)
    if "Threshold_sensitivity" in wb.sheetnames:
        del wb["Threshold_sensitivity"]
    ws = wb.create_sheet("Threshold_sensitivity")
    ws.append(["item", "value"])
    ws.append(["threshold rule",
               "module member = max|Spearman rho| against that temperature's three AUC "
               "phenotypes >= the quantile of the permutation null"])
    ws.append(["null construction",
               f"trait labels permuted {N_PERM} times; max|rho| taken over the 4,484 variable "
               "orthogroups, per temperature"])
    ws.append(["note",
               "no single threshold is privileged: a Pfam enrichment appears only inside a "
               "narrow band and is carried by one or two orthogroups"])
    ws.append([])
    ws.append(list(out.columns))
    for r in out.itertuples(index=False):
        ws.append(list(r))
    wb.save(XLSX)
    print(f"\nwrote {path} and sheet 'Threshold_sensitivity' in {os.path.basename(XLSX)}")


if __name__ == "__main__":
    main()
