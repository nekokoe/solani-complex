#!/usr/bin/env python3
"""Strain-specific ("toolkit") HOG sets for the temperature-biased groups.

The grouping comes from the per-strain Cochran-Mantel-Haenszel test in
`cmh_temperature_preference.py` (P < 0.05, stratified by inoculum dose), not from a
hard cutoff on the continuous difference in mean mortality. Under that test
B37 = {BMU01287, BMU02031} and B28 = {BMU03976, BMU12685}.

Sets, unchanged from the definition used for the previous 15/1 report:
    keep  variable HOGs carried by 3-19 of the 21 strains (the same 4,484 HOG universe
          the genome-wide association screen uses)
    sp37  carried by every 37 °C-biased strain, by <= 3 strains among B28 + neutral
    sp28  carried by every 28 °C-biased strain, by <= 4 strains among B37 + neutral

Output
    07.lethal_groups/genomewide_assoc/hog_sp37.tsv        (copy-number matrix)
    07.lethal_groups/genomewide_assoc/hog_sp28.tsv
    07.lethal_groups/genomewide_assoc/toolkit_hog_sets_summary.tsv
"""

import csv
import os
from math import comb

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PHENO = os.path.join(ROOT, "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv")
BIAS = os.path.join(ROOT, "07.lethal_groups/phenotype_analysis/temperature_bias_cmh.tsv")
HOGCOUNT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
ANNOT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
OUTDIR = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc")
CARRY_MIN, CARRY_MAX = 3, 19
NEAR_37, NEAR_28 = 3, 4


def fisher_two_sided(a, b, c, d):
    """Exact two-sided Fisher test on [[a, b], [c, d]]."""
    n = a + b + c + d
    row1, col1 = a + b, a + c
    def pmf(x):
        return comb(row1, x) * comb(n - row1, col1 - x) / comb(n, col1)
    p_obs = pmf(a)
    lo, hi = max(0, col1 - (n - row1)), min(row1, col1)
    return min(1.0, sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= p_obs * (1 + 1e-9)))


def main():
    pheno = pd.read_csv(PHENO, sep="\t", index_col=0)
    bias = pd.read_csv(BIAS, sep="\t").set_index("strain")
    strains = list(pheno.index)
    assert set(bias.index) == set(strains)

    B37 = sorted(bias.index[bias["class"] == "37C_biased"])
    B28 = sorted(bias.index[bias["class"] == "28C_biased"])
    MID = sorted(set(strains) - set(B37) - set(B28))
    print(f"B37 ({len(B37)}): {', '.join(B37)}")
    print(f"B28 ({len(B28)}): {', '.join(B28)}")
    print(f"neutral ({len(MID)})")

    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = hog[hog.iloc[:, 0].isin(strains)].pivot(
        index="HOG_ID", columns="SampleID", values="GeneCount")
    mat = mat.reindex(columns=strains).fillna(0)

    n_pres = (mat > 0).sum(axis=1)
    var = mat[strains].nunique(axis=1) > 1
    keep = mat[(n_pres >= CARRY_MIN) & (n_pres <= CARRY_MAX) & var]
    print(f"keep (variable, carried by {CARRY_MIN}-{CARRY_MAX} strains): {len(keep)} HOGs")

    both37 = keep[(keep[B37] > 0).all(axis=1)]
    sp37 = both37[(both37[B28 + MID] > 0).sum(axis=1) <= NEAR_37]
    all28 = keep[(keep[B28] > 0).all(axis=1)]
    sp28 = all28[(all28[B37 + MID] > 0).sum(axis=1) <= NEAR_28]
    print(f"sp37 (all {len(B37)} B37, <= {NEAR_37} carriers among B28+neutral): {len(sp37)} HOGs")
    print(f"sp28 (all {len(B28)} B28, <= {NEAR_28} carriers among B37+neutral): {len(sp28)} HOGs")

    for name, df in (("hog_sp37.tsv", sp37), ("hog_sp28.tsv", sp28)):
        path = os.path.join(OUTDIR, name)
        if os.path.exists(path) and not os.path.exists(path + ".bak-before-cmh"):
            os.rename(path, path + ".bak-before-cmh")
        df.reset_index().to_csv(path, sep="\t", index=False)
        print(f"wrote {path}")

    annot = pd.read_csv(ANNOT, sep="\t", dtype=str)
    cazy = dict(zip(annot.HOG_ID, annot.CAZy_Class.fillna("")))
    def is_cazy(h):
        v = cazy.get(h, "")
        return v not in ("", "nan", "-")

    rows = []
    for setname, df in (("sp37", sp37), ("sp28", sp28)):
        ids = list(df.index)
        k = sum(is_cazy(h) for h in ids)
        rows.append(dict(set=setname, n_hogs=len(ids), n_cazy=k,
                          frac_cazy=f"{k / len(ids):.3f}" if ids else "NA"))
        print(f"{setname}: {k}/{len(ids)} carry a CAZy class ({k / len(ids):.1%})")

    a, b = rows[0]["n_cazy"], rows[0]["n_hogs"] - rows[0]["n_cazy"]
    c, d = rows[1]["n_cazy"], rows[1]["n_hogs"] - rows[1]["n_cazy"]
    p = fisher_two_sided(a, b, c, d)
    rows[0]["fisher_P_cazy_vs_other"] = rows[1]["fisher_P_cazy_vs_other"] = f"{p:.4f}"
    print(f"CAZy fraction, sp37 vs sp28: Fisher two-sided P = {p:.4f}")

    summary = os.path.join(OUTDIR, "toolkit_hog_sets_summary.tsv")
    with open(summary, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["set", "n_hogs", "n_cazy", "frac_cazy",
                                           "fisher_P_cazy_vs_other"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
