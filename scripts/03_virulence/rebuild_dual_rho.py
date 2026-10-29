#!/usr/bin/env python3
# DEPRECATED PHENOTYPE (2026-10-22): reads the mortality-era association table
# (nine phenotypes, 10^5 included). The live universe is 1 - AUC at 10^3 and 10^4
# with the two-dose means and their difference -> hog_phenotype_assoc_auc_no10e5.tsv.
# Kept only to regenerate the archived artefacts; do not use for new figures.
"""Rebuild hog_dual_rho.tsv on the 21-strain, nine-phenotype analysis.

The previous file was produced ad hoc during the 15-strain era and never had a
generator, so it still carried 3,630 orthogroups. figure 7 panel e reads it for
the background cloud and for the HOG count in its title, which therefore said
"3,630 variable HOGs" while every other analysis in the paper used 4,484.

Definitions (unchanged from the 15-strain file, applied to the new universe):
  rho28   max |Spearman rho| of that HOG against mortality at 28 C
          (10^3, 10^4, 10^5 and the 28 C mean)
  rho37   the same for 37 C
  rows    the variable HOGs of the 21-strain cohort (present in 3-19 strains)

Input : 07.lethal_groups/genomewide_assoc/hog_phenotype_assoc.tsv
        05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv
Output: 07.lethal_groups/genomewide_assoc/hog_dual_rho.tsv
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ASSOC = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc/hog_phenotype_assoc.tsv")
PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
OUT = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc/hog_dual_rho.tsv")

STRAINS = ["BMU" + x.zfill(5) for x in
           ["9608", "7029", "7936", "1336", "2406", "5572", "2479", "2031", "12685",
            "12686", "1287", "1935", "1936", "12683", "12684", "2407", "3976", "2028",
            "1973", "2745", "10147"]]
T28 = ["mort_28Ce3", "mort_28Ce4", "mort_28Ce5", "mort_28_mean"]
T37 = ["mort_37Ce3", "mort_37Ce4", "mort_37Ce5", "mort_37_mean"]

assoc = pd.read_csv(ASSOC, sep="\t")
pav = pd.read_csv(PAV, sep="\t", dtype={"Present": "int8"})
wide = (pav[pav.SampleID.isin(STRAINS)]
        .pivot(index="SampleID", columns="HOG_ID", values="Present").fillna(0).loc[STRAINS])
prev = wide.sum(0)
universe = sorted(wide.columns[(prev >= 3) & (prev <= 19)])


def best(traits):
    sub = assoc[assoc.trait.isin(traits) & assoc.HOG_ID.isin(universe)]
    return sub.groupby("HOG_ID")["rho"].apply(lambda x: x.abs().max())


dual = pd.DataFrame({"HOG": universe})
dual["rho28"] = dual["HOG"].map(best(T28))
dual["rho37"] = dual["HOG"].map(best(T37))
dual.to_csv(OUT, sep="\t", index=False)

print(f"wrote {OUT}")
print(f"  rows {len(dual):,} (was 3,630) · rho28/rho37 无缺失: "
      f"{int(dual.rho28.notna().sum())}/{int(dual.rho37.notna().sum())}")
print(f"  rho28 >= 0.50: {int((dual.rho28 >= 0.50).sum())} · "
      f"rho37 >= 0.50: {int((dual.rho37 >= 0.50).sum())}")

# 核验正文里引用的具体 rho
CHECK = {"N0.HOG0017677": "amino-acid permease (text: rho = 0.86)",
         "N0.HOG0017872": "VPS1/dynamin",
         "N0.HOG0002312": "ribosomal S10"}
print("\n正文引用的 HOG 在新表里的读数:")
for h, note in CHECK.items():
    if h in set(dual.HOG):
        r = dual[dual.HOG == h].iloc[0]
        sub = assoc[(assoc.HOG_ID == h) & (assoc.rho.abs() >= 0.5)].sort_values("p")
        best_row = sub.iloc[0] if len(sub) else None
        print(f"  {h} ({note}): rho28={r.rho28:.3f} rho37={r.rho37:.3f}"
              + (f" · 最强 {best_row.trait} rho={best_row.rho:+.3f} p={best_row.p:.2g}"
                 if best_row is not None else " · 无 |rho|>=0.5 的表型"))
    else:
        print(f"  {h} ({note}): 不在新宇宙中")
