#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S27 — orthogroups enriched in 37C-biased strains (BMU01287/BMU02031).

46 variable HOGs carried by BOTH 37C-biased strains but rare (<=3/13) in the
other strains. Heatmap of copy numbers across the 15 strains; rows annotated
with a short functional tag; chitin-related and infection-relevant HOGs
highlighted.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 7,
    "legend.frameon": False,
})

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
OUT = "07.lethal_groups/genomewide_assoc"
B37 = ["BMU01287", "BMU02031"]

pheno = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)
other = [s for s in pheno.index if s not in B37]
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv",
                  sep="\t")
mat = hog[hog.iloc[:, 0].isin(pheno.index)].pivot(
    index="HOG_ID", columns="SampleID", values="GeneCount")
mat = mat.reindex(columns=pheno.index).fillna(0)
n_pres = (mat > 0).sum(axis=1)
var = mat[pheno.index].nunique(axis=1) > 1
keep = mat[(n_pres >= 3) & (n_pres <= 13) & var]

both = keep[(keep[B37] > 0).all(axis=1)]
rare = both[(both[other] > 0).sum(axis=1) <= 3]
rare = rare.sort_index()

annot = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv",
                    sep="\t")
def geth(h, col):
    r = annot[annot.HOG_ID == h]
    if len(r) == 0:
        return ""
    s = str(r[col].iloc[0])
    return s if s.lower() not in ("nan", "-") else ""

CHIT = ["Glyco_hydro_18", "Glyco_hydro_19", "Glyco_hydro_20", "Glyco_hydro_12",
        "Chitin_synth", "chitin"]
INF = ["Hce2", "CFEM", "LysM", "HET", "Peptidase_S", "Effector",
       "Pkinase", "Carboxypep"]
rows = []
for h in rare.index:
    pf = geth(h, "Pfam_Domain")
    gn = geth(h, "Gene_Name")
    pfs = pf.split(";") if pf else []
    is_chit = any(c in pf for c in CHIT) or any(c in gn for c in ["chitin", "Chi"])
    is_inf = any(c in pf for c in INF)
    tag = (gn.split(";")[0] if gn else "") or (pfs[0] if pfs else h)
    rows.append({"HOG": h, "tag": tag, "chit": is_chit, "inf": is_inf})
meta = pd.DataFrame(rows).set_index("HOG")
# sort: chitin first, then infection, then by rarity
meta["rank"] = np.where(meta.chit, 0, np.where(meta.inf, 1, 2))
meta = meta.sort_values(["rank", "tag"])

M = rare.loc[meta.index].copy()
M = M[B37 + other]

fig, (axg, axs) = plt.subplots(
    2, 1, figsize=(7.2, 9.5), dpi=300,
    gridspec_kw={"height_ratios": [8, 1], "hspace": 0.02})
C0 = "#F2F2F2"
Cmax = "#7F2704"
im = axg.imshow(np.clip(M.values, 0, 5), cmap="YlOrBr", vmin=0, vmax=5,
                aspect="auto")
axg.set_yticks(range(len(M)))
axg.set_yticklabels([f"{r.tag}" for _, r in meta.iterrows()], fontsize=5)
for i, (h, r) in enumerate(meta.iterrows()):
    col = "#C44E52" if r.chit else ("#E08E0B" if r.inf else "#555555")
    axg.get_yticklabels()[i].set_color(col)
axg.set_xticks(range(len(M.columns)))
axg.set_xticklabels([s[3:] for s in M.columns], rotation=90, fontsize=5)
for j, s in enumerate(M.columns):
    if s in B37:
        axg.get_xticklabels()[j].set_color("#C44E52")
        axg.get_xticklabels()[j].set_fontweight("bold")
axg.set_title(
    "S27  Orthogroups enriched in 37 \u00b0C-biased strains\n"
    "(carried by both BMU01287 and BMU02031, rare elsewhere)",
    fontsize=8, loc="left")
# cell values
for i in range(len(M)):
    for j in range(len(M.columns)):
        v = int(M.iloc[i, j])
        if v > 0:
            axg.text(j, i, str(v), ha="center", va="center", fontsize=3.4,
                     color="white" if v >= 3 else "#333333")
cb = fig.colorbar(im, ax=axg, fraction=0.035, pad=0.02)
cb.set_label("Copy number", fontsize=6)
cb.ax.tick_params(labelsize=5)

# legend strip
axs.axis("off")
handles = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#C44E52",
           markersize=7, label="Chitin-related"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#E08E0B",
           markersize=7, label="Infection-relevant (HET/MFS/glycoside\nhydrolase/peptidase)"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#555555",
           markersize=7, label="Other"),
]
axs.legend(handles=handles, fontsize=6, loc="center", ncol=1, frameon=False)
axs.text(0.5, 0.05, "Strain labels in red: 37 \u00b0C-biased isolates",
         transform=axs.transAxes, ha="center", fontsize=6, color="#C44E52")

fig.savefig(f"{OUT}/figure_hog_enrich_b37_S27.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_hog_enrich_b37_S27.png", dpi=300, bbox_inches="tight")
print("S27 saved:", len(M), "HOGs; chitin:", int(meta.chit.sum()),
      "; infection:", int(meta.inf.sum()))
