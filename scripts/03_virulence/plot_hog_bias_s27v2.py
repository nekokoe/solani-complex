#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S27 v2 (nature-figure, python) — temperature-bias strain-specific HOGs.

One integrated figure contrasting the 37C-biased and 28C-biased strain
toolkits:

- Top strip: strain grouping (37C-biased red / intermediate grey / 28C-biased blue)
- Heatmap: 37C-specific HOGs (warm YlOrRd) on top, 28C-specific HOGs
  (cool PuBu) below, separated by a black divider
- Row labels: two-line — HOG ID (grey small) + gene/function (black);
  chitin-related HOGs in red bold
- Left edge: per-row group tick (red/blue)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

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
B28 = ["BMU03976", "BMU12685", "BMU02028", "BMU01936"]
MID = ["BMU12684", "BMU02407", "BMU02479", "BMU01973", "BMU01935",
       "BMU12686", "BMU12683", "BMU10147", "BMU02745"]
STRAINS = B37 + MID + B28

CHIT = ["Glyco_hydro_18", "Glyco_hydro_19", "Glyco_hydro_20", "Glyco_hydro_12",
        "Chitin_synth", "chitin", "Glyco_hydro_28"]
INF = ["Hce2", "CFEM", "LysM", "HET", "Peptidase_S", "Effector"]

pheno = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv",
                  sep="\t")
mat = hog[hog.iloc[:, 0].isin(pheno.index)].pivot(
    index="HOG_ID", columns="SampleID", values="GeneCount")
mat = mat.reindex(columns=pheno.index).fillna(0)
n_pres = (mat > 0).sum(axis=1)
var = mat[pheno.index].nunique(axis=1) > 1
keep = mat[(n_pres >= 3) & (n_pres <= 13) & var]

both37 = keep[(keep[B37] > 0).all(axis=1)]
sp37 = both37[(both37[B28 + MID] > 0).sum(axis=1) <= 3]
all28 = keep[(keep[B28] > 0).all(axis=1)]
sp28 = all28[(all28[B37 + MID] > 0).sum(axis=1) <= 4]

annot = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv",
                    sep="\t")
def geth(h, col):
    r = annot[annot.HOG_ID == h]
    if len(r) == 0:
        return ""
    s = str(r[col].iloc[0])
    return s if s.lower() not in ("nan", "-") else ""

def build(hset):
    rows = []
    for h in hset.index:
        pf = geth(h, "Pfam_Domain")
        gn = geth(h, "Gene_Name")
        pfs = pf.split(";") if pf else []
        chit = any(c in pf for c in CHIT) or "chitin" in gn.lower()
        inf = any(c in pf for c in INF)
        tag = (gn.split(";")[0] if gn else "") or (pfs[0] if pfs else "")
        rows.append({"HOG": h, "tag": tag, "chit": chit, "inf": inf,
                     "n_other": int((hset.drop(columns=[c for c in hset.columns
                                                        if c in B37]) > 0).sum(axis=1)[h])})
    return pd.DataFrame(rows).set_index("HOG")

m37 = build(sp37).sort_values(["chit", "inf", "n_other"], ascending=[False, False, True])
m28 = build(sp28).sort_values(["chit", "inf", "n_other"], ascending=[False, False, True])

def rows_block(m, block, H, group_color):
    """draw one heatmap block + row group ticks + two-line labels"""
    M = mat.loc[m.index][STRAINS]
    im = block.imshow(np.clip(M.values, 0, 5), cmap=H, vmin=0, vmax=5,
                      aspect="auto", interpolation="nearest")
    block.set_xticks([])
    block.set_yticks([])
    for i, (h, r) in enumerate(m.iterrows()):
        # group tick on the left
        block.add_patch(Rectangle((-0.55, i - 0.5), 0.28, 1, color=group_color,
                                  clip_on=False))
        # two-line label: HOG id (grey, top) + gene/function (black, bottom)
        col = "#C44E52" if r.chit else ("#E08E0B" if r.inf else "#222222")
        weight = "bold" if r.chit else "normal"
        block.text(-0.62, i + 0.26, h.replace("N0.HOG", ""), fontsize=4.2,
                   color="#888888", ha="right", va="center", clip_on=False)
        block.text(-0.62, i - 0.26, (r.tag[:20] if r.tag else "\u2014"),
                   fontsize=5.2, color=col, ha="right", va="center",
                   weight=weight, clip_on=False)
        # cell values
        for j in range(len(STRAINS)):
            v = int(M.iloc[i, j])
            if v > 0:
                block.text(j, i, str(v), ha="center", va="center", fontsize=3.2,
                           color="white" if v >= 3 else "#333333")
    return im

fig = plt.figure(figsize=(10.8, 13.5), dpi=300)
gs = fig.add_gridspec(2, 1, height_ratios=[0.045, 1], hspace=0.0)
axstrip = fig.add_subplot(gs[0, 0])
axh = fig.add_subplot(gs[1, 0])

n37, n28 = len(m37), len(m28)
N = n37 + n28
# strain group strip: dedicated top axes
for j, s in enumerate(STRAINS):
    c = "#C44E52" if s in B37 else ("#3B7DD8" if s in B28 else "#CCCCCC")
    axstrip.add_patch(Rectangle((j - 0.5, 0), 1, 1, color=c))
axstrip.set_xlim(-0.5, len(STRAINS) - 0.5)
axstrip.set_ylim(0, 1)
axstrip.axis("off")

# heatmap: 37C block (warm) on top, 28C block (cool) below
im37 = rows_block(m37, axh, "YlOrRd", "#C44E52")
im28 = rows_block(m28, axh, "PuBu", "#3B7DD8")
# divider between blocks
axh.axhline(n37 - 0.5, color="black", lw=1.6)
axh.set_xlim(-0.62, len(STRAINS) - 0.5)
axh.set_ylim(-0.5, N - 0.5)
axh.invert_yaxis()

# strain labels below
axh.set_xticks(range(len(STRAINS)))
axh.set_xticklabels([s[3:] for s in STRAINS], rotation=90, fontsize=5.5)
for j, s in enumerate(STRAINS):
    if s in B37:
        axh.get_xticklabels()[j].set_color("#C44E52")
        axh.get_xticklabels()[j].set_fontweight("bold")
    elif s in B28:
        axh.get_xticklabels()[j].set_color("#3B7DD8")
        axh.get_xticklabels()[j].set_fontweight("bold")
# legend for the strip
axstrip.set_title("37 \u00b0C-biased (red)  /  intermediate (grey)  /  "
                  "28 \u00b0C-biased (blue)", fontsize=6, color="#555555",
                  loc="left", pad=1)

cb = fig.colorbar(im37, ax=axh, fraction=0.025, pad=0.015)
cb.set_label("Copy number", fontsize=6)
cb.ax.tick_params(labelsize=5)

# legends
handles = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#C44E52",
           markersize=7, label="37 \u00b0C-biased strain"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#3B7DD8",
           markersize=7, label="28 \u00b0C-biased strain"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#C44E52",
           markersize=7, label="Chitin-related HOG"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor="#E08E0B",
           markersize=7, label="Infection-relevant HOG"),
]
fig.legend(handles=handles, fontsize=6, loc="lower center",
           bbox_to_anchor=(0.55, 0.005), ncol=4, columnspacing=1.0)

fig.suptitle("Strain-specific orthogroups of temperature-biased isolates",
             fontsize=9, y=0.99)
fig.text(0.13, 0.955, "Upper block: 37 \u00b0C-biased toolkit (n = %d);  "
                      "lower block: 28 \u00b0C-biased toolkit (n = %d)"
         % (n37, n28), fontsize=6.5, color="#333333")

fig.subplots_adjust(left=0.13, right=0.96, bottom=0.085, top=0.97)
fig.savefig(f"{OUT}/figure_hog_bias_S27v2.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_hog_bias_S27v2.png", dpi=300, bbox_inches="tight")
print("S27v2 saved: 37C", n37, "+ 28C", n28, "= ", N)
print("chitin 37C:", int(m37.chit.sum()), "| chitin 28C:", int(m28.chit.sum()))
