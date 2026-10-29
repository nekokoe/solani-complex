#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S26 — literature pathogenicity genes summary (nature-figure, python).

a: copy-number heatmap with numeric values
b: best |rho| of each gene across conditions (ranked bar, condition labelled)
c: Spearman-rho heatmap (gene x condition) with numeric values
d: best-hit identity heatmap with numeric values
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
OUT = "07.lethal_groups/literature_genes"
C28, C37 = "#3B7DD8", "#C44E52"
# Same seven phenotypes as the orthogroup screen: 1 - AUC at 10^3 and 10^4, the two
# two-dose means and their difference. The 10^5 dose is excluded.
CONDS = ["vir_28Ce3", "vir_28Ce4", "vir_37Ce3", "vir_37Ce4",
         "vir_28_mean2", "vir_37_mean2", "vir_temp_effect2"]
LABELS = {"vir_28Ce3": "28\u00b0C 1e3", "vir_28Ce4": "28\u00b0C 1e4",
          "vir_37Ce3": "37\u00b0C 1e3", "vir_37Ce4": "37\u00b0C 1e4",
          "vir_28_mean2": "28\u00b0C mean", "vir_37_mean2": "37\u00b0C mean",
          "vir_temp_effect2": "\u0394 temp\u00b0C"}

copy = pd.read_csv(f"{OUT}/gene_strain_copy.tsv", sep="\t", index_col=0)
ident = pd.read_csv(f"{OUT}/gene_strain_identity.tsv", sep="\t", index_col=0)
assoc = pd.read_csv(f"{OUT}/gene_phenotype_assoc.tsv", sep="\t")
genes = copy.columns.tolist()
STRAINS = copy.index.tolist()

fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0), dpi=300)
(axa, axb), (axc, axd) = axes

# ---------- a: copy-number heatmap with values ----------
copy_g = copy.T.reindex(genes)
im = axa.imshow(copy_g.values, cmap="YlOrRd", aspect="auto")
axa.set_xticks(range(len(STRAINS)))
axa.set_xticklabels([s[3:] for s in STRAINS], rotation=90, fontsize=5)
axa.set_yticks(range(len(genes)))
axa.set_yticklabels(genes, fontsize=5.5)
for i in range(len(genes)):
    for j in range(len(STRAINS)):
        v = copy_g.iloc[i, j]
        axa.text(j, i, str(int(v)), ha="center", va="center", fontsize=3.6,
                 color="black")
axa.set_title("a  Copy number", fontsize=8, loc="left")
fig.colorbar(im, ax=axa, fraction=0.04, pad=0.02, label="Copies")

# ---------- b: best |rho| per gene (ranked) ----------
best = assoc.loc[assoc.groupby("gene").p.idxmin()]
best = best.sort_values("spearman_rho", key=lambda x: np.abs(x))
best = best.iloc[::-1]
ys = np.arange(len(best))
cols = [C28 if str(c).startswith("28") else C37 for c in best.cond]
axb.barh(ys, np.abs(best.spearman_rho), color=cols, alpha=0.8, height=0.7)
for y, (_, r) in zip(ys, best.iterrows()):
    axb.text(np.abs(r.spearman_rho) + 0.01, y,
             f"{np.abs(r.spearman_rho):.2f}  ({LABELS.get(r.cond, r.cond)})",
             va="center", fontsize=4.6)
axb.set_yticks(ys)
axb.set_yticklabels(best.gene, fontsize=5.5)
axb.set_xlim(0, 0.75)
axb.set_xlabel("Best |Spearman \u03c1|", fontsize=7)
axb.set_title("b  Best association of each gene (across conditions)",
              fontsize=8, loc="left")
axb.tick_params(labelsize=5.5)

# ---------- c: rho heatmap with values ----------
piv = assoc.pivot(index="gene", columns="cond", values="spearman_rho").reindex(genes)
piv = piv[CONDS]
# FDR-significant gene x condition pairs, computed from the table.
# An earlier version hard-coded {"velvet", "VeA", "LaeA"} here together with the
# claim "* FDR < 0.05" in the panel title. No product ever supported it: the
# minimum BH q is 0.48 in the 21-strain cohort and 0.51 in the 15-strain one, so
# the asterisks were asserting a significance that held in neither. Deriving the
# set here means the asterisks appear exactly when the data earn them.
SIG = {(r.gene, r.cond) for r in assoc.itertuples() if r.q < 0.05}
im = axc.imshow(piv.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
axc.set_xticks(range(len(CONDS)))
axc.set_xticklabels([LABELS[c] for c in CONDS], rotation=45, ha="right", fontsize=5)
axc.set_yticks(range(len(genes)))
axc.set_yticklabels(genes, fontsize=5.5)
for i in range(len(genes)):
    for j in range(len(CONDS)):
        v = piv.iloc[i, j]
        if not np.isnan(v):
            star = "*" if (genes[i], CONDS[j]) in SIG else ""
            axc.text(j, i, f"{v:.1f}{star}", ha="center", va="center", fontsize=3.8,
                     color="white" if abs(v) > 0.5 else "black")
axc.set_title("c  Spearman \u03c1 (copy vs mortality); "
              + ("* BH q < 0.05" if SIG else "no pair reaches BH q < 0.05"),
              fontsize=8, loc="left")
fig.colorbar(im, ax=axc, fraction=0.04, pad=0.02)

# ---------- d: identity heatmap with values ----------
ident_v = (ident * 100).T.reindex(genes)
im = axd.imshow(ident_v.values, cmap="viridis", vmin=0, vmax=100, aspect="auto")
axd.set_xticks(range(len(STRAINS)))
axd.set_xticklabels([s[3:] for s in STRAINS], rotation=90, fontsize=5)
axd.set_yticks(range(len(genes)))
axd.set_yticklabels(genes, fontsize=5.5)
for i in range(len(genes)):
    for j in range(len(STRAINS)):
        v = ident_v.iloc[i, j]
        axd.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=3.6,
                 color="white" if v > 55 else "black")
axd.set_title("d  Best-hit identity (%)", fontsize=8, loc="left")
fig.colorbar(im, ax=axd, fraction=0.04, pad=0.02)

fig.subplots_adjust(left=0.07, right=0.94, bottom=0.07, top=0.94,
                    wspace=0.32, hspace=0.30)
fig.savefig(f"{OUT}/figure_literature_S26.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_literature_S26.png", dpi=300, bbox_inches="tight")
print("S26 saved")
