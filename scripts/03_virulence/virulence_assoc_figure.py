#!/usr/bin/env python3
# DEPRECATED PHENOTYPE (2026-10-22): reads the mortality-era association table
# (nine phenotypes, 10^5 included). The live universe is 1 - AUC at 10^3 and 10^4
# with the two-dose means and their difference -> hog_phenotype_assoc_auc_no10e5.tsv.
# Kept only to regenerate the archived artefacts; do not use for new figures.
# -*- coding: utf-8 -*-
"""Step 4: integrative figure for literature genes + genome-wide hits."""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
OUT = "07.lethal_groups"
STRAINS = ["BMU01287","BMU01935","BMU01936","BMU01973","BMU02028","BMU02031",
           "BMU02407","BMU02479","BMU02745","BMU03976","BMU10147",
           "BMU12683","BMU12684","BMU12685","BMU12686"]

pheno = pd.read_csv(f"{OUT}/phenotype_analysis/phenotype_metrics.tsv", sep="\t", index_col=0)
copy = pd.read_csv(f"{OUT}/literature_genes/gene_strain_copy.tsv", sep="\t", index_col=0)
assoc = pd.read_csv(f"{OUT}/genomewide_assoc/hog_phenotype_assoc.tsv", sep="\t")
annot = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv", sep="\t")

fig, axes = plt.subplots(2, 2, figsize=(11, 7), dpi=300)
C28, C37 = "#C44E52", "#3B7DD8"

# a: literature gene copy vs mortality at 28Ce3 (velvet/VeA/LaeA)
ax = axes[0, 0]
for g, c, mk in [("velvet", "#984EA3", "o"), ("VeA", "#377EB8", "s"), ("LaeA", "#E69F00", "^")]:
    x = copy[g]; y = pheno["mort_28Ce3"].reindex(copy.index)
    rho, p = stats.spearmanr(x, y)
    ax.scatter(x, y, s=34, c=c, marker=mk, label=f"{g} (rho={rho:.2f}, p={p:.3g})")
ax.set_xlabel("Copies per strain")
ax.set_ylabel("Mortality 28C e3")
ax.set_title("a  Velvet-complex copies vs 28C e3 mortality")
ax.legend(fontsize=6.5, frameon=False)

# b: genome-wide top hit vs 28Ce4 (dose gradient)
ax = axes[0, 1]
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv", sep="\t")
mat = hog[hog.iloc[:,0].isin(STRAINS)].pivot(index="HOG_ID", columns="SampleID",
                                             values="GeneCount").reindex(columns=STRAINS)
x = mat.loc["N0.HOG0000412"]; y = pheno["mort_28Ce4"].reindex(STRAINS)
rho, p = stats.spearmanr(x, y)
ax.scatter(x, y, s=40, c=C28)
for s, dx, dy in [("BMU03976", -0.5, 0.05), ("BMU01287", 0.06, -0.06)]:
    ax.annotate(s.replace("BMU",""), (x[s]+dx, y[s]+dy), fontsize=6, fontweight="bold")
ax.set_xlabel("HOG0000412 copies (MFS transporter / CFEM / F-box)")
ax.set_ylabel("Mortality 28C e4")
ax.set_title(f"b  HOG0000412 vs 28C e4 (rho={rho:.2f}, p={p:.3g})")

# c: Manhattan-style top association summary per trait
ax = axes[1, 0]
traits = ["mort_28Ce3","mort_28Ce4","mort_28Ce5","mort_37Ce3","mort_37Ce4",
          "mort_37Ce5","temp_effect","conc_slope_28","conc_slope_37"]
colors = [C28]*3 + [C37]*3 + ["#555555","#E69F00","#4DAF4A"]
for i, t in enumerate(traits):
    sub = assoc[assoc.trait == t].sort_values("p").head(8)
    ax.scatter(np.full(len(sub), i) + (np.arange(len(sub)) - 4) * 0.08,
               -np.log10(sub.p), s=20, c=colors[i], alpha=0.7)
ax.set_xticks(range(len(traits)))
ax.set_xticklabels([t.replace("mort_","").replace("conc_slope_","slope_")
                    .replace("temp_effect","T_eff") for t in traits],
                   rotation=90, ha="center", fontsize=5.5)
ax.set_ylabel("-log10 p (top 8 HOGs per trait)")
ax.set_title("c  Top genome-wide associations per trait")
ax.axhline(-np.log10(0.001), ls=":", color="grey", lw=0.8)

# d: condition specificity of the temp-effect HOG
ax = axes[1, 1]
x = mat.loc["N0.HOG0017677"]
m28 = pheno[["mort_28Ce3","mort_28Ce4","mort_28Ce5"]].mean(axis=1)
m37 = pheno[["mort_37Ce3","mort_37Ce4","mort_37Ce5"]].mean(axis=1)
carry = x > 0
ax.scatter(m28[~carry], m37[~carry], s=40, c="#999999", label="no copy")
ax.scatter(m28[carry], m37[carry], s=40, c=C28, label="carry HOG0017677")
for s, dx, dy in [("BMU03976", 0.01, 0.03), ("BMU01287", 0.02, -0.05),
                    ("BMU02031", 0.01, 0.03), ("BMU12685", -0.04, 0.02)]:
    ax.annotate(s.replace("BMU",""), (m28[s]+dx, m37[s]+dy), fontsize=6)
ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=0.8)
ax.set_xlabel("Mean mortality 28C"); ax.set_ylabel("Mean mortality 37C")
ax.set_title("d  HOG0017677 carriers shift toward 28C virulence")
ax.legend(fontsize=6.5, frameon=False)

for axx in axes.flat:
    for s in ["top", "right"]:
        axx.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_virulence_assoc.pdf")
fig.savefig(f"{OUT}/figure_virulence_assoc.png")
print("saved ->", f"{OUT}/figure_virulence_assoc.pdf")
