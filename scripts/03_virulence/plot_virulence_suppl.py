#!/usr/bin/env python3
# DEPRECATED PHENOTYPE (2026-10-22): reads the mortality-era association table
# (nine phenotypes, 10^5 included). The live universe is 1 - AUC at 10^3 and 10^4
# with the two-dose means and their difference -> hog_phenotype_assoc_auc_no10e5.tsv.
# Kept only to regenerate the archived artefacts; do not use for new figures.
# -*- coding: utf-8 -*-
"""Supplementary figure variants — genome-wide association overview.
Three variants (user picks the best):
  V1: rho heatmap of variable HOGs x traits (clustered)
  V2: Manhattan-style -log10P across HOGs, 28C/37C panels
  V3: representative HOGs — per-strain copies vs mortality strips
"""

import os
import numpy as np
import pandas as pd
from scipy import stats, cluster
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
OUT = "07.lethal_groups/genomewide_assoc"
STRAINS = ["BMU01287","BMU01935","BMU01936","BMU01973","BMU02028","BMU02031",
           "BMU02407","BMU02479","BMU02745","BMU03976","BMU10147",
           "BMU12683","BMU12684","BMU12685","BMU12686"]
assoc = pd.read_csv(f"{OUT}/hog_phenotype_assoc.tsv", sep="\t")
traits = ["mort_28Ce3","mort_28Ce4","mort_28Ce5","mort_37Ce3","mort_37Ce4",
          "mort_37Ce5","mort_28_mean","mort_37_mean","temp_effect",
          "conc_slope_28","conc_slope_37"]
C28, C37 = "#3B7DD8", "#C44E52"

# rho matrix: HOG x trait
rho_mat = assoc.pivot(index="HOG_ID", columns="trait", values="rho").reindex(columns=traits)

# ---------------- V1: clustered rho heatmap ----------------
rho_mat = rho_mat.dropna()
# keep HOGs with at least one |rho| >= 0.5
sel = rho_mat[(rho_mat.abs() >= 0.5).any(axis=1)]
# cluster rows (ward on 1-correlation distance)
D = 1 - np.corrcoef(sel.values)
D = np.nan_to_num(D, nan=0.0)
D = (D + D.T) / 2
link = cluster.hierarchy.ward(D)
order = cluster.hierarchy.leaves_list(link)
sel = sel.iloc[order]

fig, ax = plt.subplots(figsize=(4.6, 6.0), dpi=300)
im = ax.imshow(sel.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(traits)))
ax.set_xticklabels([t.replace("mort_","").replace("conc_slope_","slp_")
                    .replace("temp_effect","T_eff") for t in traits],
                   rotation=90, fontsize=5.5)
ax.set_yticks([])
cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label("Spearman \u03c1 (copy vs trait)", fontsize=6)
ax.set_title(f"V1  Variable HOG x trait \u03c1 ({len(sel)} HOGs, |\u03c1|\u22650.5)",
             fontsize=7.5)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_gwa_variant1_heatmap.pdf")
print("V1 saved:", len(sel), "HOGs")

# ---------------- V2: Manhattan-style ----------------
g28 = assoc[assoc.trait.isin(["mort_28Ce3","mort_28Ce4","mort_28Ce5",
                              "mort_28_mean","conc_slope_28"])]
g37 = assoc[assoc.trait.isin(["mort_37Ce3","mort_37Ce4","mort_37Ce5",
                              "mort_37_mean","conc_slope_37"])]
p28 = g28.groupby("HOG_ID").p.min()
p37 = g37.groupby("HOG_ID").p.min()
hogs = sorted(set(p28.index) | set(p37.index))
x = np.arange(len(hogs))
p28 = p28.reindex(hogs).fillna(1)
p37 = p37.reindex(hogs).fillna(1)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 3.6), dpi=300,
                               sharex=True, gridspec_kw={"hspace": 0.08})
for ax, p, col, lab in [(ax1, p28, C28, "28 \u00b0C"), (ax2, p37, C37, "37 \u00b0C")]:
    yv = -np.log10(np.clip(p.values, 1e-8, 1))
    ax.scatter(x, yv, s=2.5, c=col, alpha=0.45, rasterized=True)
    ax.axhline(-np.log10(0.001), ls=":", color="grey", lw=0.8)
    ax.text(0.008, -np.log10(0.001) + 0.18, "P = 0.001",
            ha="left", fontsize=5.5, color="grey")
    ax.set_ylabel(f"{lab}\n\u2212log\u2081\u2080 P")
    ax.tick_params(labelbottom=False)
    ax.set_ylim(0, 5.2)
ax2.tick_params(labelbottom=True)
ax2.set_xlabel("Variable orthogroups (n = {:,})".format(len(hogs)))
ax1.set_title("V2  Manhattan-style genome-wide scan (min P per HOG)", fontsize=7.5)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_gwa_variant2_manhattan.pdf")
print("V2 saved:", len(hogs), "HOGs")

# ---------------- V3: representative HOG strips ----------------
rep = pd.read_excel("00.manuscript/supplementary/suppl_tables.virulence_assoc.xlsx",
                      sheet_name="S3_representative_HOGs")
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv", sep="\t")
mat = hog[hog.iloc[:,0].isin(STRAINS)].pivot(index="HOG_ID", columns="SampleID",
                                             values="GeneCount").reindex(columns=STRAINS)
pheno = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)
rows = rep[rep.Temperature_mode == "28C"].head(4)
rows = pd.concat([rows, rep[rep.Temperature_mode == "37C"].head(4)])

fig, axes = plt.subplots(1, len(rows), figsize=(9.0, 2.8), dpi=300)
for ax, (_, r) in zip(axes, rows.iterrows()):
    h, trait = r.Representative_HOG, r.Best_trait
    x = mat.loc[h]
    y = pheno[trait].reindex(STRAINS)
    rho, p = stats.spearmanr(x, y)
    for s in STRAINS:
        col = C28 if r.Temperature_mode == "28C" else C37
        ax.scatter(x[s], y[s], s=34, c=col, alpha=0.85,
                   edgecolors="white", linewidths=0.3, zorder=3)
    ax.set_title(f"{h.split('N0.')[-1]}\n{r.Pfam_domain} ({r.Temperature_mode})\n"
                 f"\u03c1={rho:.2f}", fontsize=5.2)
    ax.set_xlabel("copies", fontsize=5.5)
    ax.tick_params(labelsize=5)
    ax.set_ylim(-0.05, 1.1)
fig.suptitle("V3  Representative temperature-associated HOGs: copy number vs trait",
             fontsize=8, y=1.02)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_gwa_variant3_rep_hogs.pdf")
print("V3 saved:", len(rows), "HOGs")
