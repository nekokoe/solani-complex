#!/usr/bin/env python3
"""Supplementary Fig. S28 — cytokine induction, both macrophage models.

Rebuilt on the 2026-10-21 source data (`parse_cell_source_2026_10_21.py`). Two things
changed relative to the earlier version, which read the 2026-10-12 workbook:

  * there are now TWO models (THP-1 and PBMC-derived), so every panel carries a model
    axis instead of assuming THP-1;
  * ELISA exists for IL-1beta and TNF-alpha only (the old workbook had TNF-alpha and
    CSF2), so the readout list is four mRNAs + two proteins.

The low-concentration block has no Ctrl/LPS of its own; the shared control pair is
carried over from the 10^5 block by the parser. Because that carry-over is an
assumption, the source comparison (panel e) is shown at 10^5 only, where the control
was measured.

Panels
  a THP-1, 10^5      b THP-1, 10^3      c PBMC, 10^5      d PBMC, 10^3   (fold vs control)
  e induction by isolate source, Mann-Whitney, 10^5, both models
  f Spearman rho of induction (10^5) against the Galleria phenotypes, both models

Outputs  07.lethal_groups/cell_assay/figure_cell_S28.{pdf,png}
         07.lethal_groups/cell_assay/cell_clinical_vs_env.tsv
         07.lethal_groups/cell_assay/cell_vs_phenotype.tsv
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
OUT = "07.lethal_groups/cell_assay"
PHENO = "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv"
META = "05.pan-genome/bmu_ref_groups.211.txt"
C_ENV = "#BBBBBB"; C_CLIN = "#C44E52"; C_THP = "#4D8AC1"; C_PBMC = "#C92626"
ORDER = ["IL1b", "TNFa", "cxcl1", "csf2"]
LBL = {"IL1b": "IL-1\u03b2", "TNFa": "TNF-\u03b1", "cxcl1": "CXCL1", "csf2": "CSF2"}
MODEL_LBL = {"THP1": "THP-1-derived", "PBMC": "PBMC-derived"}

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "pdf.fonttype": 42, "svg.fonttype": "none"})

s = pd.read_csv(f"{OUT}/cell_assay_summary.tsv", sep="\t")
meta = pd.read_csv(META, sep="\t").set_index("BMU_ID")
pheno = pd.read_csv(PHENO, sep="\t", index_col=0)
s["group"] = np.where(meta.MedEnv.reindex(s.strain).values == "med",
                      "clinical", "environmental")

READOUTS = [(c, "mRNA") for c in ORDER] + [("IL1b", "protein"), ("TNFa", "protein")]


def label(cy, asy):
    return LBL[cy] + (" protein" if asy == "protein" else " mRNA")


def matrix(model, conc):
    cols = {}
    for cy, asy in READOUTS:
        sub = s[(s.model == model) & (s.cytokine == cy) & (s.assay == asy)
                & (s.concentration == conc)]
        if len(sub):
            cols[label(cy, asy)] = sub.set_index("strain").fold_vs_ctrl
    m = pd.DataFrame(cols).T
    return m[sorted(m.columns)]


fig = plt.figure(figsize=(13.4, 10.6), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.15, 1.0],
                      hspace=0.42, wspace=0.62)

# ---- a-d: induction heatmaps, model x concentration ----
panels = [("a", "THP1", 1e5), ("b", "THP1", 1e3), ("c", "PBMC", 1e5), ("d", "PBMC", 1e3)]
for k, (letter, model, conc) in enumerate(panels):
    ax = fig.add_subplot(gs[k // 2, k % 2])
    M = matrix(model, conc)
    key = label("TNFa", "mRNA")
    order = M.loc[key].sort_values(ascending=False).index if key in M.index else M.columns
    M = M[order]
    norm = TwoSlopeNorm(vmin=0.5, vcenter=1.0, vmax=max(3.0, np.nanmax(M.values)))
    im = ax.imshow(M.values, aspect="auto", cmap="RdYlBu_r", norm=norm)
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels([f"{t}  (mean {M.loc[t].mean():.1f}\u00d7)" for t in M.index],
                       fontsize=6.5)
    ax.set_xticks(range(len(M.columns)))
    ax.set_xticklabels([t.replace("BMU", "") for t in M.columns], rotation=90,
                       fontsize=4.6)
    for tick, st in zip(ax.get_xticklabels(), M.columns):
        if meta.MedEnv.get(st) != "med":
            tick.set_color("#2E8B57")
    ax.set_title(f"{letter}  {MODEL_LBL[model]}, "
                 f"{'10$^5$' if conc == 1e5 else '10$^3$'} conidia", fontsize=8, loc="left")
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.035)
    cb.set_label("fold vs control", fontsize=7, labelpad=6)
    cb.ax.tick_params(labelsize=6)

# ---- e: clinical vs environmental, 10^5 ----
axe = fig.add_subplot(gs[2, 0])
rows = []
for model in ("THP1", "PBMC"):
    for cy, asy in READOUTS:
        sub = s[(s.model == model) & (s.cytokine == cy) & (s.assay == asy)
                & (s.concentration == 1e5)]
        a = sub[sub.group == "clinical"].fold_vs_ctrl.dropna()
        b = sub[sub.group == "environmental"].fold_vs_ctrl.dropna()
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"model": model, "readout": label(cy, asy), "conc": 1e5,
                     "n_clin": len(a), "n_env": len(b),
                     "median_clin": a.median(), "median_env": b.median(), "U": u, "p": p})
cw = pd.DataFrame(rows)
cw.to_csv(f"{OUT}/cell_clinical_vs_env.tsv", sep="\t", index=False)
x = np.arange(len(cw))
axe.bar(x - 0.2, cw.median_clin, 0.38, color=C_CLIN, label="clinical (17)")
axe.bar(x + 0.2, cw.median_env, 0.38, color=C_ENV, label="environmental (4)")
axe.set_yscale("log")
axe.axhline(1, ls=":", color="grey", lw=0.8)
for i, r in enumerate(cw.itertuples()):
    axe.text(i, max(r.median_clin, r.median_env) * 1.15, f"{r.p:.2f}",
             ha="center", fontsize=5.5)
axe.set_xticks(x)
short = {"THP1": "THP-1", "PBMC": "PBMC"}
axe.set_xticklabels([f"{r.readout}\n{short[r.model]}" for r in cw.itertuples()],
                    rotation=90, fontsize=5)
axe.set_ylabel("median fold change vs control", fontsize=7)
axe.set_title("e  Induction by isolate source (Mann-Whitney P, 10$^5$)", fontsize=8, loc="left")
axe.legend(fontsize=6, frameon=False)
for sp in ("top", "right"):
    axe.spines[sp].set_visible(False)

# ---- f: correlation with the Galleria phenotypes, 10^5 ----
axf = fig.add_subplot(gs[2, 1])
cors = []
for model in ("THP1", "PBMC"):
    for cy, asy in READOUTS:
        sub = (s[(s.model == model) & (s.cytokine == cy) & (s.assay == asy)
                 & (s.concentration == 1e5)].set_index("strain").fold_vs_ctrl)
        for trait, tl in (("mort_28_mean", "28 \u00b0C"), ("mort_37_mean", "37 \u00b0C"),
                          ("temp_effect", "temp. pref.")):
            y = pheno[trait].reindex(sub.index)
            ok = sub.notna() & y.notna()
            if ok.sum() < 5:
                continue
            rho, p = stats.spearmanr(sub[ok], y[ok])
            cors.append({"model": model, "readout": label(cy, asy), "trait": tl,
                         "conc": 1e5, "rho": rho, "p": p, "n": int(ok.sum())})
cd = pd.DataFrame(cors)
cd.to_csv(f"{OUT}/cell_vs_phenotype.tsv", sep="\t", index=False)
traits = ["28 \u00b0C", "37 \u00b0C", "temp. pref."]
rl = [label(c, a) for c, a in READOUTS]
grid = np.full((len(rl), len(traits) * 2), np.nan)
for r in cd.itertuples():
    grid[rl.index(r.readout), traits.index(r.trait) + (0 if r.model == "THP1" else 3)] = r.rho
im = axf.imshow(grid, cmap="RdBu_r", vmin=-0.6, vmax=0.6)
axf.set_yticks(range(len(rl))); axf.set_yticklabels(rl, fontsize=6)
axf.set_xticks(range(len(traits) * 2))
SHORT = {"28 \u00b0C": "28 \u00b0C", "37 \u00b0C": "37 \u00b0C", "temp. pref.": "\u0394T"}
axf.set_xticklabels([f"{SHORT[t]}\n{'THP-1' if m == 'THP1' else 'PBMC'}"
                     for m in ("THP1", "PBMC") for t in traits], fontsize=5.2)
for i in range(len(rl)):
    for j in range(len(traits) * 2):
        if not np.isnan(grid[i, j]):
            axf.text(j, i, f"{grid[i, j]:+.2f}", ha="center", va="center", fontsize=5.5)
axf.axvline(len(traits) - 0.5, color="white", lw=1.8)
axf.set_title("f  Spearman \u03c1, induction vs Galleria phenotype (10$^5$)",
              fontsize=8, loc="left")
cbf = fig.colorbar(im, ax=axf, fraction=0.030, pad=0.035)
cbf.set_label("\u03c1", fontsize=7, labelpad=6); cbf.ax.tick_params(labelsize=6)

fig.suptitle("Cytokine induction by 21 F. solani isolates in THP-1- and PBMC-derived "
             "macrophages", fontsize=9.5, x=0.01, ha="left")
fig.savefig(f"{OUT}/figure_cell_S28.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_cell_S28.png", bbox_inches="tight", dpi=600)
print(f"source comparison: {len(cw)} tests (10^5, 2 models); min P = {cw.p.min():.3f}; "
      f"range {cw.p.min():.3f}-{cw.p.max():.3f}")
print(f"phenotype correlations: {len(cd)} tests; |rho| {cd.rho.abs().min():.2f}-"
      f"{cd.rho.abs().max():.2f}; min P = {cd.p.min():.3f}")
print(f"wrote {OUT}/figure_cell_S28.pdf/.png + cell_clinical_vs_env.tsv + cell_vs_phenotype.tsv")
