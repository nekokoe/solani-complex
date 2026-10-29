#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S25 — phenotype summary figure (nature-figure, python/matplotlib).

Top two rows: survival curves for all six conditions (28Ce3/e4/e5, 37Ce3/e4/e5;
15 strains each + bold mean).
Bottom row: (g) all-condition overlay of mean curves, (h) mortality by
condition (explicit labels, e.g. '28 °C, 1x10^3'), (i) concentration response
with per-temperature statistics.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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
OUT = "07.lethal_groups"
C28, C37 = "#3B7DD8", "#C44E52"
CONDS = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
LABELS = {"28Ce3": "28 \u00b0C, 1\u00d710\u00b3", "28Ce4": "28 \u00b0C, 1\u00d710\u2074",
          "28Ce5": "28 \u00b0C, 1\u00d710\u2075", "37Ce3": "37 \u00b0C, 1\u00d710\u00b3",
          "37Ce4": "37 \u00b0C, 1\u00d710\u2074", "37Ce5": "37 \u00b0C, 1\u00d710\u2075"}

d = pd.read_csv(f"{OUT}/caterpillar_37C_narrow.txt", sep="\t")
d["cond"] = d.Measure.str.extract(r"(\d+C[eE]\d+)")[0]
d["day"] = d.Measure.str.extract(r"day(\d+)")[0].astype(int)
mort = pd.read_csv(f"{OUT}/phenotype_analysis/phenotype_mortality.tsv",
                    sep="\t", index_col=0)
STRAINS = mort.index.tolist()
days = sorted(d.day.unique())

fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.2), dpi=300)
surv_axes = axes[0:2].flatten()
bot = axes[2]

# ---------- top: survival curves per condition ----------
for i, (ax, cond) in enumerate(zip(surv_axes, CONDS)):
    sub = d[d.cond == cond]
    for s in STRAINS:
        dd = sub[sub.Sample == s].set_index("day").Percent
        ax.plot(dd.index, dd.values, "o-", ms=2, lw=0.6, alpha=0.35,
                color="grey")
    mean = sub.groupby("day").Percent.mean()
    col = C28 if cond.startswith("28") else C37
    ax.plot(mean.index, mean.values, "-o", color=col, lw=1.8, ms=3)
    ax.set_title(f"{'abcdef'[i]}  {LABELS[cond]}", fontsize=7.5, loc="left")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.02)
    ax.tick_params(labelsize=5.5)
    if ax not in surv_axes[3:]:
        ax.set_xticklabels([])
    if ax not in surv_axes[[0, 3]]:
        ax.set_yticklabels([])
for i, ax in enumerate(surv_axes):
    if i >= len(CONDS):
        ax.axis("off")
surv_axes[3].set_ylabel("Survival", fontsize=7)
surv_axes[3].set_xlabel("Day", fontsize=7)

# ---------- g: overlay of all condition means ----------
ax = bot[0]
for cond in CONDS:
    mean = d[d.cond == cond].groupby("day").Percent.mean()
    col = C28 if cond.startswith("28") else C37
    ls = ["-", "--", ":"][int(cond[-1]) - 3]
    ax.plot(mean.index, mean.values, ls, color=col, lw=1.4, label=LABELS[cond])
ax.set_xlim(0, 10); ax.set_ylim(0, 1.02)
ax.set_xlabel("Day", fontsize=7)
ax.set_ylabel("Mean survival", fontsize=7)
ax.set_title("g  All-condition mean survival", fontsize=8, loc="left")
ax.legend(fontsize=5.2, ncol=2, loc="upper right")
ax.tick_params(labelsize=5.5)

# ---------- h: mortality by condition (explicit labels) ----------
# clinical (med) vs environmental (env) strain colouring
ENV_STRAINS = {"BMU12683", "BMU12684", "BMU12685", "BMU12686"}
ax = bot[1]
xs = np.arange(len(CONDS))
for i, cond in enumerate(CONDS):
    vals = mort[f"mort_{cond}"].dropna()
    col = C28 if cond.startswith("28") else C37
    # clinical: filled dots in temp colour; environmental: open grey circles
    for s, v in vals.items():
        env = s in ENV_STRAINS
        ax.scatter(i + np.random.uniform(-0.16, 0.16), v, s=20 if not env else 26,
                   c="grey" if env else col, alpha=0.55, zorder=3,
                   marker="o", facecolors="none" if env else col,
                   edgecolors="grey" if env else "none", linewidths=0.8)
    ax.plot([i - 0.25, i + 0.25], [vals.mean(), vals.mean()], color="black",
            lw=1.6, zorder=5)
ax.set_xticks(xs)
ax.set_xticklabels([LABELS[c] for c in CONDS], rotation=45, ha="right", fontsize=6)
ax.set_ylabel("Day-10 mortality", fontsize=7)
ax.set_title("h  Mortality by condition", fontsize=8, loc="left")
ax.set_ylim(0, 1.05)
ax.tick_params(labelsize=5.5)
# legend: clinical vs environmental
from matplotlib.lines import Line2D
leg_h = [Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
                markersize=5, label="Clinical"),
         Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                markeredgecolor="grey", markersize=5, label="Environmental")]
ax.legend(handles=leg_h, fontsize=5.2, loc="upper right", frameon=False)
# concentration effect p values
for temp, xl in [("28", [0, 1, 2]), ("37", [3, 4, 5])]:
    groups = [mort[f"mort_{CONDS[i]}"].dropna().values for i in xl]
    p = stats.f_oneway(*groups).pvalue
    ax.text(np.mean(xl), 1.03, f"P = {p:.0e}", ha="center", fontsize=5.5)

# ---------- i: concentration response ----------
ax = bot[2]
conc = [3, 4, 5]
for temp, col, cs in [("28", C28, ["28Ce3", "28Ce4", "28Ce5"]),
                      ("37", C37, ["37Ce3", "37Ce4", "37Ce5"])]:
    for s in STRAINS:
        y = mort[[f"mort_{c}" for c in cs]].loc[s].values
        ax.plot(conc, y, "-o", color=col, lw=0.5, ms=2, alpha=0.3)
    m = mort[[f"mort_{c}" for c in cs]].mean()
    ax.plot(conc, m.values, "-o", color=col, lw=2, ms=6,
            label=f"{temp} \u00b0C")
ax.set_xticks(conc)
ax.set_xticklabels(["1\u00d710\u00b3", "1\u00d710\u2074", "1\u00d710\u2075"])
ax.set_xlabel("Inoculum (conidia per larva)", fontsize=7)
ax.set_ylabel("Day-10 mortality", fontsize=7)
ax.set_title("i  Concentration response", fontsize=8, loc="left")
ax.legend(fontsize=6, loc="upper left")
ax.set_ylim(0, 1.05)
ax.tick_params(labelsize=5.5)
for temp, cs in [("28", ["28Ce3", "28Ce4", "28Ce5"]),
                 ("37", ["37Ce3", "37Ce4", "37Ce5"])]:
    groups = [mort[f"mort_{c}"].dropna().values for c in cs]
    p = stats.f_oneway(*groups).pvalue
    ax.text(5.0, 0.05 if temp == "28" else 0.13, f"{temp} \u00b0C P = {p:.0e}",
            ha="right", fontsize=5.5, color=C28 if temp == "28" else C37)

fig.subplots_adjust(left=0.055, right=0.945, bottom=0.06, top=0.94,
                    wspace=0.35, hspace=0.45)
fig.savefig(f"{OUT}/figure_phenotype_S25.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_phenotype_S25.png", dpi=300, bbox_inches="tight")
print("S25 saved")
