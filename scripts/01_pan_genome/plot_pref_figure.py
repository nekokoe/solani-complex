#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone figure of the FUNCTIONAL PREFERENCE between clinical and
non-clinical genomes, shown as horizontal BOX PLOTS of pairwise differences.

For each functional category (Pfam class / KEGG level-2) within a set
(RF-selected genes or 0-20% accessory genes):
  - for every isolate, compute the fraction of its annotated in-scope genes
    that belong to the category (per-isolate fraction);
  - take ALL pairwise differences (clinical-isolate fraction minus
    non-clinical-isolate fraction, in percentage points);
  - plot the distribution of these pairwise differences as a horizontal box.

Categories are ordered by the median pairwise difference; the 20 most
significant categories from each side (Mann-Whitney U test on the per-isolate
fractions, from the summary table) are shown. 'No Pfam'/'No KEGG'/'Other Pfam'/
'Unmapped' are excluded.

Inputs : suppl_table.rf_clinical_preference.tsv (summary: p-values, direction),
         suppl_table.rf_preference_per_isolate.tsv (per-isolate fractions)
Output : figure_pref.clinical_vs_nonclinical.pdf

Usage:
    python3 plot_pref_figure.py \
        --pref-table 09.rf_clinical_genes/suppl_table.rf_clinical_preference.tsv \
        --per-isolate 09.rf_clinical_genes/suppl_table.rf_preference_per_isolate.tsv \
        --output 09.rf_clinical_genes/figure_pref.clinical_vs_nonclinical.pdf
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

COLOR_CLIN = "#C44E52"
COLOR_NONCLIN = "#55A868"

PANELS = [("rf_selected", "Functional preference - RF-selected genes"),
          ("pan_0_20", "Functional preference - 0-20% accessory genes")]
TOP_N = 20          # most significant categories taken from each side
EXCLUDE = {"No Pfam", "No KEGG", "Other Pfam", "Unmapped"}


def parse_arguments():
    p = argparse.ArgumentParser(description="Functional-preference box-plot figure")
    p.add_argument("--pref-table", required=True)
    p.add_argument("--per-isolate", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def main():
    args = parse_arguments()
    summary = pd.read_csv(args.pref_table, sep="\t")
    per_iso = pd.read_csv(args.per_isolate, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
    })

    # per-panel category lists + pairwise-difference arrays
    panel_data = []
    n_bars_max = 0
    for set_name, _ in PANELS:
        sub = summary[(summary["Set"] == set_name)
                      & (~summary["Category"].isin(EXCLUDE))
                      & summary["pvalue"].notna()]
        clin = sub[sub["Effect_pp"] > 0].sort_values("pvalue").head(TOP_N)
        non = sub[sub["Effect_pp"] < 0].sort_values("pvalue").head(TOP_N)
        cats = list(clin["Category"]) + list(non["Category"])
        p_map = pd.concat([clin, non]).set_index("Category")["pvalue"].to_dict()
        iso = per_iso[per_iso["Set"] == set_name]
        diffs, meds = [], []
        for cat in cats:
            sub_iso = iso[iso["Category"] == cat]
            cf = sub_iso.loc[sub_iso["Group"] == "med", "Fraction"].values
            nf = sub_iso.loc[sub_iso["Group"] == "env", "Fraction"].values
            d = (cf[:, None] - nf[None, :]).ravel() * 100
            diffs.append(d)
            meds.append(float(np.median(d)))
        order = np.argsort(meds)                    # most negative at bottom
        cats_o = [cats[i] for i in order]
        diffs_o = [diffs[i] for i in order]
        p_o = [p_map[c] for c in cats_o]
        panel_data.append((cats_o, diffs_o, p_o))
        n_bars_max = max(n_bars_max, len(cats_o))

    fig_w = 8.0
    fig_h = 0.4 + n_bars_max * 0.10
    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    plt.subplots_adjust(wspace=0.45, top=0.90, bottom=0.13)

    for ax, (_, title), (cats, diffs, ps) in zip(axes, PANELS, panel_data):
        y_pos = np.arange(len(cats))
        bp = ax.boxplot(diffs, vert=False, positions=y_pos, widths=0.62,
                        patch_artist=True, showfliers=False,
                        medianprops={"color": "black", "linewidth": 0.7},
                        boxprops={"linewidth": 0.5},
                        whiskerprops={"linewidth": 0.4},
                        capprops={"linewidth": 0.4})
        for box, d in zip(bp["boxes"], diffs):
            box.set_facecolor(COLOR_CLIN if np.median(d) > 0 else COLOR_NONCLIN)
            box.set_alpha(0.85)
        XLIM = 15.0                       # unified x range for both panels
        star_off = 0.15
        for yi, (d, p) in enumerate(zip(diffs, ps)):
            st = stars(p)
            if st:
                med = float(np.median(d))
                if med > 0:
                    sx = min(float(np.percentile(d, 99)) + star_off, XLIM - 1.2)
                    ax.text(sx, yi, st, va="center", ha="left", color="black",
                            fontsize=5.0, fontweight="bold")
                else:
                    sx = max(float(np.percentile(d, 1)) - star_off, -XLIM + 1.2)
                    ax.text(sx, yi, st, va="center", ha="right", color="black",
                            fontsize=5.0, fontweight="bold")
        ax.axvline(0, color="black", linewidth=0.4)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(cats, fontsize=5.0)
        ax.set_xlim(-XLIM, XLIM)
        ax.set_ylim(-0.5, len(cats) - 0.5)
        ax.set_xlabel("Pairwise difference (clinical \u2212 non-clinical, % of genes)", fontsize=6.5)
        ax.set_title(title, fontsize=7.6, pad=1)
        ax.spines[["top", "right"]].set_visible(False)
    # one shared legend, single horizontal row below both panels
    fig.legend(handles=[Patch(facecolor=COLOR_CLIN, label="Higher in clinical"),
                        Patch(facecolor=COLOR_NONCLIN, label="Higher in non-clinical")],
               loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=2,
               fontsize=6.5, frameon=False)

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"saved -> {args.output} (categories per panel: {[len(c) for c, _, _ in panel_data]})")


if __name__ == "__main__":
    main()
