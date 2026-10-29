#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RF-selected clinical-associated HOGs: composition pies split by species
(F. falciformis, F. solani), in the same style as figure_subpanel.020_genes.
For each species: KEGG pathway pie, Pfam functional-class pie and an
EffectorP-3 class donut (all proportions of the species' RF-selected HOG total).

Input : rf_genes.hog_stats.tsv (stats_rf_genes.py)
Output: figure_rf.composition.pdf

Usage:
    python3 plot_rf_pies.py \
        --rf-hog-stats 09.rf_clinical_genes/rf_genes.hog_stats.tsv \
        --output 09.rf_clinical_genes/figure_rf.composition.pdf
"""

import argparse
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
from matplotlib.patches import Patch

PALETTE = plt.get_cmap("tab20").colors
COLOR_OTHERS = "#B0B0B0"
SHORT = {
    "Metabolism of cofactors and vitamins": "Cofactor/vitamin met.",
    "Metabolism of other amino acids": "Other amino acids",
    "Xenobiotics biodegradation and metabolism": "Xenobiotic met.",
    "Biosynthesis of other secondary metabolites": "Other sec. met.",
    "Metabolism of terpenoids and polyketides": "Terpenoid/polyketide",
    "Glycan biosynthesis and metabolism": "Glycan met.",
    "Carbohydrate metabolism": "Carbohydrate met.",
    "Amino acid metabolism": "Amino acid met.",
    "Nucleotide metabolism": "Nucleotide met.",
    "Energy metabolism": "Energy met.",
    "Lipid metabolism": "Lipid met.",
    "Signal transduction": "Signal transduc.",
    "Transcription factor": "Transcription factor",
    "Oxidoreductase": "Oxidoreductase",
    "Protein kinase": "Protein kinase",
    "Methyltransferase": "Methyltransferase",
    "Ankyrin repeat": "Ankyrin",
    "HET domain": "HET",
    "WD40/TPR repeat": "WD40/TPR",
    "Peptidase": "Peptidase",
    "Transporter": "Transporter",
}


def shorten(lab):
    return SHORT.get(lab, lab)


def parse_arguments():
    p = argparse.ArgumentParser(description="RF-selected HOG composition pies by species")
    p.add_argument("--rf-hog-stats", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def top_categories(counts, n_top=8, merge_into_other=(), color_map=None):
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    display = [it for it in items if it[0] not in merge_into_other][:n_top]
    labels = [it[0] for it in display]
    values = [it[1] for it in display]
    rest = sum(counts.values()) - sum(values)
    if rest > 0:
        labels.append("Others")
        values.append(rest)
    colors = []
    for i, lab in enumerate(labels):
        if lab == "Others":
            colors.append(COLOR_OTHERS)
        elif color_map and lab in color_map:
            colors.append(color_map[lab])
        else:
            colors.append(PALETTE[i % len(PALETTE)])
    return labels, values, colors


def build_color_map(dfs, col, exclude=()):
    cats = []
    for df in dfs:
        cats.extend(df[df[col] != "No Pfam"][col].tolist())
    order = [c for c, _ in Counter(cats).most_common()]
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(order)}


def pie_panel(ax, df, col, exclude, title, color_map, merge_into_other=()):
    s = df[df[col] != exclude][col]
    counts = s.value_counts().to_dict()
    labels, vals, colors = top_categories(counts, n_top=8,
                                          merge_into_other=merge_into_other,
                                          color_map=color_map)
    _, _, autotexts = ax.pie(
        vals, colors=colors, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.1f}%" if p >= 6.0 else "",
        pctdistance=0.72, textprops={"fontsize": 7.8, "color": "black"},
        wedgeprops={"edgecolor": "white", "linewidth": 0.6})
    for t in autotexts:
        t.set_fontweight("bold")
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8.8, pad=3)
    ax.legend([Patch(facecolor=c) for c in colors], [shorten(l) for l in labels],
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              fontsize=7.8, frameon=False)


def effector_donut(ax, df, title):
    vc = df["EffectorClass"].value_counts().to_dict()
    order = ["non-effector", "cytoplasmic", "apoplastic", "both"]
    colors = {"non-effector": "#E8E8E8", "cytoplasmic": "#E64B35",
              "apoplastic": "#4C72B0", "both": "#984EA3"}
    vals = [float(vc.get(c, 0)) for c in order]
    total = sum(vals)
    n_eff = total - vals[0]
    _, _, autotexts = ax.pie(
        vals, colors=[colors[c] for c in order], startangle=90,
        counterclock=False, autopct=lambda p: f"{p:.1f}%" if p >= 8.0 else "",
        pctdistance=0.72, textprops={"fontsize": 7.8, "color": "black"},
        wedgeprops={"edgecolor": "white", "linewidth": 0.6, "width": 0.34})
    for t in autotexts:
        t.set_fontweight("bold")
    ax.text(0, 0, f"{n_eff / total * 100:.1f}%\neffector", ha="center",
            va="center", fontsize=7.8, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8.8, pad=3)
    ax.legend([Patch(facecolor=colors[c]) for c in order],
              [f"{c} ({vals[i] / total * 100:.1f}%)" for i, c in enumerate(order)],
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              fontsize=7.8, frameon=False)


def main():
    args = parse_arguments()
    rf = pd.read_csv(args.rf_hog_stats, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
    })

    sol = rf[rf["Species"].isin(["solani", "both"])]
    fal = rf[rf["Species"].isin(["falciformis", "both"])]

    # global colour maps per annotation source, across both species
    kegg_map = build_color_map([sol, fal], "KEGG_Category",
                               exclude={"No KEGG"})
    pfam_map = build_color_map([sol, fal], "PfamClass",
                               exclude={"No Pfam", "Other Pfam"})

    fig = plt.figure(figsize=(14.4, 8.6))
    gs = GridSpec(2, 3, hspace=0.75, wspace=0.5,
                  top=0.94, bottom=0.05, left=0.04, right=0.97)
    axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(2)]

    for row, (sp_df, sp_name) in enumerate([(fal, "F. falciformis"),
                                            (sol, "F. solani")]):
        pie_panel(axes[row][0], sp_df, "KEGG_Category", "No KEGG",
                  f"KEGG - {sp_name}", kegg_map)
        pie_panel(axes[row][1], sp_df, "PfamClass", "No Pfam",
                  f"Pfam - {sp_name}", pfam_map, merge_into_other={"Other Pfam"})
        effector_donut(axes[row][2], sp_df, f"Effector - {sp_name}")

    for row, start in [(0, "a"), (1, "d")]:
        for c, letter in enumerate(["abcdefghij"[ord(start) - 97 + c] for c in range(3)]):
            axes[row][c].text(-0.14, 1.1, letter,
                              transform=axes[row][c].transAxes, fontsize=10,
                              fontweight="bold", va="top", ha="right")

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
