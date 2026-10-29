#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main-figure sub-panel: functional composition of the 0-20% (strain-specific /
rare) FSSC pan-genome genes. Contains ONLY 0-20% content (the RF-selected
composition pies are provided separately in 09.rf_clinical_genes/).

Layout (3 rows):
  row 1: a = exclusivity, b = effector (0-20%), c = effector (clinical-only),
         d = effector (non-clinical-only)
  row 2: e = KEGG (0-20%), g = KEGG (clinical-only), h = KEGG (non-clinical-only)
  row 3: f = Pfam (0-20%), i = Pfam (clinical-only), j = Pfam (non-clinical-only)

Colours are globally consistent: each functional category has a fixed colour
across all pies of the same annotation source.

Inputs : figure_data.020.tsv (stats_020_genes.py)
Output : figure_subpanel.020_genes.pdf

Usage:
    python3 plot_020_figure.py \
        --figure-data 08.pan_020_genes/figure_data.020.tsv \
        --output 08.pan_020_genes/figure_subpanel.020_genes.pdf
"""

import argparse
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

COLOR_CLIN = "#C44E52"
COLOR_ENV = "#55A868"
COLOR_SHARED = "#8C8C8C"
COLOR_CYTO = "#E64B35"
COLOR_APO = "#4C72B0"
COLOR_BOTH = "#984EA3"
COLOR_NONEFF = "#E8E8E8"
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
    p = argparse.ArgumentParser(description="Plot 0-20% gene composition sub-panel")
    p.add_argument("--figure-data", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def draw_pie(ax, values, colors, autopct_min=4.0, pctdistance=0.72,
             donut=False, center_text=None):
    wedgeprops = {"edgecolor": "white", "linewidth": 0.6}
    if donut:
        wedgeprops["width"] = 0.34
    _, _, autotexts = ax.pie(
        values, colors=colors, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.1f}%" if p >= autopct_min else "",
        pctdistance=pctdistance, textprops={"fontsize": 6.8, "color": "black"},
        wedgeprops=wedgeprops)
    for t in autotexts:
        t.set_fontweight("bold")
    if center_text:
        ax.text(0, 0, center_text, ha="center", va="center", fontsize=6.8,
                fontweight="bold")
    ax.set_aspect("equal")


def build_color_map(fd, panels, exclude=()):
    """Global category -> colour map, ordered by total frequency across pies."""
    cats = []
    for p in panels:
        sub = fd[fd["panel"] == p]
        sub = sub[~sub["category"].isin(exclude)]
        cats.extend(sub["category"].tolist())
    order = [c for c, _ in Counter(cats).most_common()]
    return {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(order)}


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


def composition_pie(ax, rows, n_top=8, merge_into_other=(), color_map=None,
                    autopct_min=6.0, title=None):
    counts = dict(zip(rows["category"], rows["n_HOGs"].astype(float)))
    labels, vals, colors = top_categories(counts, n_top=n_top,
                                          merge_into_other=merge_into_other,
                                          color_map=color_map)
    draw_pie(ax, vals, colors, autopct_min=autopct_min)
    ax.set_title(title, fontsize=7.8, pad=3)
    ax.legend([Patch(facecolor=c) for c in colors], [shorten(l) for l in labels],
              loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              fontsize=7.6, frameon=False)


def effector_donut(ax, rows, center_text, title):
    vals = {r["category"]: float(r["n_HOGs"]) for _, r in rows.iterrows()}
    total = sum(vals.values())
    colors = {"non-effector": COLOR_NONEFF, "cytoplasmic": COLOR_CYTO,
              "apoplastic": COLOR_APO, "both": COLOR_BOTH}
    order = ["non-effector", "cytoplasmic", "apoplastic", "both"]
    draw_pie(ax, [vals[c] for c in order], [colors[c] for c in order],
             autopct_min=8.0, donut=True, center_text=center_text)
    ax.set_title(title, fontsize=7.8, pad=3)
    ax.legend([Patch(facecolor=colors[c]) for c in order],
              [f"{c} ({vals[c] / total * 100:.1f}%)" for c in order],
              loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              fontsize=7.6, frameon=False)


def main():
    args = parse_arguments()
    fd = pd.read_csv(args.figure_data, sep="\t").fillna("")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
    })

    # global colour maps per annotation source
    kegg_map = build_color_map(fd, ["kegg_primary_020", "kegg_clinonly_020",
                                    "kegg_nonclinonly_020"],
                               exclude={"no_KEGG"})
    pfam_map = build_color_map(fd, ["pfam_primary_020", "pfam_clinonly_020",
                                    "pfam_nonclinonly_020"],
                               exclude={"no_Pfam", "Other Pfam"})

    fig = plt.figure(figsize=(12.8, 9.4))
    gs = GridSpec(3, 4, hspace=1.0, wspace=0.45,
                  top=0.95, bottom=0.05, left=0.04, right=0.97)
    # row 0: a b c d ; row 1: e g h ; row 2: f i j
    ax = {
        "a": fig.add_subplot(gs[0, 0]), "b": fig.add_subplot(gs[0, 1]),
        "c": fig.add_subplot(gs[0, 2]), "d": fig.add_subplot(gs[0, 3]),
        "e": fig.add_subplot(gs[1, 0]), "g": fig.add_subplot(gs[1, 1]),
        "h": fig.add_subplot(gs[1, 2]),
        "f": fig.add_subplot(gs[2, 0]), "i": fig.add_subplot(gs[2, 1]),
        "j": fig.add_subplot(gs[2, 2]),
    }

    # ---------------- a: 0-20% exclusivity pie ----------------
    a = ax["a"]
    e = fd[(fd["panel"] == "exclusivity_by_interval") & (fd["interval"] == "0%-20%")]
    vals = [int(e.loc[e["category"] == c, "n_HOGs"].iloc[0]) for c in
            ["clinical-only", "environmental-only", "shared"]]
    draw_pie(a, vals, [COLOR_CLIN, COLOR_ENV, COLOR_SHARED], autopct_min=5.0)
    a.set_title("Carrier exclusivity of 0-20% HOGs", fontsize=8.2, pad=3)
    a.legend([Patch(facecolor=COLOR_CLIN), Patch(facecolor=COLOR_ENV),
              Patch(facecolor=COLOR_SHARED)],
             ["Clinical-only", "Non-clinical-only", "Both"],
             loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3,
             fontsize=7.8, frameon=False)

    # ---------------- b-d: effector composition donuts ----------------
    eff = fd[(fd["panel"] == "effector_class_by_interval") & (fd["interval"] == "0%-20%")]
    n_eff_total = sum(float(eff.loc[eff["category"] == c, "pct"].iloc[0]) for c in
                      ["cytoplasmic", "apoplastic", "both"])
    rows = eff[["category", "pct"]].copy()
    rows = pd.concat([rows, pd.DataFrame([{"category": "non-effector",
                                           "pct": 100.0 - n_eff_total}])],
                     ignore_index=True)
    rows["n_HOGs"] = rows["pct"]
    effector_donut(ax["b"], rows, f"{n_eff_total:.2f}%\neffector",
                   "Effector - 0-20% genes")
    for tag, lab in [("clinonly", "Clinical-only"), ("nonclinonly", "Non-clinical-only")]:
        key = "c" if tag == "clinonly" else "d"
        sub = fd[(fd["panel"] == f"effector_{tag}_020")]
        n_eff = 100.0 - float(sub.loc[sub["category"] == "non-effector", "pct"].iloc[0])
        effector_donut(ax[key], sub, f"{n_eff:.1f}%\neffector",
                       f"Effector - {lab} 0-20% HOGs")

    # ---------------- e/g/h: KEGG pies ----------------
    kegg020 = fd[fd["panel"] == "kegg_primary_020"]
    kegg020 = kegg020[kegg020["category"] != "no_KEGG"]
    composition_pie(ax["e"], kegg020, n_top=8, color_map=kegg_map,
                    autopct_min=6.0, title="KEGG - 0-20% HOGs")
    for key, tag, lab in [("g", "clinonly", "Clinical-only"),
                          ("h", "nonclinonly", "Non-clinical-only")]:
        sub = fd[(fd["panel"] == f"kegg_{tag}_020") & (fd["category"] != "no_KEGG")]
        composition_pie(ax[key], sub, n_top=8, color_map=kegg_map,
                        autopct_min=8.0, title=f"KEGG - {lab} 0-20% HOGs")

    # ---------------- f/i/j: Pfam pies ----------------
    pfam020 = fd[fd["panel"] == "pfam_primary_020"]
    pfam020 = pfam020[pfam020["category"] != "no_Pfam"]
    composition_pie(ax["f"], pfam020, n_top=8, merge_into_other={"Other Pfam"},
                    color_map=pfam_map, autopct_min=6.0, title="Pfam - 0-20% HOGs")
    for key, tag, lab in [("i", "clinonly", "Clinical-only"),
                          ("j", "nonclinonly", "Non-clinical-only")]:
        sub = fd[(fd["panel"] == f"pfam_{tag}_020") & (fd["category"] != "no_Pfam")]
        composition_pie(ax[key], sub, n_top=8, merge_into_other={"Other Pfam"},
                        color_map=pfam_map, autopct_min=8.0,
                        title=f"Pfam - {lab} 0-20% HOGs")

    # panel letters
    for i, key in enumerate("abcd"):
        ax[key].text(-0.14, 1.1, key, transform=ax[key].transAxes, fontsize=10,
                     fontweight="bold", va="top", ha="right")
    for key in "eghfij":
        ax[key].text(-0.14, 1.1, key, transform=ax[key].transAxes, fontsize=10,
                     fontweight="bold", va="top", ha="right")

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
