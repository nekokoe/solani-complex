#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main-figure sub-panel evaluating the practical significance of the
RF-selected clinical-associated HOGs, computed separately for F. solani
(319 HOGs) and F. falciformis (196 HOGs).

Five panels (single row):
  a: fraction of the pan-genome represented by the selected HOGs (per species)
  b: per-genome representation - selected HOGs carried per genome, clinical vs
     non-clinical isolates within each training species (Mann-Whitney U)
  c: per-genome effector GENE fraction: paired within each genome, selected
     genes vs all genes of that genome (Wilcoxon signed-rank)
  d: per-genome protein-kinase GENE fraction, paired as in c
  e: clinical vs non-clinical preference: mean presence rates per species
     (Wilcoxon signed-rank) and direction counts

(The functional-composition preference panels and the multiple-testing panel
are provided as separate figures: figure_pref.clinical_vs_nonclinical.pdf.)

Inputs : figure_data.rf.tsv, rf_genes.per_genome.tsv (stats_rf_genes.py)
Output : figure_subpanel.rf_genes.pdf

Usage:
    python3 plot_rf_figure.py --figure-data 09.rf_clinical_genes/figure_data.rf.tsv \
        --per-genome 09.rf_clinical_genes/rf_genes.per_genome.tsv \
        --output 09.rf_clinical_genes/figure_subpanel.rf_genes.pdf
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
COLOR_SEL = "#C44E52"
COLOR_ALL = "#BDBDBD"
COLOR_SOLANI = "#4C72B0"
COLOR_FALC = "#DD8452"
COLOR_GRAY = "#8C8C8C"

SPECIES = ["solani", "falciformis"]
SPECIES_LABEL = {"solani": "F. solani", "falciformis": "F. falciformis"}


def parse_arguments():
    p = argparse.ArgumentParser(description="Plot RF-selected HOG evaluation sub-panel")
    p.add_argument("--figure-data", required=True)
    p.add_argument("--per-genome", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main():
    args = parse_arguments()
    fd = pd.read_csv(args.figure_data, sep="\t")
    pg = pd.read_csv(args.per_genome, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
    })

    def val(panel, species, category):
        r = fd[(fd["panel"] == panel) & (fd["species"] == species)
               & (fd["category"] == category)]
        return r["value"].iloc[0] if len(r) else np.nan

    fig, axes = plt.subplots(1, 5, figsize=(13.6, 2.9))
    plt.subplots_adjust(wspace=0.62)

    # ---------------- a: fraction of the pan-genome (per species) ----------------
    ax = axes[0]
    n_tot = int(fd[(fd["panel"] == "fraction") & (fd["species"] == "combined")
                   & (fd["category"] == "total")]["n"].iloc[0])
    species_frac = [val("fraction", sp, "selected") for sp in SPECIES]
    species_n = [int(fd[(fd["panel"] == "fraction") & (fd["species"] == sp)
                        & (fd["category"] == "selected")]["n"].iloc[0]) for sp in SPECIES]
    x = np.arange(2)
    ax.bar(x, species_frac, width=0.5, color=[COLOR_SOLANI, COLOR_FALC],
           edgecolor="white", linewidth=0.4)
    for xi, (frac, n) in enumerate(zip(species_frac, species_n)):
        ax.text(xi, frac + 0.04, f"{frac:.2f}%\n(n={n})", ha="center", va="bottom", fontsize=6.2)
    ax.set_xticks(x)
    ax.set_xticklabels([SPECIES_LABEL[s] for s in SPECIES], fontsize=6.4,
                              rotation=45, ha="right")
    ax.set_ylabel("% of pan-genome HOGs")
    ax.set_ylim(0, 1.9)
    ax.set_title(f"Selected HOGs vs pan-genome\n(total {n_tot:,} HOGs)", fontsize=7.2, pad=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.16, 1.1, "a", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")

    # ---------------- b: per-genome representation ----------------
    ax = axes[1]
    positions, labels, data, colors = [], [], [], []
    for sp in SPECIES:
        sub = pg[pg["Species"] == sp]
        for src, color, lab in [("med", COLOR_CLIN, "Clinical"),
                                ("env", COLOR_NONCLIN, "Non-clinical")]:
            positions.append(len(positions) + 1)
            labels.append(lab)
            data.append(sub.loc[sub["Source"] == src, "n_hogs_selected"].values)
            colors.append(color)
    bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                    medianprops={"color": "black", "linewidth": 0.9},
                    boxprops={"linewidth": 0.7}, whiskerprops={"linewidth": 0.7},
                    capprops={"linewidth": 0.7}, showfliers=False)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    for sp, idx in [("solani", 0), ("falciformis", 2)]:
        p = val("per_genome", sp, "clin_vs_nonclin_p")
        ytop = max(data[idx].max(), data[idx + 1].max())
        ax.plot([positions[idx], positions[idx + 1]], [ytop + 6] * 2, color="k", lw=0.6)
        ax.text(positions[idx] + 0.5, ytop + 10, stars(p), ha="center",
                fontsize=7, fontweight="bold")
        ax.text(positions[idx] + 0.5, 356, SPECIES_LABEL[sp], ha="center",
                fontsize=6.4, fontstyle="italic")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("Selected HOGs per genome")
    ax.set_ylim(100, 372)
    ax.set_title(f"Per-genome representation\n(~{val('per_genome', 'combined', 'frac_genes'):.1f}% of genome genes)",
                 fontsize=7.2, pad=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.16, 1.1, "b", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")

    # ---------------- c/d: per-genome paired effector / kinase (gene level) ----------------
    for panel_idx, (metric, panel_name, ylab) in enumerate(
            [(("SelectedEffectorFrac", "AllEffectorFrac"), "effector_paired", "Effector genes (%)"),
             (("SelectedKinaseFrac", "AllKinaseFrac"), "kinase_paired", "Kinase genes (%)")]):
        ax = axes[2 + panel_idx]
        positions, labels, data, colors = [], [], [], []
        for sp in SPECIES:
            sub = pg[pg["Species"] == sp]
            for col, color, lab in [(metric[0], COLOR_SEL, "Selected"),
                                    (metric[1], COLOR_ALL, "All genes")]:
                positions.append(len(positions) + 1)
                labels.append(lab)
                data.append(sub[col].dropna().values * 100)
                colors.append(color)
        bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                        medianprops={"color": "black", "linewidth": 0.9},
                        boxprops={"linewidth": 0.7}, whiskerprops={"linewidth": 0.7},
                        capprops={"linewidth": 0.7}, showfliers=False)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ymax_all = max(d.max() for d in data)
        ylim_top = ymax_all * 1.55
        for sp, idx in [("solani", 0), ("falciformis", 2)]:
            p = val(panel_name, sp, "p")
            ytop = max(data[idx].max(), data[idx + 1].max())
            ax.plot([positions[idx], positions[idx + 1]], [ytop + 0.15] * 2, color="k", lw=0.6)
            ax.text(positions[idx] + 0.5, ytop + 0.38, stars(p), ha="center",
                    fontsize=7, fontweight="bold")
            ax.text(positions[idx] + 0.5, ylim_top - 0.55, SPECIES_LABEL[sp], ha="center",
                    fontsize=6.4, fontstyle="italic")
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, ylim_top)
        ax.set_title(panel_name.replace("_", " ").title(), fontsize=7.2, pad=3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.16, 1.1, "cd"[panel_idx], transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="right")
        if panel_idx == 0:
            ax.legend(handles=[Patch(facecolor=COLOR_SEL, label="Selected genes"),
                               Patch(facecolor=COLOR_ALL, label="All genes of genome")],
                      loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
                      fontsize=5.4, frameon=False)

    # ---------------- e: preference per species ----------------
    ax = axes[4]
    positions, labels, data, colors = [], [], [], []
    for sp in SPECIES:
        for rate_col, color, lab in [("mean_clinical_rate", COLOR_CLIN, "Clinical"),
                                     ("mean_nonclinical_rate", COLOR_NONCLIN, "Non-clinical")]:
            positions.append(len(positions) + 1)
            labels.append(lab)
            data.append(val("preference", sp, rate_col) * 100)
            colors.append(color)
    ax.bar(positions, data, width=0.5, color=colors, edgecolor="white", linewidth=0.4)
    for pi, d in zip(positions, data):
        ax.text(pi, d + 1.5, f"{d:.0f}%", ha="center", fontsize=5.8)
    for sp, idx in [("solani", 0), ("falciformis", 2)]:
        p = val("preference", sp, "signedrank_p")
        cb = int(val("preference", sp, "clinical_biased"))
        nb = int(val("preference", sp, "nonclinical_biased"))
        ax.plot([positions[idx], positions[idx + 1]], [62] * 2, color="k", lw=0.6)
        ax.text(positions[idx] + 0.5, 65.5, stars(p), ha="center", fontsize=7, fontweight="bold")
        ax.text(positions[idx] + 0.5, 74, SPECIES_LABEL[sp], ha="center",
                fontsize=6.4, fontstyle="italic")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=6, rotation=45, ha="right")
    ax.set_ylabel("Mean presence rate (%)")
    ax.set_ylim(0, 80)
    ax.set_title("Clinical preference\n(clinical- / non-clinical-biased)", fontsize=7.0, pad=3)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.16, 1.1, "e", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
