#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Redesigned main-figure sub-panel for the RF-selected clinical-associated HOGs
(F. solani: 319, F. falciformis: 196), computed per species.

Panels (2 x 2):
  a: two pie charts showing the selected HOGs as a tiny fraction of the
     pan-genome (29,080 HOGs) for each species
  b: histogram of the per-HOG prevalence (fraction of CLINICAL isolates of the
     species cohort carrying the HOG), F. falciformis vs F. solani
  c: same histogram for NON-CLINICAL isolates
  d: functional composition (effector / kinase / transporter / regulatory) of
     the clinical-dominant vs non-clinical-dominant selected HOGs, per species

Inputs : suppl_table.rf_prevalence.tsv, suppl_table.rf_dom_functions.tsv,
         figure_data.rf_v2.tsv (stats_rf_prevalence.py)
Output : figure_subpanel.rf_genes.pdf

Usage:
    python3 plot_rf_v2.py \
        --prevalence 09.rf_clinical_genes/suppl_table.rf_prevalence.tsv \
        --dom-functions 09.rf_clinical_genes/suppl_table.rf_dom_functions.tsv \
        --figure-data 09.rf_clinical_genes/figure_data.rf_v2.tsv \
        --output 09.rf_clinical_genes/figure_subpanel.rf_genes.pdf
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLOR_FALC = "#DD8452"
COLOR_SOLANI = "#4C72B0"
COLOR_CLIN = "#C44E52"
COLOR_NONCLIN = "#55A868"
FUNC_ORDER = ["Effector", "Kinase", "Transporter", "Regulatory", "Others"]
FUNC_LABEL = {"Effector": "Effector", "Kinase": "Kinase",
              "Transporter": "Transporter", "Regulatory": "Regulatory",
              "Others": "Others"}


def parse_arguments():
    p = argparse.ArgumentParser(description="Plot redesigned RF-selected HOG sub-panel")
    p.add_argument("--prevalence", required=True)
    p.add_argument("--dom-functions", required=True)
    p.add_argument("--figure-data", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()


def main():
    args = parse_arguments()
    prev = pd.read_csv(args.prevalence, sep="\t")
    dom = pd.read_csv(args.dom_functions, sep="\t")
    fd = pd.read_csv(args.figure_data, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 7.4))
    axes = axes.ravel()
    # symmetric square figure + symmetric margins => ~1:1 square subplots
    plt.subplots_adjust(left=0.10, right=0.97, bottom=0.09, top=0.95,
                        wspace=0.45, hspace=0.45)

    # ---------------- a: two pies (selected fraction of pan-genome) ----------------
    ax = axes[0]
    for i, sp in enumerate(["falciformis", "solani"]):
        sel = fd[(fd["panel"] == "pie") & (fd["species"] == sp) & (fd["category"] == "selected")]
        tot = fd[(fd["panel"] == "pie") & (fd["species"] == sp) & (fd["category"] == "total")]
        n_sel, n_tot = int(sel["n"].iloc[0]), int(tot["n"].iloc[0])
        pct = n_sel / n_tot * 100
        cx = -0.55 if sp == "falciformis" else 0.55
        ax.pie([n_sel, n_tot - n_sel],
               colors=[COLOR_FALC if sp == "falciformis" else COLOR_SOLANI, "#EDEDED"],
               startangle=90, counterclock=False,
               center=(cx, 0), radius=0.38,
               wedgeprops={"edgecolor": "white", "linewidth": 0.6})
        ax.text(cx, 0, f"{pct:.2f}%", ha="center", va="center",
                fontsize=9.5, fontweight="bold")
        ax.text(cx, -0.52, f"F. {'falciformis' if sp=='falciformis' else 'solani'}\n"
                f"{n_sel:,} / {n_tot:,} HOGs", ha="center", va="top", fontsize=7.2)
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-0.72, 0.6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Selected HOGs as fraction\nof the pan-genome", fontsize=8.6, pad=3)
    ax.text(-0.12, 1.08, "a", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")

    # ---------------- b/c: prevalence histograms (within-species comparison) ----------------
    # b = F. falciformis, c = F. solani; each panel compares the clinical vs
    # non-clinical prevalence distribution of the selected HOGs within that
    # species, using the species' colour family (dark = clinical, light =
    # non-clinical; same palette as panel e). Bins are fixed to 0-0.1, 0.1-0.2,
    # ..., 0.9-1.0 and labelled with the bin ranges.
    bins = np.arange(0, 1.05, 0.1)
    sp_palette = {
        "falciformis": {"clinical": "#DD8452", "non-clinical": "#F5C9A0"},
        "solani": {"clinical": "#4C72B0", "non-clinical": "#A9C4E8"},
    }
    for panel_idx, sp in enumerate(["falciformis", "solani"]):
        ax = axes[1 + panel_idx]
        data_clin = prev.loc[prev["Species"] == sp, "Clinical_rate"].values
        data_nonclin = prev.loc[prev["Species"] == sp, "NonClinical_rate"].values
        ax.hist([data_clin, data_nonclin], bins=bins, histtype="bar",
                color=[sp_palette[sp]["clinical"], sp_palette[sp]["non-clinical"]],
                edgecolor="white", linewidth=0.4,
                label=["Clinical isolates", "Non-clinical isolates"])
        for data, color in [(data_clin, sp_palette[sp]["clinical"]),
                            (data_nonclin, sp_palette[sp]["non-clinical"])]:
            ax.axvline(data.mean(), color=color, linestyle="--", linewidth=0.9)
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_xlabel("Prevalence (fraction of isolates)")
        ax.set_ylabel("Selected HOGs")
        ax.set_ylim(0, 100)
        ax.set_title(f"Prevalence in clinical vs non-clinical\nF. {sp}",
                     fontsize=7.2, pad=3)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.12, 1.08, "bc"[panel_idx], transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="top", ha="right")
        if panel_idx == 0:
            ax.legend(loc="upper right", fontsize=5.4, frameon=False)

    # ---------------- e: functional composition by dominance group ----------------
    # Absolute HOG counts; species use different colour families (falciformis
    # orange, solani blue, as in panels b/c) with clinical/non-clinical as
    # paired shades within each family.
    ax = axes[3]
    x = np.arange(len(FUNC_ORDER))
    bar_w = 0.2
    sp_colors = {
        "falciformis": {"clinical": "#DD8452", "non-clinical": "#F5C9A0"},
        "solani": {"clinical": "#4C72B0", "non-clinical": "#A9C4E8"},
    }
    for i, sp in enumerate(["falciformis", "solani"]):
        for j, domg in enumerate(["clinical", "non-clinical"]):
            row = dom[(dom["Species"] == sp) & (dom["Dominance"] == domg)]
            vals = [int(row.loc[row["Function"] == f, "n_Function"].iloc[0]) if len(row.loc[row["Function"] == f]) else 0
                    for f in FUNC_ORDER]
            off = (i * 2 + j - 1.5) * bar_w
            ax.bar(x + off, vals, width=bar_w, color=sp_colors[sp][domg],
                   edgecolor="white", linewidth=0.4,
                   label=f"F. {sp} {'clinical' if domg=='clinical' else 'non-clinical'}-dominant")
    ax.set_xticks(x)
    ax.set_xticklabels([FUNC_LABEL[f] for f in FUNC_ORDER], fontsize=6.4)
    ax.set_ylabel("Selected HOGs (count)")
    ax.set_ylim(0, 150)
    ax.set_title("Functional composition of clinical- vs\nnon-clinical-dominant HOGs",
                 fontsize=7.2, pad=3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2,
              fontsize=5.0, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.12, 1.08, "d", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="right")

    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"saved -> {args.output}")


if __name__ == "__main__":
    main()
