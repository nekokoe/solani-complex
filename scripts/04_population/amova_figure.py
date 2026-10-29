#!/usr/bin/env python3
"""Supplementary Figure S38 - the AMOVA variance partition, in the dot style of Fig. 2G.

The main text reports two different statistics for the same question and is careful not
to compare them numerically: FST gives the magnitude of the source partition, AMOVA its
significance. This figure keeps that separation and makes it readable.

  a  AMOVA: source and error variance components per species, dot size proportional to
     the component and dot colour to permutation significance; phi_ST and P annotated
  b  FST: the between-source value (filled dot) against the between-subpopulation
     values (open dots) for the same species, on one axis -- the comparison the text
     draws, made with a single estimator throughout

  /home/miniconda3/envs/torch/bin/python 04/population/amova_figure.py
  (reads suppl_table.S17_population_differentiation.xlsx)
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
XLSX = os.path.join(ROOT, "00.manuscript", "supplementary",
                    "suppl_table.S17_population_differentiation.xlsx")
OUTDIR = os.path.join(ROOT, "00.manuscript", "supplementary")

SPECIES = ["falciformis", "solani", "keratoplastica"]
NICE = {"falciformis": "F. falciformis", "solani": "F. solani",
        "keratoplastica": "F. keratoplastica"}

plt.rcParams.update({
    "font.family": ["Nimbus Sans", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def main():
    amova = pd.read_excel(XLSX, sheet_name="amova").set_index("species")
    src = pd.read_excel(XLSX, sheet_name="source_fst").set_index("species")
    sub = pd.read_excel(XLSX, sheet_name="subpopulation_fst").set_index("species")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.2, 2.9),
                                     gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.34})

    # ---------------- a: AMOVA variance components ----------------
    for r, sp in enumerate(SPECIES):
        p = amova.loc[sp, "p"]
        for comp, off in (("sigma2_source", -0.16), ("sigma2_error", 0.16)):
            v = amova.loc[sp, comp]
            source = comp == "sigma2_source"
            # only the source component carries a test, so only it is coloured by significance
            face = ("#C8442E" if p < 0.05 else "white") if source else "#e8e8e8"
            edge = "#7a2b1c" if source else "#999999"
            ax_a.scatter(off, -r, s=abs(v) * 5200, zorder=3, facecolor=face,
                         edgecolor=edge, linewidth=0.9)
            ax_a.text(off, -r - 0.26, f"{v:+.4f}", ha="center", va="top", fontsize=6.0,
                      color="#333333")
        ph = amova.loc[sp, "phi_ST"]
        pv = amova.loc[sp, "p"]
        ptxt = "P < 0.001" if pv < 0.001 else f"P = {pv:.2f}"
        ax_a.text(0.70, -r, f"φ$_{{ST}}$ = {ph:+.3f}\n{ptxt}", va="center", fontsize=6.6,
                  color="#222222")

    ax_a.set_xlim(-0.45, 1.15)
    ax_a.set_ylim(-len(SPECIES) + 0.45, 0.45)
    ax_a.set_xticks([-0.16, 0.16])
    ax_a.set_xticklabels([r"$\sigma^2_{source}$", r"$\sigma^2_{error}$"], fontsize=7)
    ax_a.set_yticks([-i for i in range(len(SPECIES))])
    ax_a.set_yticklabels([NICE[s] for s in SPECIES], fontsize=7.5)
    ax_a.set_title("AMOVA partition of core-SNP variance", fontsize=8, pad=6)
    ax_a.text(0.02, 0.97, "a", transform=ax_a.transAxes, fontsize=11, fontweight="bold",
              va="top")

    # ---------------- b: FST, source vs subpopulation ----------------
    for r, sp in enumerate(SPECIES):
        s = src.loc[sp, "mean_of_window_MEAN_FST"]
        lo, hi = sub.loc[sp, "fst_min"], sub.loc[sp, "fst_max"]
        ax_b.plot([lo, hi], [-r, -r], color="#bbbbbb", lw=1.0, zorder=1,
                  solid_capstyle="round")
        pts = [lo] if abs(hi - lo) < 1e-9 else [lo, hi]
        ax_b.scatter(pts, [-r] * len(pts), s=26, facecolor="white", edgecolor="#2C6FBB",
                     linewidth=0.9, zorder=3)
        ax_b.scatter([s], [-r], s=34, facecolor="#2C6FBB", edgecolor="#1b4a7d",
                     linewidth=0.9, zorder=4)
        span = f"{lo:.3f}" if abs(hi - lo) < 1e-9 else f"{lo:.3f}–{hi:.3f}"
        ax_b.text(max(hi, s) + 0.04, -r, f"{s:.3f}  vs  {span}", va="center", fontsize=6.4,
                  color="#222222")

    ax_b.axvline(0, color="#999999", lw=0.6, ls=(0, (3, 3)), zorder=0)
    ax_b.set_xlim(-0.12, 1.05)
    ax_b.set_ylim(-len(SPECIES) + 0.45, 0.45)
    ax_b.set_yticks([-i for i in range(len(SPECIES))])
    ax_b.set_yticklabels([])
    ax_b.set_xlabel("Weir–Cockerham FST")
    ax_b.set_title("Between-source vs between-subpopulation FST", fontsize=8, pad=6)
    ax_b.text(0.02, 0.97, "b", transform=ax_b.transAxes, fontsize=11, fontweight="bold",
              va="top")
    ax_b.scatter([], [], s=34, facecolor="#2C6FBB", edgecolor="#1b4a7d", label="between source groups")
    ax_b.scatter([], [], s=26, facecolor="white", edgecolor="#2C6FBB", label="between subpopulations")
    ax_b.legend(frameon=False, fontsize=5.6, loc="lower right", handletextpad=0.4,
                borderaxespad=0.2)

    fig.subplots_adjust(left=0.12, right=0.985, top=0.86, bottom=0.17)
    for ext in ("pdf", "png"):
        path = os.path.join(OUTDIR, f"Figure-S38-AMOVA.{ext}")
        fig.savefig(path, dpi=400 if ext == "png" else None)
        print("wrote", path)


if __name__ == "__main__":
    main()
