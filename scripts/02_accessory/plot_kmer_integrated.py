#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive 3x2 integrated figure for the Q1-Q6 kmer-based AC distribution
analysis. Designed for publication presentation (Nature-style multi-panel).

Panels:
  a - Method flowchart (compact)
  b - Q1: alluvial / connection plot: kmer-type x gene-transfer class independence
  c - Q2: dual-domain scatter with marginal histograms + type2 subtype inset
  d - Q3/Q4: clinical association (subtype bars + per-AC lollipop)
  e - Q5: RF HOG location (grouped bars with OR annotation)
  f - Q6: BMU02483 AC carrier profile (horizontal bar + inset RF enrichment)
"""

import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.patheffects import withStroke
import numpy as np
import pandas as pd
from scipy import stats

K = "06.accessory-chrs/kmer_analysis"
SH = "06.accessory-chrs/tgs_assemblies"
GT = "06.accessory-chrs/gene_transfer_211"
OUT = K

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                     "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
C_T2, C_T1, C_CORE, C_MED, C_ENV = "#D7432D", "#3B7DD8", "#6B6B6B", "#C44E52", "#55A868"
C_GOLD, C_LIGHT, C_WARN = "#F5A623", "#E8E8E8", "#8B4513"


def load():
    mat = pd.read_csv(f"{K}/kmer_matrix.ACx211.tsv", sep="\t", index_col=0)
    v3 = pd.read_csv(f"{K}/kmer_type.classification.v3.tsv", sep="\t")
    dd = pd.read_csv(f"{K}/kmer_type.dual_domain.tsv", sep="\t")
    sub = pd.read_csv(f"{K}/kmer_type2.subtypes.tsv", sep="\t")
    gt = pd.read_csv(f"{GT}/suppl_table.ac_transfer_profile.tsv", sep="\t")
    ca = pd.read_csv(f"{K}/kmer_clinical_ac.assessment.tsv", sep="\t")
    med = {}
    for line in open(f"{SH}/bmu_ref_groups.211.txt"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "BMU_ID":
            med[p[0]] = p[1]
    sp_map = {}
    for line in open(f"{SH}/species.211.map.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            sp_map[p[0]] = p[1]
    return mat, v3, dd, sub, gt, ca, sp_map, med


def panel_a(ax):
    """Compact method flowchart."""
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    def box(x, y, w, h, text, fc, fs=5.8, lw=0.7):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.08",
                           linewidth=lw, edgecolor="#444", facecolor=fc, zorder=2)
        ax.add_patch(b)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                zorder=3, linespacing=1.3)
    def arr(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=10, linewidth=0.8,
                                     color="#666", zorder=1))
    C_IN, C_WF, C_OUT = "#DCE9F7", "#FFF9E6", "#E8F4E1"
    box(1.5, 8.0, 7.0, 1.2, "78 ACs (13 TGS) \u00d7 211 isolates\nshared-kmer proportion matrix", C_IN, fs=6.2)
    arr(5, 8.0, 5, 7.5)
    box(1.5, 6.2, 7.0, 1.2, "1  kmer Type I/II classification\n  (species-coherence: single-species vs mosaic)", C_WF, fs=6.2)
    arr(5, 6.2, 2.5, 5.7); arr(5, 6.2, 7.5, 5.7)
    box(0.2, 4.4, 4.6, 1.2, "2  dual-domain: intra-species\n  vs whole-SC sharing breadth", C_WF, fs=6.0)
    box(5.2, 4.4, 4.6, 1.2, "3  type2 subtypes + clinical\n  association (med vs env carriers)", C_WF, fs=6.0)
    arr(2.5, 4.4, 2.5, 3.9); arr(7.5, 4.4, 7.5, 3.9)
    box(0.2, 2.6, 4.6, 1.2, "4  clinical AC identification\n  (burden + species-mix corrected)", C_WF, fs=6.0)
    box(5.2, 2.6, 4.6, 1.2, "5  RF clinical gene location\n  (type2 vs core: both contribute)", C_WF, fs=6.0)
    arr(2.5, 2.6, 5, 2.1); arr(7.5, 2.6, 5, 2.1)
    box(1.5, 0.8, 7.0, 1.2,
        "key findings: kmer types are orthogonal to gene-transfer classes"
        "\ntype2 ACs carry clinical-associated cargo; ~45% RF genes core-resident"
        "\nBMU02483 acquired near-intact clinical type2 ACs (5.9\u00d7 RF enrichment)", C_OUT, fs=5.8)
    ax.text(0.05, 0.97, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


def panel_b(ax, v3, gt):
    """Q1: connection/alluvial showing kmer-type vs gene-transfer type independence."""
    m = v3[["AC", "kmer_type"]].merge(gt[["AC", "AC_class"]], on="AC")
    ct = pd.crosstab(m.kmer_type, m.AC_class)
    labels_k = ["kmer Type I\n(29 ACs)", "kmer Type II\n(49 ACs)"]
    labels_g = ["block-transfer\n(15 ACs)", "scattered-transfer\n(63 ACs)"]
    counts = np.array(ct)
    n_k = counts.sum(axis=1)
    n_g = counts.sum(axis=0)

    # left bars (kmer types)
    left_x = 0; bar_w = 0.8
    y_cum_k = [0]
    for i in range(len(n_k)):
        ya = y_cum_k[-1] + 0.3
        y_cum_k.append(ya + n_k[i])
        ax.bar(left_x, n_k[i], bar_w, bottom=ya, color=[C_T1, C_T2][i],
               edgecolor="white", linewidth=0.8, label=labels_k[i].split("\n")[0] if i == 0 else None)
    # right bars (gene transfer types)
    right_x = 2.0
    y_cum_g = [0]
    for j in range(len(n_g)):
        ya = y_cum_g[-1] + 0.3
        y_cum_g.append(ya + n_g[j])
        ax.bar(right_x, n_g[j], bar_w, bottom=ya,
               color=[C_CORE, "#B0B0B0"][j], edgecolor="white", linewidth=0.8, alpha=0.7)

    # connecting ribbons (flows)
    from matplotlib.patches import Polygon
    for i in range(len(n_k)):
        for j in range(len(n_g)):
            v = counts[i, j]
            if v == 0: continue
            yl = y_cum_k[i] + counts[i, :j].sum()
            yr = y_cum_g[j] + counts[:i, j].sum()
            # proportionally scale the ribbon width
            frac = v / max(1, n_k[i] if i == 0 else n_k[i])
            yl0 = yl; yl1 = yl + v
            yr0 = yr; yr1 = yr + v
            ax.fill_between([left_x + bar_w, right_x], [yl0, yr0], [yl1, yr1],
                            color="#AAAAAA", alpha=0.45, linewidth=0)
    # labels
    ax.text(left_x + bar_w / 2, 85, labels_k[0], ha="center", fontsize=6, fontweight="bold")
    ax.text(left_x + bar_w / 2, 35, labels_k[1], ha="center", fontsize=6, fontweight="bold")
    ax.text(right_x + bar_w / 2, 70, "block", ha="center", fontsize=6)
    ax.text(right_x + bar_w / 2, 30, "scattered", ha="center", fontsize=6)
    ax.text(left_x + bar_w / 2, -3, "kmer\nclassification", ha="center", fontsize=6.5, fontweight="bold")
    ax.text(right_x + bar_w / 2, -3, "gene-transfer\nclassification", ha="center", fontsize=6.5, fontweight="bold")
    ax.set_xlim(-0.5, 3.5); ax.set_ylim(-10, 95)
    ax.axis("off")
    ax.text(0.02, 0.97, "b  Q1: kmer vs gene-transfer classes: independent (kappa = 0.09, P = 0.58)",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")
    ax.text(0.02, 0.89, f"blue = kmer type I  |  red = kmer type II\ndark bars = block-transfer  |  light = scattered",
            transform=ax.transAxes, fontsize=5.8, va="top")


def panel_c(ax_main, dd, sub):
    """Q2: dual-domain scatter + marginal histograms + type2 subtype inset."""
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    # main scatter
    ax = ax_main
    colorm = {"type1_species_specific": C_T1, "type2_SC_wide": C_T2,
              "rare_isolate_specific": "#B0B0B0"}
    for dm, c in colorm.items():
        s = dd[dd.dominant_mode == dm]
        ax.scatter(s["intra_breadth_0.4"], s["inter_breadth_0.4"], s=55, c=c,
                   alpha=0.82, edgecolors="white", linewidths=0.3, zorder=3)
    ax.plot([0, 1], [0.3, 0.3], color="black", linestyle=":", linewidth=1.0)
    ax.set_xlim(-0.02, 1.05); ax.set_ylim(-0.02, 0.7)
    ax.set_xlabel("Intra-species sharing breadth", fontsize=7.5)
    ax.set_ylabel("Cross-species sharing breadth", fontsize=7.5)
    ax.text(0.6, 0.32, "SC-wide threshold", fontsize=6, color="black", va="bottom")
    # marginal histograms (built as small insets)
    # inset: type2 subtype species bars
    ax_ins = ax.inset_axes([0.52, 0.52, 0.45, 0.45])
    spec = ["solani", "falciformis", "keratoplastica", "other"]
    x = np.arange(3); w = 0.18
    cols = {"solani": C_T1, "falciformis": C_T2, "keratoplastica": "#55A868", "other": "#B0B0B0"}
    for i, sp in enumerate(spec):
        vals = [sub[sub.subtype == s][f"br_{sp}"].mean() for s in ["subtype1", "subtype2", "subtype3"]]
        ax_ins.bar(x + (i - 1.5) * w, vals, w, color=cols[sp], label=sp)
    ax_ins.set_xticks(x)
    ax_ins.set_xticklabels(["1","2","3"], fontsize=5.5)
    ax_ins.set_ylabel("sharing breadth", fontsize=5.5)
    ax_ins.set_title("type2 subtypes\nper-species breadth", fontsize=6.2, pad=2)
    ax_ins.legend(fontsize=4.5, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.3))
    ax_ins.tick_params(labelsize=5); ax_ins.spines[["top", "right"]].set_visible(False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.text(0.02, 0.97, "c  Q2: within-species vs whole-SC sharing continuum (78 ACs)",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")
    ax.legend(handles=[plt.Line2D([0], [0], color=C_T1, lw=4, label="species-specific"),
                       plt.Line2D([0], [0], color=C_T2, lw=4, label="SC-wide"),
                       plt.Line2D([0], [0], color="#B0B0B0", lw=4, label="rare")],
              fontsize=5.5, frameon=False, loc="upper left",
              bbox_to_anchor=(0.0, -0.08), ncol=3)


def panel_d(ax, mat, sub, kt, sp_map, med):
    """Q3+Q4: composite clinical association panel."""
    iso211 = [l.strip() for l in open("02.annotation/sample.solaniSC.211.list") if l.strip()]
    valid = [i for i in iso211 if med.get(i) in ("med", "env")]
    n_med_bg = sum(1 for i in valid if med[i] == "med")
    bg_frac = n_med_bg / len(valid)

    # --- left half: Q3 subtype bars ---
    ax_left = ax.inset_axes([0.02, 0.12, 0.46, 0.85])
    rows = []
    for st in ["subtype1", "subtype2", "subtype3"]:
        acs = sub[sub.subtype == st].AC
        cars = set()
        for a in acs: cars |= {c for c in valid if mat.loc[a, c] >= 0.4}
        cars = list(cars)
        rows.append({"subtype": st, "n": len(cars), "med_frac": sum(1 for c in cars if med[c] == "med") / len(cars)})
    r = pd.DataFrame(rows)
    x = np.arange(3)
    ax_left.bar(x, r.med_frac * 100, 0.55, color=[C_T2, "#F07B6C", "#F5C9A0"], edgecolor="white", linewidth=0.6)
    ax_left.axhline(bg_frac * 100, color="black", linestyle="--", linewidth=1.0)
    for i, (_, rr) in enumerate(r.iterrows()):
        ax_left.text(x[i], rr.med_frac * 100 + 1.2, f"{rr.med_frac*100:.1f}%", ha="center", fontsize=7.5)
        p = stats.fisher_exact([[int(rr.n * rr.med_frac), int(rr.n * (1 - rr.med_frac))],
                                [n_med_bg, len(valid) - n_med_bg]])[1]
        sig = "**" if p < 0.01 else ("*" if p < 0.05 else "")
        ax_left.text(x[i], rr.med_frac * 100 + 7.0, sig, ha="center", fontsize=9)
    ax_left.set_xticks(x)
    ax_left.set_xticklabels(["deep SC-wide\n(n=4)", "broad sol-falc\n(n=6)", "pair-restricted\n(n=6)"], fontsize=5.8)
    ax_left.set_ylabel("clinical carriers (%)", fontsize=7)
    ax_left.set_ylim(0, 100)
    ax_left.text(0.98, 0.95, "Q3  subtype\nclinical\nassociation", transform=ax_left.transAxes,
                 fontsize=6.8, fontweight="bold", va="top", ha="right")
    ax_left.axhline(bg_frac * 100, color="black", linestyle="--", linewidth=1.0)
    ax_left.text(2.35, bg_frac * 100 + 1, f"background {bg_frac*100:.0f}%", fontsize=5.5, ha="right")
    ax_left.tick_params(labelsize=5.5)
    ax_left.spines[["top", "right"]].set_visible(False)
    ax_left.grid(axis="y", color="#E8E8E8", linewidth=0.5)
    ax_left.set_axisbelow(True)

    # --- right half: Q4 per-AC clinical lollipop ---
    ax_right = ax.inset_axes([0.52, 0.12, 0.46, 0.85])
    sp_med = {}
    for sp in ["solani", "falciformis", "keratoplastica"]:
        iso = [i for i in valid if sp_map.get(i) == sp]
        sp_med[sp] = sum(1 for i in iso if med[i] == "med") / len(iso)
    t2 = set(kt[kt.kmer_type == "type2_mosaic"].AC)
    rows2 = []
    for ac in mat.index:
        v = mat.loc[ac]
        cars = [c for c in valid if v[c] >= 0.4]
        if len(cars) < 5: continue
        spc = Counter(sp_map.get(c) for c in cars)
        exp = sum(sp_med.get(s, 0.65) * n for s, n in spc.items()) / len(cars)
        n_med = sum(1 for c in cars if med[c] == "med")
        rows2.append({"AC": ac, "obs": n_med / len(cars), "exp": exp, "t2": ac in t2,
                      "excess": n_med / len(cars) - exp})
    a = pd.DataFrame(rows2).sort_values("obs")
    y = np.arange(len(a))
    ax_right.scatter(a.exp, y, s=10, c="#B0B0B0", zorder=2, label="species-adjusted\nexpectation")
    colors = [C_T2 if t else C_T1 for t in a.t2]
    sizes = [40 + a.excess.iloc[i] * 60 for i in range(len(a))]
    ax_right.scatter(a.obs, y, s=sizes, c=colors, alpha=0.85, edgecolors="white",
                     linewidths=0.3, zorder=3)
    for i in range(len(a)):
        if a.excess.iloc[i] > 0.25 and a.excess.iloc[i] < a.excess.max() + 0.1:
            ax_right.plot([a.exp.iloc[i], a.obs.iloc[i]], [y[i], y[i]], color="#555", linewidth=0.4)
    ax_right.set_yticks([]); ax_right.set_xlim(0.25, 1.05)
    ax_right.set_xlabel("clinical fraction of carrier isolates", fontsize=7)
    ax_right.text(0.98, 0.95, "Q4  per-AC clinical\nassessment\n(size = clinical excess)", transform=ax_right.transAxes,
                  fontsize=6.8, fontweight="bold", va="top", ha="right")
    ax_right.legend(fontsize=5.2, frameon=False, loc="lower right")
    ax_right.spines[["top", "right"]].set_visible(False)
    ax_right.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax_right.set_axisbelow(True)
    ax.tick_params(labelsize=5.5)

    ax.text(0.02, 0.97, "d  Q3+Q4: clinical association of type2 ACs",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")
    ax.axis("off")


def panel_e(ax):
    """Q5: RF HOG location - grouped bars."""
    # hardcoded from earlier analysis
    cats = ["type2 AC", "type1 AC", "core chr", "no 13-TGS"]
    clin = [133, 18, 130, 15]
    nclin = [58, 13, 83, 20]
    x = np.arange(4); w = 0.35
    b1 = ax.bar(x - w / 2, clin, w, color=C_T2, alpha=0.9, label="clinical-biased\n(n = 296)", edgecolor="white", linewidth=0.5)
    b2 = ax.bar(x + w / 2, nclin, w, color="#F5C9A0", alpha=0.8, label="non-clinical-biased\n(n = 174)", edgecolor="white", linewidth=0.5)
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=6.5)
    ax.set_ylabel("number of RF-selected HOGs", fontsize=7.5)
    ax.text(0.02, 0.97, "e  Q5: clinical-associated genes live on both type2 ACs and core chromosomes\n(type2: 44.9% vs 33.3%, OR = 1.63, P = 0.015)",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")
    ax.legend(fontsize=6, frameon=False, loc="upper right")
    ax.tick_params(labelsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    # add OR annotation
    


def panel_f(ax, mat, v3, ca):
    """Q6: BMU02483 AC carrier profile."""
    v = mat["BMU02483"]
    s = v.sort_values(ascending=False)
    s = s[s >= 0.35]
    t2 = set(v3[v3.kmer_type == "type2_mosaic"].AC)
    clin_set = set(ca[ca.FDR < 0.05].AC)
    y = np.arange(len(s))
    colors = []
    alphas = []
    for ac in s.index:
        if ac in t2 and ac in clin_set:
            colors.append(C_T2); alphas.append(1.0)
        elif ac in t2:
            colors.append("#F07B6C"); alphas.append(0.8)
        else:
            colors.append(C_T1); alphas.append(0.7)
    import matplotlib.colors as mcolors
    colors_rgba = [(*mcolors.to_rgb(c), a) for c, a in zip(colors, alphas)]
    ax.barh(y, s.values, 0.7, color=colors_rgba, edgecolor="white", linewidth=0.3)
    ax.axvline(0.4, color="black", linestyle=":", linewidth=0.9)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=0.6)
    # label top bars
    for i in range(min(10, len(s))):
        ax.text(s.values[i] + 0.012, y[i], s.index[i].replace("BMU", "B"), fontsize=5.5, va="center")
    ax.set_yticks([]); ax.set_xlim(0, 1.05)
    ax.set_xlabel("kmer sharing of BMU02483 with the AC", fontsize=7.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    # legend + stats inset
    ax.legend(handles=[plt.Line2D([0], [0], color=C_T2, lw=4, label="clinical type2"),
                       plt.Line2D([0], [0], color="#F07B6C", lw=4, label="type2"),
                       plt.Line2D([0], [0], color=C_T1, lw=4, label="type1")],
              fontsize=5.5, frameon=False, loc="upper right",
              bbox_to_anchor=(0.98, 0.98))
    ax.text(0.02, 0.97, "f  Q6: BMU02483 carries near-intact clinical type2 ACs\n   (1,319 acquired genes; 85 RF HOGs = 5.9\u00d7 enrichment)",
            transform=ax.transAxes, fontsize=7.5, fontweight="bold", va="top")


def main():
    mat, v3, dd, sub, gt, ca, sp_map, med = load()
    fig = plt.figure(figsize=(13.5, 16.5))
    gs = fig.add_gridspec(3, 2, hspace=0.42, wspace=0.35, top=0.98, bottom=0.02,
                          left=0.04, right=0.96)

    ax_a = fig.add_subplot(gs[0, 0])
    panel_a(ax_a)

    ax_b = fig.add_subplot(gs[0, 1])
    panel_b(ax_b, v3, gt)

    ax_c = fig.add_subplot(gs[1, 0])
    panel_c(ax_c, dd, sub)

    ax_d = fig.add_subplot(gs[1, 1])
    panel_d(ax_d, mat, sub, v3, sp_map, med)

    ax_e = fig.add_subplot(gs[2, 0])
    panel_e(ax_e)

    ax_f = fig.add_subplot(gs[2, 1])
    panel_f(ax_f, mat, v3, ca)

    fig.savefig(f"{OUT}/figure_kmer_integrated.pdf", dpi=300)
    fig.savefig(f"{OUT}/figure_kmer_integrated.png", dpi=200)
    print("saved figure_kmer_integrated.pdf/.png")


if __name__ == "__main__":
    main()
