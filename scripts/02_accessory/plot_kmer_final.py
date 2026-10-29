#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Final publication figures for the Q1-Q6 kmer AC distribution analysis.

3 composite figures + 1 small method chart (Nature-family figure contract):
  F1 (Q1+Q2): classification landscape — quadrant scatter + type2 subtype bars
  F2 (Q3+Q4): clinical association (hero) — per-AC lollipop + subtype bars + volcano
  F3 (Q5+Q6): clinical consequence — RF HOG location bars + BMU02483 profile
  Suppl: compact method flowchart

Color policy: one neutral (gray #6B6B6B), one signal (blue #3B7DD8 = type1),
one accent (red #D7432D = type2/clinical). All panels share the same palette.
Every panel annotates n.

Output: 06.accessory-chrs/kmer_analysis/figure_kmer_F{1,2,3}.pdf/.png
"""

import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy import stats

K, SH, GT = "06.accessory-chrs/kmer_analysis", "06.accessory-chrs/tgs_assemblies", "06.accessory-chrs/gene_transfer_211"
OUT = K
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Helvetica","Arial","DejaVu Sans"],
                     "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
C_T1, C_T2, C_GRAY = "#3B7DD8", "#D7432D", "#6B6B6B"
C_T2_LIGHT, C_KERA, C_BG = "#F5C9A0", "#55A868", "#E8E8E8"


def load():
    mat = pd.read_csv(f"{K}/kmer_matrix.ACx211.tsv", sep="\t", index_col=0)
    v3 = pd.read_csv(f"{K}/kmer_type.classification.v3.tsv", sep="\t")
    dd = pd.read_csv(f"{K}/kmer_type.dual_domain.tsv", sep="\t")
    sub = pd.read_csv(f"{K}/kmer_type2.subtypes.tsv", sep="\t")
    gt = pd.read_csv(f"{GT}/suppl_table.ac_transfer_profile.tsv", sep="\t")
    ca = pd.read_csv(f"{K}/kmer_clinical_ac.assessment.tsv", sep="\t")
    med = {}; sp = {}
    for line in open(f"{SH}/bmu_ref_groups.211.txt"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "BMU_ID": med[p[0]] = p[1]
    for line in open(f"{SH}/species.211.map.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2: sp[p[0]] = p[1]
    return mat, v3, dd, sub, gt, ca, sp, med


def letter(ax, s, x=-0.08, y=1.06):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


# ===================== F1: classification landscape =====================
def f1(mat, v3, dd, gt, sub):
    m = v3[["AC","kmer_type"]].merge(dd[["AC","specificity"]], on="AC")
    m = m.merge(gt[["AC","AC_class","genes_per_partner"]], on="AC")
    fig = plt.figure(figsize=(13, 5.6))
    gs = fig.add_gridspec(1, 2, wspace=0.38)

    # ---- a: quadrant scatter ----
    ax = fig.add_subplot(gs[0, 0])
    X_B, Y_B = 2.0, 0.2
    markers = {"Type2_scattered":"o", "Type1_block":"^"}
    for typ, c in [("type1_stable",C_T1), ("type2_mosaic",C_T2)]:
        for gcls, mk in markers.items():
            s = m[(m.kmer_type==typ)&(m.AC_class==gcls)]
            ax.scatter(s.genes_per_partner, s.specificity, s=55, c=c, marker=mk,
                       alpha=0.78, edgecolors="white", linewidths=0.3, zorder=3)
    ax.axvline(X_B, color="black", linestyle="--", linewidth=0.9)
    ax.axhline(Y_B, color="black", linestyle="--", linewidth=0.9)
    ax.set_xscale("log"); ax.set_xlim(0.1, 60); ax.set_ylim(-0.45, 0.85)
    ax.set_xlabel("Gene-transfer concentration (genes per partner, log)", fontsize=7.5)
    ax.set_ylabel("Kmer species-specificity (intra \u2212 inter sharing)", fontsize=7.5)
    ax.set_title("a  Two orthogonal classifications of 78 accessory chromosomes", fontsize=8.8, pad=3)
    ax.text(0.12, 0.78, "block\n(scatter)", fontsize=6, va="top", ha="center", color=C_T1)
    ax.text(30, 0.78, "block\n(bulk)", fontsize=6, va="top", ha="center", color=C_T1)
    ax.text(0.12, -0.38, "mosaic\n(scatter)\n(classic)", fontsize=6, ha="center", color=C_T2)
    ax.text(30, -0.38, "mosaic\n(bulk)", fontsize=6, ha="center", color=C_T2)
    ax.text(0.95, 0.92, f"kappa = 0.09  P = 0.58  n = 78 ACs", fontsize=6.5, va="top", ha="right", color="#333",
            fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="#ccc", pad=3))
    ax.scatter([],[],s=40,c=C_T1,marker="o",label="kmer block"); ax.scatter([],[],s=40,c=C_T2,marker="o",label="kmer mosaic")
    ax.scatter([],[],s=40,c="grey",marker="o",label="scatter transfer"); ax.scatter([],[],s=40,c="grey",marker="^",label="bulk transfer")
    ax.legend(fontsize=5.5, frameon=False, ncol=2, loc="lower left")
    ax.spines[["top","right"]].set_visible(False); ax.grid(color=C_BG, linewidth=0.5); ax.set_axisbelow(True)

    # ---- b: type2 subtype bars ----
    ax = fig.add_subplot(gs[0, 1])
    spec = ["solani","falciformis","keratoplastica","outgroup"]
    labs_map = {"solani":"F. solani","falciformis":"F. falciforme","keratoplastica":"F. keratoplastica","outgroup":"outgroup"}
    x = np.arange(3); w = 0.18
    colors_sp = [C_T1, C_T2, C_KERA, C_GRAY]
    for i, sp_name in enumerate(spec):
        vals = [sub[sub.subtype==s][f"br_{sp_name if sp_name!='outgroup' else 'other'}"].mean()
                for s in ["subtype1","subtype2","subtype3"]]
        ax.bar(x+(i-1.5)*w, vals, w, color=colors_sp[i], label=labs_map[sp_name])
    ax.set_xticks(x); ax.set_xticklabels(["deep SC-wide\n(n = 4)","broad sol-falc\n(n = 6)","pair-restricted\n(n = 6)"], fontsize=6.5)
    ax.set_ylabel("sharing breadth per species", fontsize=7.5)
    ax.set_title("b  Mosaic subtypes by cross-species depth (n = 16 mosaic ACs)", fontsize=8.8, pad=3)
    ax.legend(fontsize=5.5, ncol=2, frameon=False)
    ax.tick_params(labelsize=6); ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", color=C_BG, linewidth=0.5); ax.set_axisbelow(True)

    fig.savefig(f"{OUT}/figure_kmer_F1_landscape.pdf", dpi=300)
    fig.savefig(f"{OUT}/figure_kmer_F1_landscape.png", dpi=200)
    plt.close(fig); print("F1 saved")


# ===================== F2: clinical association (hero) =====================
def f2(mat, v3, dd, sub, gt, ca, sp, med):
    iso211 = [l.strip() for l in open("02.annotation/sample.solaniSC.211.list") if l.strip()]
    valid = [i for i in iso211 if med.get(i) in ("med","env")]
    n_med_bg = sum(1 for i in valid if med[i]=="med"); bg_frac = n_med_bg/len(valid)
    sp_med = {}
    for s in ["solani","falciformis","keratoplastica"]:
        iso = [i for i in valid if sp.get(i)==s]
        sp_med[s] = sum(1 for i in iso if med[i]=="med")/len(iso)
    t2 = set(v3[v3.kmer_type=="type2_mosaic"].AC)

    fig = plt.figure(figsize=(13, 7.5))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.38, height_ratios=[1.5, 1])

    # ---- a (hero, spans top-left+top-center): per-AC clinical lollipop ----
    ax = fig.add_subplot(gs[0, :2])
    rows = []
    for ac in mat.index:
        v = mat.loc[ac]
        cars = [c for c in valid if v[c]>=0.4]
        if len(cars)<5: continue
        spc = Counter(sp.get(c) for c in cars)
        exp = sum(sp_med.get(s,0.65)*n for s,n in spc.items())/len(cars)
        n_med = sum(1 for c in cars if med[c]=="med")
        rows.append({"AC":ac, "obs":n_med/len(cars), "exp":exp, "t2":ac in t2, "excess":n_med/len(cars)-exp})
    a = pd.DataFrame(rows).sort_values("obs")
    y = np.arange(len(a))
    ax.scatter(a.exp, y, s=8, c=C_GRAY, zorder=2)
    colors = [C_T2 if t else C_T1 for t in a.t2]
    sizes = [35 + max(0, a.excess.iloc[i])*80 for i in range(len(a))]
    ax.scatter(a.obs, y, s=sizes, c=colors, alpha=0.85, edgecolors="white", linewidths=0.2, zorder=3)
    top = a.sort_values("excess", ascending=False).head(6)
    stagger = [0, -5, 3, -2, 6, 0]
    for j, (_, r) in enumerate(top.iterrows()):
        i = a.index.get_loc(r.name)
        ax.plot([r.exp, r.obs], [i, i], color="#555", linewidth=0.4)
        ax.text(r.obs+0.02, i+stagger[j], r.AC.split("_")[-1], fontsize=5, va="center")
    ax.set_yticks([]); ax.set_xlim(0.2, 1.05)
    ax.set_xlabel("clinical fraction of carrier isolates (kmer \u2265 0.4)", fontsize=7.5)
    ax.set_title("a  Clinical AC identification: observed vs species-adjusted expectation (n = 76 ACs)", fontsize=8.8, pad=3)
    ax.text(0.02, 0.95, f"clinical isolates carry 32 ACs vs 10 (env)\nP = 1.8\u00d710\u207b\u00b9\u00b2",
            transform=ax.transAxes, fontsize=6, va="top", ha="left", color="#555",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="#ccc", pad=3))
    ax.spines[["top","right"]].set_visible(False); ax.grid(axis="x", color=C_BG, linewidth=0.5); ax.set_axisbelow(True)

    # ---- b: subtype clinical bars (top-right) ----
    ax = fig.add_subplot(gs[0, 2])
    rows2 = []
    for st in ["subtype1","subtype2","subtype3"]:
        acs = sub[sub.subtype==st].AC
        cars = set(); [cars.update({c for c in valid if mat.loc[a,c]>=0.4}) for a in acs]
        cars = list(cars)
        rows2.append({"subtype":st, "n":len(cars), "med_frac":sum(1 for c in cars if med[c]=="med")/len(cars)})
    r = pd.DataFrame(rows2)
    x = np.arange(3)
    ax.bar(x, r.med_frac*100, 0.55, color=[C_T2,C_T2_LIGHT,"#FCE4D6"], edgecolor="white", linewidth=0.5)
    ax.axhline(bg_frac*100, color="black", linestyle="--", linewidth=1.0)
    for i,(_,rr) in enumerate(r.iterrows()):
        ax.text(x[i], rr.med_frac*100+1, f"{rr.med_frac*100:.1f}%", ha="center", fontsize=7)
        ax.text(x[i], 3, f"n={rr.n}", ha="center", fontsize=5.5, color="white")
    ax.set_xticks(x); ax.set_xticklabels(["deep SC-wide","broad sol-falc","pair-restricted"], fontsize=6)
    ax.set_ylabel("% clinical carriers", fontsize=7); ax.set_ylim(0, 95)
    ax.set_title("b  Mosaic subtypes: clinical association", fontsize=7.8, pad=3)
    ax.tick_params(labelsize=5.5); ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", color=C_BG, linewidth=0.5); ax.set_axisbelow(True)

    # ---- c: functional volcano (bottom, full width) ----
    ax = fig.add_subplot(gs[1, :])
    func = pd.read_csv(f"{GT}/suppl_table.intra_species_function.tsv", sep="\t")
    cats_show = ["Transporter","Oxidoreductase","HET","Effector","CAZyme","Protein kinase",
                 "Transcription factor","Peptidase","Methyltransferase","Ankyrin repeat"]
    func["neg_log_p"] = -np.log10(func.p.clip(lower=1e-300))
    func = func[func.Mode=="intra_falciformis"]
    func = func.set_index("Category")
    ax.axhline(-np.log10(0.05), color="black", linestyle=":", linewidth=0.8)
    ax.axvline(1, color=C_GRAY, linestyle="--", linewidth=0.6)
    for j, c in enumerate(cats_show):
        r = func.loc[f"Pfam_{c}"] if f"Pfam_{c}" in func.index else None
        if r is not None:
            ax.scatter(r.Odds_ratio, r.neg_log_p, s=80, c=C_T2 if r.Odds_ratio>1 else C_T1,
                       alpha=0.8, edgecolors="white", linewidths=0.3, zorder=3)
            ax.annotate(c, (r.Odds_ratio, r.neg_log_p), fontsize=6.2,
                        textcoords="offset points", xytext=(10, j*5-10))
    ax.set_xscale("log"); ax.set_xlim(0.3, 10)
    ax.set_xlabel("odds ratio (clinical subtype vs others)", fontsize=7.5)
    ax.set_ylabel("\u2212log10(P)", fontsize=7.5)
    ax.set_title("c  Functional enrichment: clinical-associated type-II cargo (F. falciformis, n = 1,447 genes)", fontsize=8.8, pad=3)
    ax.text(0.95, 0.95, "enriched: HET (OR 4.6)\ntransporter (2.0)\noxidoreductase (1.7)", transform=ax.transAxes,
            fontsize=6, va="top", ha="right", color="#555")
    ax.spines[["top","right"]].set_visible(False); ax.grid(color=C_BG, linewidth=0.5); ax.set_axisbelow(True)
    ax.tick_params(labelsize=6)

    fig.savefig(f"{OUT}/figure_kmer_F2_clinical.pdf", dpi=300)
    fig.savefig(f"{OUT}/figure_kmer_F2_clinical.png", dpi=200)
    plt.close(fig); print("F2 saved")


# ===================== F3: clinical consequence =====================
def f3(mat, v3, ca):
    # ---- Q5 RF location data (hardcoded from earlier analysis) ----
    cats = ["type2 AC", "type1 AC", "core chr", "no 13-TGS"]
    clin  = [133, 18, 130, 15]
    nclin = [58, 13, 83, 20]

    fig = plt.figure(figsize=(12, 5.8))
    gs = fig.add_gridspec(1, 2, wspace=0.38)

    # ---- a: RF HOG location grouped bars ----
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(4); w = 0.35
    ax.bar(x-w/2, clin, w, color=C_T2, alpha=0.9, label=f"clinical-biased (n = 296)", edgecolor="white", linewidth=0.5)
    ax.bar(x+w/2, nclin, w, color=C_T2_LIGHT, alpha=0.8, label=f"non-clinical-biased (n = 174)", edgecolor="white", linewidth=0.5)
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=7)
    ax.set_ylabel("number of RF-selected HOGs", fontsize=7.5)
    ax.set_title("a  Where do clinical-associated genes live? (512 RF HOGs)", fontsize=8.8, pad=3)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.tick_params(labelsize=6); ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", color=C_BG, linewidth=0.5); ax.set_axisbelow(True)
    ax.annotate("OR = 1.63\nP = 0.015", xy=(0, 135), xytext=(0.5, 155), fontsize=6.5,
                arrowprops=dict(arrowstyle="->", color="#333", lw=0.8), color=C_T2, fontweight="bold")

    # ---- b: BMU02483 profile ----
    ax = fig.add_subplot(gs[0, 1])
    v = mat["BMU02483"]
    s = v.sort_values(ascending=False); s = s[s>=0.35]
    t2_acs = set(v3[v3.kmer_type=="type2_mosaic"].AC)
    clin_set = set(ca[ca.FDR<0.05].AC)
    y = np.arange(len(s))
    colors_rgba = []
    for ac in s.index:
        if ac in t2_acs and ac in clin_set:
            colors_rgba.append((*mcolors.to_rgb(C_T2), 1.0))
        elif ac in t2_acs:
            colors_rgba.append((*mcolors.to_rgb(C_T2_LIGHT), 0.85))
        else:
            colors_rgba.append((*mcolors.to_rgb(C_T1), 0.7))
    ax.barh(y, s.values, 0.7, color=colors_rgba, edgecolor="white", linewidth=0.2)
    ax.axvline(0.4, color="black", linestyle=":", linewidth=0.9)
    ax.axvline(0.5, color="black", linestyle="--", linewidth=0.6)
    for i in range(min(10,len(s))):
        ax.text(s.values[i]+0.012, y[i], s.index[i].replace("BMU","B"), fontsize=5.5, va="center")
    ax.set_yticks([]); ax.set_xlim(0, 1.05)
    ax.set_xlabel("kmer sharing of BMU02483 with the AC", fontsize=7.2)
    ax.set_title("b  BMU02483: AC-acquired clinical genes (n = 28 carried ACs \u2265 0.35)", fontsize=8.8, pad=3)
    ax.text(0.02, 0.5, "1,319 shared genes\n819 HOGs\n85 RF candidates (5.9\u00d7 enrichment)\n53 clinical-biased RF",
            transform=ax.transAxes, fontsize=6.5, va="center", color="#333",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="#ccc", pad=4))
    ax.legend(handles=[plt.Line2D([0],[0],color=C_T2,lw=4,label="clinical type2"),
                       plt.Line2D([0],[0],color=C_T2_LIGHT,lw=4,label="type2"),
                       plt.Line2D([0],[0],color=C_T1,lw=4,label="type1")],
              fontsize=5.5, frameon=False, loc="upper left", bbox_to_anchor=(0.02, 0.98))
    ax.spines[["top","right"]].set_visible(False); ax.grid(axis="x", color=C_BG, linewidth=0.5); ax.set_axisbelow(True)

    fig.savefig(f"{OUT}/figure_kmer_F3_consequence.pdf", dpi=300)
    fig.savefig(f"{OUT}/figure_kmer_F3_consequence.png", dpi=200)
    plt.close(fig); print("F3 saved")


# ===================== Suppl: method flowchart =====================
def suppl_flowchart():
    fig, ax = plt.subplots(figsize=(8.2, 10.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 14); ax.axis("off")
    def box(x, y, w, h, text, fc, fs=6.2):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.1",
                           linewidth=0.8, edgecolor="#444", facecolor=fc, zorder=2)
        ax.add_patch(b); ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, zorder=3, linespacing=1.3)
    def arr(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2), arrowstyle="-|>", mutation_scale=10, linewidth=0.8, color="#666", zorder=1))
    ax.text(5, 13.6, "Kmer-based accessory chromosome distribution analysis", ha="center", fontsize=11, fontweight="bold")
    box(1, 12.2, 8, 0.9, "Input: 78 ACs (13 TGS) \u00d7 211 FSSC isolates shared-kmer matrix", "#DCE9F7"); arr(5,12.2,5,11.8)
    box(1, 10.5, 8, 1.2, "Q1  kmer block/mosaic classification\nspecies-coherence: single-species (block) vs mosaic", "#FFF9E6"); arr(5,10.5,5,10.1)
    box(1, 8.8, 8, 1.2, "Q2  dual-domain: intra-species vs whole-SC sharing breadth\n-> type2 subtypes (3 groups by cross-species depth)", "#FFF9E6"); arr(5,8.8,5,8.4)
    box(1, 7.1, 8, 1.2, "Q3  clinical association of type2 subtypes\ncarrier-isolate med% + functional preference of cargo", "#FFF9E6"); arr(5,7.1,5,6.7)
    box(1, 5.4, 8, 1.2, "Q4  clinical AC identification\nspecies-mix-adjusted observed vs expected clinical fraction", "#FFF9E6"); arr(5,5.4,5,5.0)
    box(1, 3.7, 8, 1.2, "Q5  RF clinical candidate gene location\n(type2 ACs enriched, but core chromosomes contribute equally)", "#FFF9E6"); arr(5,3.7,5,3.3)
    box(1, 2.0, 8, 1.2, "Q6  BMU02483 case study\nnear-intact clinical type2 ACs acquired (9.2\u00d7 RF enrichment)", "#FFF9E6"); arr(5,2.0,5,1.6)
    box(1, 0.3, 8, 1.2, "Conclusions: kmer types orthogonal to gene-transfer classes; clinical ACs exist (falci + solani);\ntype2 ACs carry clinical-associated genes, but core chromosomes contribute ~45% of the signal", "#E8F4E1")
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_suppl_method.pdf", dpi=300)
    plt.close(fig); print("suppl flowchart saved")


def main():
    mat, v3, dd, sub, gt, ca, sp_map, med = load()
    f1(mat, v3, dd, gt, sub)
    f2(mat, v3, dd, sub, gt, ca, sp_map, med)
    f3(mat, v3, ca)
    suppl_flowchart()
    print("all figures done")


if __name__ == "__main__":
    main()
