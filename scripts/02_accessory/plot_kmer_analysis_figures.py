#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Panels for the k-mer based accessory-chromosome (AC) distribution analysis.

The six panels used to be called Q1-Q6, after the six questions the analysis was
originally framed around. Those codes say nothing about the content, so the
functions and the output files now carry descriptive names; the old codes are
kept in this table only, so that notes, commit messages and the lab notebook
remain searchable:

    Q1  -> plot_class_orthogonality        k-mer block/mosaic class vs gene-transfer class
    Q2  -> plot_sharing_continuum          within-species vs whole-complex sharing; mosaic subtypes
    Q3  -> plot_subtype_clinical_bias      clinical bias of the mosaic subtypes
    Q4  -> plot_clinical_ac_identification observed vs species-mix-adjusted clinical fraction
    Q5  -> plot_rf_candidate_location      where the RF candidate HOGs sit (mosaic AC / block AC / core)
    Q6  -> plot_bmu02483_sharing           BMU02483 carrier profile against the clinical ACs

Outputs (all in 06.accessory-chrs/kmer_analysis/):
    figure_kmer_flowchart.pdf
    figure_kmer_class_orthogonality.pdf
    figure_kmer_sharing_continuum.pdf
    figure_kmer_subtype_clinical_bias.pdf
    figure_kmer_clinical_ac_identification.pdf
    figure_kmer_rf_candidate_location.pdf
    figure_kmer_bmu02483_sharing.pdf
    rf_candidate_location_by_ac_class.tsv   (the numbers behind the RF-location panel)
"""

import os
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy import stats

D = "06.accessory-chrs/kmer_analysis"
SH = "06.accessory-chrs/tgs_assemblies"
OUT = D
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42})
C_T2 = "#C44E52"
C_T1 = "#4C72B0"
C_RARE = "#B0B0B0"


def load():
    mat = pd.read_csv(f"{D}/kmer_matrix.ACx211.tsv", sep="\t", index_col=0)
    kt = pd.read_csv(f"{D}/kmer_type.classification.v3.tsv", sep="\t")
    dd = pd.read_csv(f"{D}/kmer_type.dual_domain.tsv", sep="\t")
    sub = pd.read_csv(f"{D}/kmer_type2.subtypes.tsv", sep="\t")
    gt = pd.read_csv("06.accessory-chrs/gene_transfer_211/suppl_table.ac_transfer_profile.tsv",
                     sep="\t")
    ca = pd.read_csv(f"{D}/kmer_clinical_ac.assessment.tsv", sep="\t")
    sp_map, med = {}, {}
    for line in open(f"{SH}/species.211.map.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            sp_map[p[0]] = p[1]
    for line in open(f"{SH}/bmu_ref_groups.211.txt"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0] != "BMU_ID":
            med[p[0]] = p[1]
    return mat, kt, dd, sub, gt, ca, sp_map, med


def flowchart():
    fig, ax = plt.subplots(figsize=(7.6, 9.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 13.5)
    ax.axis("off")

    def box(x, y, w, h, text, fc, fs=7.4):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                           linewidth=1.0, edgecolor="#333333", facecolor=fc, zorder=2)
        ax.add_patch(b)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
                zorder=3, linespacing=1.4)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=12, linewidth=1.1,
                                     color="#555555", zorder=1))

    C_IN, C_PR, C_RES = "#DCE9F7", "#FFF2CC", "#E4F2E1"
    ax.text(5, 13.2, "k-mer-based distribution of accessory chromosomes\nacross the Fusarium solani species complex",
            ha="center", fontsize=10.5, fontweight="bold")
    box(1.0, 11.4, 8.0, 1.2,
        "Input: 78 accessory chromosomes (13 TGS) x 211 isolates\nshared-kmer proportion matrix (kmer length 19, kmer-db)", C_IN)
    arrow(5, 11.4, 5, 10.9)
    box(1.0, 9.5, 8.0, 1.35,
        "Step 1: classify ACs by k-mer sharing species coherence\nType I = sharing confined to one species (clade/species-specific)\nType II = cross-species mosaic", C_PR)
    arrow(5, 9.5, 5, 9.0)
    box(1.0, 7.7, 8.0, 1.35,
        "Step 2: dual-domain qualifier\nintra-species breadth vs whole-SC breadth per AC\n-> type2_SC_wide split into 3 subtypes by cross-species depth", C_PR)
    arrow(5, 7.7, 5, 7.2)
    box(1.0, 5.9, 8.0, 1.35,
        "Step 3: clinical association\ncarrier isolates (kmer >= 0.4) of each subtype/AC\nvs 63% med background + species-mix adjustment", C_PR)
    arrow(5, 5.9, 5, 5.4)
    box(1.0, 4.1, 8.0, 1.35,
        "Step 4: do truly clinical ACs exist?\ncontrol for isolate AC-burden (32 vs 10) and species mix", C_PR)
    arrow(5, 4.1, 5, 3.6)
    box(1.0, 2.3, 8.0, 1.35,
        "Step 5: integrate with gene transfer and the pan-genome RF screen\nwhere do the clinical-associated (RF) genes sit?\ncase study: BMU02483 accessory-chromosome acquisition", C_PR)
    arrow(5, 2.3, 5, 1.8)
    box(1.0, 0.4, 8.0, 1.3,
        "Output: AC distribution map + clinical ACs + gene-function bias\n- kmer types are orthogonal to gene-transfer classes\n- clinical type2 ACs carry transporters/HET; RF genes ~45% core\n- BMU02483 carries near-intact clinical mosaic ACs", C_RES)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_flowchart.pdf", dpi=300)
    plt.close(fig)
    print("flowchart done")


def plot_class_orthogonality(mat, kt, gt, dd):
    m = kt[["AC", "kmer_type"]].merge(gt[["AC", "AC_class", "genes_per_partner"]], on="AC")
    m = m.merge(dd[["AC", "specificity"]], on="AC")
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    X_B = 2.0  # block/scattered boundary (genes per partner)
    Y_B = 0.2  # species-specific/mosaic boundary (specificity index)
    for typ, c in [("type1_stable", C_T1), ("type2_mosaic", C_T2)]:
        for gcls, mk in [("Type2_scattered", "o"), ("Type1_block", "^")]:
            s = m[(m.kmer_type == typ) & (m.AC_class == gcls)]
            ax.scatter(s.genes_per_partner, s.specificity, s=55, c=c, marker=mk,
                       alpha=0.75, edgecolors="black", linewidths=0.4, zorder=3)
    ax.axvline(X_B, color="black", linestyle="--", linewidth=0.9)
    ax.axhline(Y_B, color="black", linestyle="--", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 60)
    ax.set_ylim(-0.45, 0.85)
    ax.text(0.12, 0.78, "kmer Type I\nspecies-specific\n+ scattered transfer\n(the solani AC family)",
            fontsize=6, va="top")
    ax.text(30, 0.78, "kmer Type I +\nblock transfer", fontsize=6, ha="right", va="top")
    ax.text(0.12, -0.4, "kmer Type II mosaic\n+ scattered transfer\n(classic type II)",
            fontsize=6, va="bottom")
    ax.text(30, -0.4, "kmer Type II +\nblock transfer", fontsize=6, ha="right", va="bottom")
    ax.set_xlabel("Gene-transfer concentration (genes per partner; block right / scattered left)",
                  fontsize=7.5)
    ax.set_ylabel("Kmer species-specificity index (intra \u2212 inter sharing)", fontsize=7.5)
    ax.set_title("k-mer type vs gene-transfer class\nquadrant view - the two classifications are orthogonal (kappa = 0.09)",
                 fontsize=8.8, pad=3)
    ax.scatter([], [], s=55, c=C_T1, marker="o", label="kmer Type I (species-specific)")
    ax.scatter([], [], s=55, c=C_T2, marker="o", label="kmer Type II (mosaic)")
    ax.scatter([], [], s=40, c="grey", marker="o", label="scattered-transfer (gene class 2)")
    ax.scatter([], [], s=40, c="grey", marker="^", label="block-transfer (gene class 1)")
    ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_class_orthogonality.pdf", dpi=300)
    plt.close(fig)
    print("class vs gene-transfer done")


def plot_sharing_continuum(mat, kt, dd, sub, sp_map):
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    colorm = {"type1_species_specific": C_T1, "type2_SC_wide": C_T2,
              "rare_isolate_specific": C_RARE}
    for dm, c in colorm.items():
        s = dd[dd.dominant_mode == dm]
        ax.scatter(s["intra_breadth_0.4"], s["inter_breadth_0.4"], s=45, c=c,
                   alpha=0.8, edgecolors="black", linewidths=0.3, label=dm, zorder=3)
    for ac in ["BMU12683_0013", "SX13_20", "BMU03267_0024"]:
        r = dd[dd.AC == ac]
        ax.annotate(ac.replace("BMU", "B"), (r["intra_breadth_0.4"].iloc[0],
                    r["inter_breadth_0.4"].iloc[0]), fontsize=5.5,
                    textcoords="offset points", xytext=(3, 3))
    ax.plot([0, 1], [0.3, 0.3], color="black", linestyle=":", linewidth=0.8)
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=0.6)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("Intra-species sharing breadth (frac of own-species isolates \u2265 0.4)", fontsize=7.5)
    ax.set_ylabel("Cross-species sharing breadth (frac of other-species isolates \u2265 0.4)", fontsize=7.5)
    ax.set_title("within-species vs whole-complex sharing (78 ACs)\ntype1/type2 is a continuum, not a dichotomy", fontsize=8.8, pad=3)
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)

    ax = axes[1]
    spec = ["solani", "falciformis", "keratoplastica", "other"]
    labels = {"subtype1": "subtype1\ndeep SC-wide", "subtype2": "subtype2\nbroad sol-falc",
              "subtype3": "subtype3\npair-restricted"}
    x = np.arange(3)
    w = 0.2
    cols = {"solani": "#4C72B0", "falciformis": "#DD8452", "keratoplastica": "#55A868",
            "other": "#B0B0B0"}
    for i, sp in enumerate(spec):
        vals = [sub[sub.subtype == s].iloc[0][f"br_{sp}"] if (sub[sub.subtype == s]["br_" + sp] if False else True) else 0
                for s in ["subtype1", "subtype2", "subtype3"]]
        vals = []
        for s in ["subtype1", "subtype2", "subtype3"]:
            r = sub[sub.subtype == s]
            vals.append(r[f"br_{sp}"].mean() if f"br_{sp}" in r.columns else 0)
        ax.bar(x + (i - 1.5) * w, vals, w, color=cols[sp], label=sp)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[s] for s in ["subtype1", "subtype2", "subtype3"]], fontsize=6)
    ax.set_ylabel("Sharing breadth per species", fontsize=7)
    ax.set_title("Type-II subtypes\ncross-species depth", fontsize=8.8, pad=3)
    ax.legend(fontsize=5.5, frameon=False)
    ax.tick_params(labelsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_sharing_continuum.pdf", dpi=300)
    plt.close(fig)
    print("sharing continuum + subtypes done")


def plot_subtype_clinical_bias(mat, sub, med, sp_map):
    iso211 = [l.strip() for l in open("02.annotation/sample.solaniSC.211.list") if l.strip()]
    n_med_bg = sum(1 for i in iso211 if med.get(i) == "med")
    rows = []
    for st in ["subtype1", "subtype2", "subtype3"]:
        acs = sub[sub.subtype == st].AC
        cars = set()
        for ac in acs:
            v = mat.loc[ac]
            cars |= {c for c in v.index if v[c] >= 0.4}
        cars = [c for c in iso211 if c in cars]
        n_med = sum(1 for c in cars if med.get(c) == "med")
        or_, p = stats.fisher_exact([[n_med, len(cars) - n_med],
                                     [n_med_bg, len(iso211) - n_med_bg]])
        rows.append({"subtype": st, "n": len(cars), "med_frac": n_med / len(cars),
                     "p": p})
    r = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    x = np.arange(3)
    ax.bar(x, r.med_frac * 100, 0.55, color=[C_T2, "#E8836B", "#F5B9AD"])
    ax.axhline(n_med_bg / len(iso211) * 100, color="black", linestyle="--", linewidth=1.0)
    ax.text(2.35, n_med_bg / len(iso211) * 100 + 0.8, f"genome background\n{n_med_bg/len(iso211)*100:.0f}% med",
            fontsize=6, ha="right", va="bottom")
    for i, (_, rr) in enumerate(r.iterrows()):
        ax.text(x[i], rr.med_frac * 100 + 1, f"{rr.med_frac*100:.1f}%", ha="center", fontsize=8)
        ax.text(x[i], 4, f"n={rr.n}", ha="center", fontsize=6, color="white")
    ax.set_xticks(x)
    ax.set_xticklabels(["subtype1\ndeep SC-wide", "subtype2\nbroad sol-falc",
                        "subtype3\npair-restricted"], fontsize=6.5)
    ax.set_ylabel("Carrier isolates that are clinical (%)", fontsize=7.5)
    ax.set_title("mosaic (type-II) subtypes are clinical-associated\n(most: subtype3, OR = 2.4, P < 0.001)", fontsize=8.8, pad=3)
    ax.set_ylim(0, 95)
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_subtype_clinical_bias.pdf", dpi=300)
    plt.close(fig)
    print("subtype clinical bias done")


def plot_clinical_ac_identification(mat, kt, med, sp_map):
    iso211 = [l.strip() for l in open("02.annotation/sample.solaniSC.211.list") if l.strip()]
    valid = [i for i in iso211 if med.get(i) in ("med", "env")]
    sp_med = {}
    for sp in ["solani", "falciformis", "keratoplastica"]:
        iso = [i for i in valid if sp_map.get(i) == sp]
        sp_med[sp] = sum(1 for i in iso if med[i] == "med") / len(iso)
    t2 = set(kt[kt.kmer_type == "type2_mosaic"].AC)
    rows = []
    for ac in mat.index:
        v = mat.loc[ac]
        cars = [c for c in valid if v[c] >= 0.4]
        if len(cars) < 5:
            continue
        from collections import Counter
        spc = Counter(sp_map.get(c) for c in cars)
        exp = sum(sp_med.get(s, 0.65) * n for s, n in spc.items()) / len(cars)
        n_med = sum(1 for c in cars if med[c] == "med")
        rows.append({"AC": ac, "obs": n_med / len(cars), "exp": exp,
                     "type2": ac in t2})
    a = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    a = a.sort_values("obs")
    y = np.arange(len(a))
    ax.scatter(a.obs, y, s=20, c=[C_T2 if t else C_T1 for t in a.type2],
               alpha=0.8, edgecolors="black", linewidths=0.2, zorder=3, label="type2")
    ax.scatter(a.exp, y, s=6, c="#B0B0B0", zorder=2, label="species-adjusted expectation")
    n_lab = 0
    for i in range(len(a)):
        if a.obs.iloc[i] - a.exp.iloc[i] > 0.2:
            ax.plot([a.exp.iloc[i], a.obs.iloc[i]], [y[i], y[i]], color="black",
                    linewidth=0.5, alpha=0.6)
        if a.obs.iloc[i] - a.exp.iloc[i] > 0.25 and n_lab < 10:
            ax.annotate(a.AC.iloc[i].split("_")[-1], (a.obs.iloc[i] + 0.01, y[i]),
                        fontsize=4.2, va="center")
            n_lab += 1
    ax.set_yticks([])
    ax.set_xlim(0.3, 1.05)
    ax.set_xlabel("Clinical fraction of carrier isolates (kmer \u2265 0.4)", fontsize=7.5)
    ax.set_title("clinical ACs after species-mix adjustment\n(black ties: carriers \u226515 pp more clinical than expected)", fontsize=8.8, pad=3)
    ax.legend(fontsize=6.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_clinical_ac_identification.pdf", dpi=300)
    plt.close(fig)
    print("clinical AC identification done")


HOG_AC = "06.accessory-chrs/hog_to_ac.tsv"
# Candidate set for the RF-location panel. The PAV screen is what the manuscript reports; the
# count-based 09.rf_clinical_genes/rf_genes.hog_stats.tsv is the superseded one
# (it is why 21 orthogroups per species had a non-zero feature importance yet an
# identical clinical / non-clinical presence rate).
RF_STATS = "09.rf_clinical_genes/pav_rerun/rf_pav_hog_stats.tsv"


def load_hog_ac(path=HOG_AC):
    """HOG -> set of accessory chromosomes its genes sit on.

    Read straight from the join table instead of re-typing every contig here:
    N0.gene_vs_HOG.tsv (transcript -> HOG) joined to
    all_solaniSC_isolates.gff3.affliation (transcript -> contig), where the contig
    name IS the AC name. Note the id columns do not line up naively: the column
    called GeneID in N0.gene_vs_HOG.tsv holds transcript ids (-T1 suffix), so it
    must be matched to TranscriptID, not to the affiliation file's GeneID.
    """
    tab = pd.read_csv(path, sep="\t", usecols=["HOG_ID", "AC"])
    return tab.groupby("HOG_ID")["AC"].apply(set).to_dict()


def classify_loc(acs, t2a, t1a):
    """A HOG is placed on the highest-ranking class any of its genes occupies.

    This 'any member' rule is the one that reproduces the published numbers. The
    earlier implementation typed each contig and then did `hog_loc[h] = ct` gene
    by gene, so a HOG with genes on both a mosaic AC and the core was labelled by
    whichever gene happened to be visited last; typing a genuinely mosaic HOG as
    core is what produced 22.0%/13.8% instead of 44.9%/33.3%.
    """
    if acs & t2a:
        return "type2_AC"
    if acs & t1a:
        return "type1_AC"
    return "core_or_no13"


def plot_rf_candidate_location(rf_delta, hog_ac, t2a, t1a):
    clin = [h for h, d in rf_delta.items() if d > 0]
    nclin = [h for h, d in rf_delta.items() if d < 0]
    cats = ["type2_AC", "type1_AC", "core_or_no13"]
    loc = {h: classify_loc(hog_ac.get(h, set()), t2a, t1a)
           for h in list(clin) + list(nclin)}

    rows = []
    for c in cats:
        a = sum(1 for h in clin if loc[h] == c)
        b = sum(1 for h in nclin if loc[h] == c)
        rows.append(dict(location=c, clinical_n=a, clinical_total=len(clin),
                         clinical_pct=a / len(clin) * 100,
                         nonclinical_n=b, nonclinical_total=len(nclin),
                         nonclinical_pct=b / len(nclin) * 100))
    # the headline contrast is mosaic AC versus everything else; computed here so the
    # title can never go stale again
    a = rows[0]["clinical_n"]
    b = rows[0]["nonclinical_n"]
    tab2 = [[a, len(clin) - a], [b, len(nclin) - b]]
    or_ = (a / (len(clin) - a)) / (b / (len(nclin) - b))
    p_ = stats.fisher_exact(tab2)[1]

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    x = np.arange(2)
    w = 0.25
    for i, r in enumerate(rows):
        ax.bar(x[0] + (i - 1) * w, r["clinical_pct"], w, color=[C_T2, C_T1, "#8C8C8C"][i])
        ax.bar(x[1] + (i - 1) * w, r["nonclinical_pct"], w,
               color=[C_T2, C_T1, "#8C8C8C"][i], alpha=0.35)
    ax.set_xticks(x)
    ax.set_xticklabels(["Clinical-biased RF HOGs\n(n = %d)" % len(clin),
                        "Non-clinical-biased RF HOGs\n(n = %d)" % len(nclin)], fontsize=7)
    ax.set_ylabel("% of RF HOGs with 13-TGS members on the class", fontsize=7)
    ax.set_title("where do the clinical-associated genes live?\n"
                 "(type2 ACs enriched in clinical-biased RF HOGs, "
                 "OR = %.2f, P = %.3g)" % (or_, p_), fontsize=8.8, pad=3)
    ax.legend(handles=[plt.Line2D([0], [0], color=C_T2, lw=4, label="type2 AC"),
                       plt.Line2D([0], [0], color=C_T1, lw=4, label="type1 AC"),
                       plt.Line2D([0], [0], color="#8C8C8C", lw=4, label="core")],
              fontsize=6, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.08))
    # headroom must follow the data: core_or_no13 reaches 71.4% on the PAV set
    ax.set_ylim(0, max(max(r["clinical_pct"], r["nonclinical_pct"])
                       for r in rows) * 1.12)
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_rf_candidate_location.pdf", dpi=300)
    plt.close(fig)

    out = pd.DataFrame(rows)
    out["type2_odds_ratio"] = or_
    out["type2_fisher_p"] = p_
    out.to_csv(f"{OUT}/rf_candidate_location_by_ac_class.tsv", sep="\t", index=False)
    print("RF candidate location done · " + " · ".join(
        f"{r['location']} {r['clinical_pct']:.1f}% vs {r['nonclinical_pct']:.1f}%"
        for r in rows) + f" · type2 OR={or_:.2f} P={p_:.3g}")



def plot_bmu02483_sharing(mat, kt, ca):
    v = mat["BMU02483"]
    s = v.sort_values(ascending=False)
    s = s[s >= 0.35]
    t2 = set(kt[kt.kmer_type == "type2_mosaic"].AC)
    clin = set(ca[ca.FDR < 0.05].AC)
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    y = np.arange(len(s))
    colors = []
    for ac in s.index:
        if ac in t2 and ac in clin:
            colors.append(C_T2)
        elif ac in t2:
            colors.append("#E8836B")
        else:
            colors.append(C_T1)
    ax.barh(y, s.values, 0.7, color=colors, edgecolor="black", linewidth=0.2)
    ax.axvline(0.4, color="black", linestyle=":", linewidth=0.8)
    ax.text(0.405, len(s) - 0.4, "carrier threshold 0.4", fontsize=6, color="black")
    for i, ac in enumerate(s.index[:8]):
        ax.text(s.values[i] + 0.01, y[i], ac.replace("BMU", "B"), fontsize=5.5, va="center")
    ax.set_yticks([])
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Kmer sharing of BMU02483 with each AC", fontsize=7.5)
    ax.set_title("BMU02483 carries near-intact clinical mosaic ACs\n(red = clinical mosaic AC; bar length = k-mer sharing with that AC)", fontsize=8.8, pad=3)
    ax.legend(handles=[plt.Line2D([0], [0], color=C_T2, lw=4, label="clinical type2 AC"),
                       plt.Line2D([0], [0], color="#E8836B", lw=4, label="type2 AC"),
                       plt.Line2D([0], [0], color=C_T1, lw=4, label="type1 AC")],
              fontsize=6.5, frameon=False, ncol=3,
              bbox_to_anchor=(0.5, -0.06), loc="upper center")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(f"{OUT}/figure_kmer_bmu02483_sharing.pdf", dpi=300)
    plt.close(fig)
    print("BMU02483 sharing profile done")


def main():
    mat, kt, dd, sub, gt, ca, sp_map, med = load()
    flowchart()
    plot_class_orthogonality(mat, kt, gt, dd)
    plot_sharing_continuum(mat, kt, dd, sub, sp_map)
    plot_subtype_clinical_bias(mat, sub, med, sp_map)
    plot_clinical_ac_identification(mat, kt, med, sp_map)
    # RF candidate location: clinical-bias split from the candidate table,
    # chromosome class from the HOG-to-AC join
    rf_stats = pd.read_csv(RF_STATS, sep="\t")
    rf_delta = rf_stats.groupby("HOG_ID")["Rate_delta"].max()
    hog_ac = load_hog_ac()
    t2a = set(kt[kt.kmer_type == "type2_mosaic"].AC)
    t1a = set(kt[kt.kmer_type == "type1_stable"].AC)
    plot_rf_candidate_location(rf_delta.to_dict(), hog_ac, t2a, t1a)
    # BMU02483 carrier profile
    ca = pd.read_csv(f"{D}/kmer_clinical_ac.assessment.tsv", sep="\t")
    plot_bmu02483_sharing(mat, kt, ca)
    print("all figures done")


if __name__ == "__main__":
    main()
