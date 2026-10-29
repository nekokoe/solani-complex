#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kmer AC distribution analysis — figures (Fig 1-2 for Q1-Q6 story).

Fig 1 (figure_kmer_flowchart.pdf): method overview flowchart.
Fig 2 (figure_kmer_q1_dual_scatter.pdf): dual-domain scatter (kmer
    intra-species vs inter-species sharing breadth) with gene-transfer
    class overlay, demonstrating the independence of the two dimensions.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def fig1_flowchart(outdir):
    """Method abstract flowchart for the kmer AC distribution analysis."""
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 8, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(8.2, 9.8))
    ax.set_xlim(0, 10); ax.set_ylim(0, 14.5)
    ax.axis("off")

    def box(x, y, w, h, t, fc, fs=7.5, lw=1.0):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.10",
                           linewidth=lw, edgecolor="#444", facecolor=fc, zorder=2)
        ax.add_patch(b)
        ax.text(x + w/2, y + h/2, t, ha="center", va="center", fontsize=fs, zorder=3, linespacing=1.3)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                                     linewidth=1.0, color="#555", zorder=1))

    CI, CP, CR = "#DCE9F7", "#FFF5CC", "#E4F2E1"

    ax.text(5, 14.1, "Kmer-content AC distribution analysis",
            ha="center", va="center", fontsize=12, fontweight="bold")

    # row 1: inputs
    box(0.3, 12.3, 4.5, 1.15, "Inputs\n78 AC kmer profiles x 211 isolates\n(common_table series)", CI)
    box(5.0, 12.3, 4.7, 1.15, "Annotations\nspecies map (211), Med/Env (133/73),\nAC gene sets, Pfam/KEGG/CAZy", CI)
    arrow(4.8, 12.85, 4.95, 12.85)

    # row 2: Q1
    box(0.3, 10.3, 9.4, 1.4, "Q1. Kmer type vs gene-transfer type\nspecies-coherence threshold 0.4 -> 29 type1_stable / 49 type2_mosaic\n"
        + "cross-test with gene-transfer block/scattered class: 62% agreement, "
        + "\u03ba=0.09, independent dimensions", CP, fs=7.0)
    arrow(5, 12.3, 5, 11.75)

    # row 3: Q2
    box(0.3, 8.2, 9.4, 1.6, "Q2. Dual-domain quantifier + type2 subtypes\n"
        + "intra-species vs whole-SC sharing breadth; 47 species-specific, 16 SC-wide, 15 rare\n"
        + "within type2_SC_wide: 3 subtypes (deep-SC / broad sol-falc / pair-restricted)", CP, fs=7.0)
    arrow(5, 10.3, 5, 9.85)

    # row 4: Q3+Q4
    box(0.3, 5.8, 4.5, 1.9, "Q3. Clinical association per subtype\n"
        + "subtype3 (pair-restricted) = 80.2% med\n"
        + "vs 63% bg (OR 2.4, P<0.001)\n"
        + "cargo: transporters + HET domains", CP, fs=6.9)
    box(5.0, 5.8, 4.7, 1.9, "Q4. Clinical ACs (species-stratified)\n"
        + "keratoplastica excluded (small n);\n"
        + "solani/falc: BMU03267_0019 (+29pp),\n"
        + "BMU12683_0014 etc (FDR<0.05)", CP, fs=6.9)
    arrow(5, 8.2, 5, 7.75)

    # row 5: Q5+Q6
    box(0.3, 3.3, 4.5, 2.0, "Q5. RF candidate HOG location\n"
        + "clinical-biased RF: type2 44.9%, core 49%\n"
        + "(OR 1.63 vs non-clinical, P=0.015)\n"
        + "type1 ACs carry only 6-8%", CP, fs=6.9)
    box(5.0, 3.3, 4.7, 2.0, "Q6. BMU02483 case study\n"
        + "carries 28 ACs (24 type2); top 3 near-intact\n"
        + "819 shared HOGs; 85 RF candidates (5.9x)\n"
        + "supports AC-mediated clinical adaptation", CP, fs=6.9)
    arrow(5, 5.8, 5, 5.35)

    # bottom: conclusion
    box(0.3, 1.2, 9.4, 1.5, "Key conclusion\n"
        + "Type2 ACs are distinct from gene-transfer classes; clinical association "
        + "inversely related to cross-species breadth;\n"
        + "clinical ACs exist (beyond burden effect); "
        + "RF candidates split ~50/50 type2-AC / core chromosomes", CR, fs=7.2)
    arrow(5, 3.3, 5, 2.75)

    fig.savefig(f"{outdir}/figure_kmer_flowchart.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {outdir}/figure_kmer_flowchart.pdf")


def fig2_q1_scatter(outdir):
    """Dual-domain scatter: intra vs inter kmer breadth, with gene-transfer overlay."""
    dual = pd.read_csv("06.accessory-chrs/tgs_assemblies/kmer_type.dual_domain.tsv", sep="\t")
    gt = pd.read_csv("06.accessory-chrs/gene_transfer_211/suppl_table.ac_transfer_profile.tsv", sep="\t")
    m = dual.merge(gt[["AC", "AC_class", "genes_per_partner"]], on="AC")

    m["gt_class"] = m["AC_class"].map({"Type1_block": "gene-type1 (block)", "Type2_scattered": "gene-type2 (scattered)"})
    m["kmer_class"] = m["dominant_mode"].map({"type1_species_specific": "kmer-type1", "type2_SC_wide": "kmer-type2",
                                               "rare_isolate_specific": "rare"})

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 8, "pdf.fonttype": 42})
    fig, ax = plt.subplots(figsize=(5.8, 5.2))
    colors = {"kmer-type1": "#4C72B0", "kmer-type2": "#DD8452", "rare": "#AAAAAA"}
    markers = {"gene-type1 (block)": "s", "gene-type2 (scattered)": "o"}
    for _, r in m.iterrows():
        kc = r["kmer_class"]; gc = r["gt_class"]
        ax.scatter(r["intra_breadth_0.4"], r["inter_breadth_0.4"],
                   c=colors[kc], marker=markers[gc], s=28, edgecolor="black", linewidths=0.3,
                   alpha=0.8, zorder=2)

    # label key ACs
    key = {"BMU01933_0021": (2, 6), "BMU02474_0015": (6, 2), "BMU03267_0019": (8, -8),
           "BMU12683_0014": (4, -10), "SX13_10": (6, -4), "BMU01933_0007": (2, -10),
           "BMU02028_0014": (2, -4), "BMU12683_0016": (4, 4)}
    for ac, (dx, dy) in key.items():
        r = m[m.AC == ac].iloc[0]
        ax.annotate(ac.split("_", 1)[-1] if "_" in ac else ac,
                    (r["intra_breadth_0.4"], r["inter_breadth_0.4"]),
                    textcoords="offset points", xytext=(dx, dy), fontsize=5.2,
                    arrowprops=dict(arrowstyle="-", color="#777", lw=0.5))
    ax.axvline(0.3, color="gray", linestyle="--", linewidth=0.7)
    ax.axhline(0.3, color="gray", linestyle="--", linewidth=0.7)
    ax.text(0.31, 0.90, "\u2190 intra threshold", fontsize=6, va="bottom", color="gray")
    ax.text(0.85, 0.31, "inter threshold \u2193", fontsize=6, ha="right", color="gray")
    ax.set_xlabel("Intra-species sharing breadth (kmer \u2265 0.4)", fontsize=8)
    ax.set_ylabel("Inter-species sharing breadth (kmer \u2265 0.4)", fontsize=8)
    ax.set_title("ACs in kmer-type vs gene-transfer-type space", fontsize=9, pad=3)
    ax.tick_params(labelsize=7); ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.3, linewidth=0.5)
    # custom legend: manually
    from matplotlib.lines import Line2D
    leg = [Line2D([], [], marker="o", color="w", markerfacecolor="#DD8452", markersize=8, markeredgecolor="black", markeredgewidth=0.3, label="kmer-type2 (mosaic)"),
           Line2D([], [], marker="o", color="w", markerfacecolor="#4C72B0", markersize=8, markeredgecolor="black", markeredgewidth=0.3, label="kmer-type1 (stable)"),
           Line2D([], [], marker="s", color="w", markerfacecolor="#777", markersize=8, markeredgecolor="black", markeredgewidth=0.3, label="gene-type1 (block)"),
           Line2D([], [], marker="o", color="w", markerfacecolor="#777", markersize=8, markeredgecolor="black", markeredgewidth=0.3, label="gene-type2 (scattered)")]
    ax.legend(handles=leg, fontsize=6.5, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{outdir}/figure_kmer_q1_dual_scatter.pdf", dpi=300)
    plt.close(fig)
    print(f"saved {outdir}/figure_kmer_q1_dual_scatter.pdf")


if __name__ == "__main__":
    out = "06.accessory-chrs/gene_transfer_211"
    fig1_flowchart(out)
    fig2_q1_scatter(out)
