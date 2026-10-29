#!/usr/bin/env python3
"""Figure 3 F/H candidate composition, reissued on the current data with the original style.

Data: KEGG level-2 and Pfam GO-MF majority labels per orthogroup (02.annotation/hog_function_labels.tsv)
plus the EffectorP-3 class from the RF table. Drawing replicates the original Figure 3 F/H panels:
tab20 colour map shared across the two species (one colour per category, ordered by total frequency),
grey 'Other', full pies for KEGG/Pfam, the fixed four-colour effector donut, bold outside percentage
labels, a two-column legend, and italic species names in the title. counterclock=False as in the
original.

PyMuPDF is imported after matplotlib on purpose: loading it first makes the bundled C++ runtime
shadow the one matplotlib links against and the import fails with GLIBCXX_3.4.29.

Writes 09.rf_clinical_genes/pav_rerun/panels/rfcomp.{kegg,pfam,effector}.{falciformis,solani}.pdf
        and rfcomp.combined.pdf (2 x 3, 594 x 557 pt)
"""

import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.offsetbox import HPacker, TextArea, AnchoredOffsetbox

import fitz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun/panels")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
LABELS = os.path.join(ROOT, "02.annotation/hog_function_labels.tsv")
STATS = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun/rf_pav_hog_stats.tsv")

SPECIES = [("falciformis", "F. falciformis"), ("solani", "F. solani")]
AXES = [("KEGG_level2", "KEGG"), ("go_mf_class", "Pfam")]
TOPN = 8
PANEL_W, PANEL_H = 2.6, 3.7

FAMILY = "Nimbus Sans"
FS_TITLE, FS_PCT, FS_LEGEND, FS_CENTRE = 12, 10, 7.5, 8
LABEL_R_SIDE, LABEL_R_MID, VIEW_R, PCT_MIN = 1.07, 1.17, 1.45, 6.0
COLOR_OTHERS = "#B0B0B0"
TAB20 = plt.get_cmap("tab20").colors
EFF_COLORS = {"non-effector": "#e8e8e8", "cytoplasmic": "#dc5243",
              "apoplastic": "#4f72a8", "both": "#8c56a2"}

# colours observed in the hand-assembled main Figure 3 F/H panels (RColorBrewer-Paired family,
# hand-tuned at assembly time); the old classifier's labels are translated back to these names
OLD_KEGG_COLORS = {
    "Carbohydrate metabolism": "#3177aa", "Amino acid metabolism": "#f47d34",
    "Lipid metabolism": "#9fbddf", "Xenobiotics biodegradation and metabolism": "#fab17c",
    "Signal transduction": "#866aae", "Energy metabolism": "#8fcc8e",
    "Glycan biosynthesis and metabolism": "#f59193", "Transport and catabolism": "#cb3536",
    "Metabolism of cofactors and vitamins": "#479e47", "Nucleotide metabolism": "#c9aed6",
    "Metabolism of other amino acids": "#b3a2c7",
    "Metabolism of terpenoids and polyketides": "#d8d0e8",
    "Posttranslational modification, protein turnover, chaperones": "#e8b98c",
    "Biosynthesis of other secondary metabolites": "#f7edb0",
    "Folding, sorting and degradation": "#efc9e3", "Translation": "#cde3cd",
    "Transcription": "#d8c8e0", "Membrane transport": "#f6f2e4",
    "Cell growth and death": "#a8d8a8", "Drug resistance: antineoplastic": "#c8b0b0",
}
OLD_PFAM_COLORS = {
    "Translocase / transporter": "#3177aa",
    "Oxidoreductase": "#f47d34",
    "Hydrolase": "#479e47",
    "Binding: protein": "#9fbddf",
    "Transferase": "#f59193",
    "No assignment": "#b0b0b0",
    "Catalytic activity, unspecified": "#8fcc8e",
    "Transcription regulator activity": "#866aae",
    "HET domain": "#cb3536",
    "Protein kinase": "#c9aed6",
    "WD40/TPR repeat": "#b3a2c7",
    "Ankyrin repeat": "#d8d0e8",
    "Peptidase": "#e8b98c",
    "Methyltransferase": "#f7edb0",
}


def resolve_color(label, old_map, fallback_map, fallback_order):
    """translate a new-classifier label back to an old-classifier name and its colour;
    unknown categories continue the tab20 cycle in their own frequency order"""
    base = re.sub(r"^\d{5} ", "", str(label))
    base = re.sub(r" activity \(EC \d\)$", "", base)
    key = base.casefold()
    for old, col in old_map.items():
        ok = old.casefold()
        if key == ok or ok in key or key in ok:
            return col
    if key in fallback_map:
        return fallback_map[key]
    col = TAB20[len(fallback_map) % len(TAB20)]
    fallback_map[key] = col
    return col

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
    "Translocase / transporter": "Transporter",
    "Oxidoreductase activity (EC 1)": "Oxidoreductase",
    "Hydrolase activity (EC 3)": "Hydrolase",
    "Transferase activity (EC 2)": "Transferase",
    "Posttranslational modification, protein turnover, chaperones": "PTM / turnover / chap.",
    "Secondary metabolites biosynthesis, transport and catabolism": "Secondary metab.",
    "Biosynthesis of other secondary metabolites": "Other sec. met.",
    "Catalytic activity, unspecified": "Catalytic, unspec.",
    "Transcription regulator activity": "Transcr. regulator",
    "Binding: protein": "Protein binding",
    "Binding: small molecule": "Small-mol. binding",
    "DNA-binding activity": "DNA binding",
    "Translocase / transporter activity (EC 7)": "Transporter",
}


def shorten(lab):
    lab = re.sub(r"^\d{5} ", "", str(lab))
    lab = re.sub(r" activity \(EC \d\)$", "", lab)
    return SHORT.get(lab, lab)


def build_color_map(dfs, col):
    """one colour per category, shared across the two species' panels, in
    descending order of total frequency (the original figure's scheme)"""
    from collections import Counter
    cats = []
    for df in dfs:
        cats.extend(df[col].dropna().tolist())
    order = [c for c, _ in Counter(cats).most_common()]
    return {c: TAB20[i % len(TAB20)] for i, c in enumerate(order)}


def mixed_title(ax, prefix, species):
    """'KEGG - ' regular + species italic, centred above the pie."""
    t1 = TextArea(prefix, textprops=dict(size=FS_TITLE, fontfamily=FAMILY))
    t2 = TextArea(species, textprops=dict(size=FS_TITLE, fontfamily=FAMILY, style="italic"))
    box = HPacker(children=[t1, t2], align="baseline", pad=0, sep=0)
    ax.add_artist(AnchoredOffsetbox(
        loc="lower center", child=box, pad=0, borderpad=0, frameon=False,
        bbox_to_anchor=(0.5, 1.01), bbox_transform=ax.transAxes))


def outside_pcts(ax, wedges, values):
    """Percentage labels outside the pie, bold, aligned by angle (no leader lines)."""
    total = float(sum(values))
    for w, v in zip(wedges, values):
        pct = v / total * 100
        if pct < PCT_MIN:
            continue
        ang = np.deg2rad((w.theta1 + w.theta2) / 2)
        c, s = np.cos(ang), np.sin(ang)
        ha = "left" if c > 0.15 else ("right" if c < -0.15 else "center")
        r = LABEL_R_SIDE if ha != "center" else LABEL_R_MID
        ax.text(r * c, r * s, f"{pct:.1f}%", ha=ha, va="center",
                fontsize=FS_PCT, fontweight="bold", color="black")


def pie_panel(ax, shares, colours, prefix, species):
    vals = list(shares.values)
    labels = list(shares.index)
    wedges, _ = ax.pie(vals, colors=colours, startangle=90, counterclock=False,
                       wedgeprops={"edgecolor": "white", "linewidth": 0.8})
    outside_pcts(ax, wedges, vals)
    mixed_title(ax, prefix, species)
    ax.legend([Patch(facecolor=c) for c in colours], [shorten(l) for l in labels],
              loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2,
              fontsize=FS_LEGEND, frameon=False, handlelength=1.0,
              handleheight=0.9, columnspacing=1.0, labelspacing=0.28,
              borderaxespad=0.0)


def effector_donut(ax, shares, species):
    order = ["non-effector", "cytoplasmic", "apoplastic", "both"]
    vals = [float(shares.get(c, 0.0)) for c in order]
    colours = [EFF_COLORS[c] for c in order]
    total = sum(vals)
    wedges, _ = ax.pie(vals, colors=colours, startangle=90, counterclock=False,
                       wedgeprops={"edgecolor": "white", "linewidth": 0.8, "width": 0.34})
    outside_pcts(ax, wedges, vals)
    ax.text(0, 0.10, f"{100 - vals[0]:.1f}%", ha="center", va="center",
            fontsize=FS_TITLE - 4, fontweight="bold")
    ax.text(0, -0.16, "effector", ha="center", va="center",
            fontsize=FS_CENTRE, fontweight="bold")
    mixed_title(ax, "Effector - ", species)
    ax.legend([Patch(facecolor=colours[i]) for i in range(len(order))],
              [f"{c} ({vals[i]:.1f}%)" for i, c in enumerate(order)],
              loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2,
              fontsize=FS_LEGEND, frameon=False, handlelength=1.0,
              handleheight=0.9, columnspacing=1.0, labelspacing=0.28,
              borderaxespad=0.0)


def panel(shares, colours, prefix, species, donut=False):
    """one standalone 2.6 x 3.7 in panel; the pie fills the same ~59% of the cell
    width as in the original figure, with the fonts scaled in proportion"""
    fig = plt.figure(figsize=(PANEL_W, PANEL_H))
    frac_w = (110 / 72) / PANEL_W
    frac_h = (110 / 72) / PANEL_H
    ax = fig.add_axes([(1 - frac_w) / 2, 0.30, frac_w, frac_h])
    ax.set_xlim(-VIEW_R, VIEW_R)
    ax.set_ylim(-VIEW_R, VIEW_R)
    ax.set_aspect("equal")
    ax.axis("off")
    if donut:
        effector_donut(ax, shares, species)
    else:
        pie_panel(ax, shares, colours, prefix, species)
    return fig


def main():
    os.makedirs(OUT, exist_ok=True)
    lab = pd.read_csv(LABELS, sep="\t").set_index("HOG_ID")
    st = pd.read_csv(STATS, sep="\t", usecols=["HOG_ID", "EffectorClass"])
    st = st.drop_duplicates("HOG_ID").set_index("HOG_ID")
    lab = lab.join(st)

    subs = {}
    for sp, _ in SPECIES:
        hogs = pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"), sep="\t").HOG_ID
        subs[sp] = lab.reindex(hogs).dropna(subset=["EffectorClass"], how="all")

    # shared fallback maps so the same unknown category gets the same colour in both species
    kegg_fb, pfam_fb = {}, {}

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [FAMILY, "DejaVu Sans"],
        "font.size": FS_LEGEND,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
    })

    written = []
    for sp, sp_label in SPECIES:
        sub = subs[sp]
        for col, axis_name, old_map, fb in (("KEGG_level2", "KEGG", OLD_KEGG_COLORS, kegg_fb),
                                            ("go_mf_class", "Pfam", OLD_PFAM_COLORS, pfam_fb)):
            vals = sub[col].dropna()
            top = vals.value_counts().head(TOPN)
            shares = (top / len(vals) * 100).round(1)
            rest = 100 - shares.sum()
            if rest > 0.05:
                shares.loc["Other"] = round(rest, 1)
            colours = [COLOR_OTHERS if c == "Other"
                       else resolve_color(c, old_map, fb, None) for c in shares.index]
            fig = panel(shares, colours, f"{axis_name} - ", sp_label)
            path = os.path.join(OUT, f"rfcomp.{axis_name.lower()}.{sp}.pdf")
            fig.savefig(path)
            plt.close(fig)
            written.append(path)
            print(f"{os.path.basename(path):44s} n = {len(vals):3d}  "
                  f"top = {shares.index[0]} {shares.iloc[0]:.1f}%")

        ec = sub.EffectorClass.dropna()
        counts = ec.value_counts()
        shares = (counts / len(ec) * 100).round(1)
        order = [k for k in ("non-effector", "cytoplasmic", "apoplastic", "both")
                 if k in shares.index]
        shares = shares.reindex(order)
        fig = panel(shares, None, "Effector - ", sp_label, donut=True)
        path = os.path.join(OUT, f"rfcomp.effector.{sp}.pdf")
        fig.savefig(path)
        plt.close(fig)
        written.append(path)
        print(f"{os.path.basename(path):44s} n = {len(ec):3d}  effector = {100 - shares.get('non-effector', 0):.1f}%")

    # 2 x 3 block, row 1 = panel F (falciformis) row 2 = panel H (solani)
    order = [f"rfcomp.{a.lower()}.{s}.pdf" for s in ("falciformis", "solani")
             for a in ("kegg", "pfam", "effector")]
    docs = [fitz.open(os.path.join(OUT, f)) for f in order]
    W, H, PAD = PANEL_W * 72, PANEL_H * 72, 8
    page_w, page_h = W * 3 + PAD * 4, H * 2 + PAD * 3
    out = fitz.open()
    pg = out.new_page(width=page_w, height=page_h)
    for i, d in enumerate(docs):
        r, c = divmod(i, 3)
        x = PAD + c * (W + PAD)
        y = PAD + r * (H + PAD)
        pg.show_pdf_page(fitz.Rect(x, y, x + W, y + H), d, 0)
        d.close()
    combined = os.path.join(OUT, "rfcomp.combined.pdf")
    out.save(combined, deflate=True, garbage=4)
    out.close()
    print(f"\nwrote {os.path.relpath(combined, ROOT)} "
          f"({os.path.getsize(combined)/1024:.1f} KB, {page_w:.0f}x{page_h:.0f}pt)")


if __name__ == "__main__":
    main()
