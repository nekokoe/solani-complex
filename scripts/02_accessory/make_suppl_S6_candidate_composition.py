#!/usr/bin/env python3
"""Supplementary Figure S6 - functional composition of the random-forest candidate orthogroups.

The candidates are shown against the orthogroups of the species they were selected in, on the
three axes the Results subsection discusses: KEGG level-2, Pfam molecular function, and the
EffectorP-3 class.

Both axes are read from 02.annotation/hog_function_labels.tsv, i.e. from the same
one-orthogroup-one-vote labelling the rest of the pan-genome section uses, so this figure cannot
drift from Supplementary Table S20 or from the text. The earlier version of this figure applied a
"first matching keyword over the union of a HOG's domains" rule, which reported 50.0% and 62.8%
transporters; that rule tracks orthogroup size rather than function and is not used here.

Writes 00.manuscript/supplementary/Supplementary-Fig-S6-rf-HOG-composition.{pdf,png}
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "00.manuscript/supplementary")
BASE = "Supplementary-Fig-S6-rf-HOG-composition"

LABELS = os.path.join(ROOT, "02.annotation/hog_function_labels.tsv")
EFFECTOR = os.path.join(ROOT, "08.pan_020_genes/effector_master.hog_level.tsv")
# the candidate EffectorP call is read from the RF table, which carries it for all 655 candidates;
# effector_master covers 439 of them. The two agree on every orthogroup they share (292 effector /
# 147 non-effector), so the candidate bars use the complete column and the background uses
# effector_master across the species.
CAND_EFFECTOR = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun/rf_pav_hog_stats.tsv")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")

FONT = "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf"
FONT_B = "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf"
C_CAND = "#C92626"
C_BG = "#BEBEBE"

SPECIES = [("falciformis", "F. falciformis"), ("solani", "F. solani")]
AXES = [("KEGG_level2", "KEGG level-2"), ("go_mf_class", "Pfam molecular function")]
TOPN = 7


def short(s):
    s = str(s)
    for a, b in ((" activity", ""), ("Metabolism of ", ""), ("metabolism", "metab."),
                 ("Translocase / transporter", "Translocase"),
                 ("Biosynthesis of other secondary metabolites", "Secondary metab.")):
        s = s.replace(a, b)
    return s if len(s) <= 34 else s[:31] + "..."


def shares(series, n_total):
    s = series.value_counts()
    return (s / n_total * 100).round(2)


def main():
    lab = pd.read_csv(LABELS, sep="\t").set_index("HOG_ID")
    eff = (pd.read_csv(EFFECTOR, sep="\t", usecols=["HOG_ID", "DominantClass"])
           .drop_duplicates("HOG_ID").set_index("HOG_ID"))
    lab = lab.join(eff)
    cand_eff = (pd.read_csv(CAND_EFFECTOR, sep="\t", usecols=["HOG_ID", "EffectorClass"])
                .drop_duplicates("HOG_ID").set_index("HOG_ID").EffectorClass)

    cand = {sp: pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"), sep="\t").HOG_ID
            for sp, _ in SPECIES}

    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                         "ytick.major.width": 0.6})
    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.4))
    letters = "abcdef"

    for ri, (sp, sp_label) in enumerate(SPECIES):
        sub = lab[lab.index.isin(
            pd.read_csv(os.path.join(ROOT, "05.stats_pan_by_species_shorter/per_species",
                                     f"{sp}.HOG.interval.tsv"), sep="\t").HOG_ID)]
        in_cand = sub.index.isin(cand[sp])
        for ci, (col, title) in enumerate(AXES + [("DominantClass", "EffectorP-3 class")]):
            ax = axes[ri][ci]
            if col == "DominantClass":
                a = cand_eff.reindex(sub.index[in_cand]).dropna()
                b = sub.loc[~in_cand, col].dropna()
            else:
                a = sub.loc[in_cand, col].dropna()
                b = sub.loc[~in_cand, col].dropna()
            if not len(a) or not len(b):
                ax.axis("off")
                continue
            top = a.value_counts().head(TOPN).index.tolist()
            if col == "DominantClass":
                top = ["cytoplasmic", "apoplastic", "non-effector"]
            pa = (a.value_counts() / len(a) * 100).reindex(top).fillna(0)
            pb = (b.value_counts() / len(b) * 100).reindex(top).fillna(0)
            y = range(len(top))
            ax.barh([i + 0.20 for i in y], pa.values, height=0.38, color=C_CAND,
                    label="candidates")
            ax.barh([i - 0.20 for i in y], pb.values, height=0.38, color=C_BG,
                    label="other orthogroups of the species")

            for i, v in zip(y, pa.values):
                ax.text(v + 0.6, i + 0.20, f"{v:.1f}", va="center", fontsize=6, color=C_CAND)
            ax.set_yticks(list(y))
            ax.set_yticklabels([short(t) for t in top], fontsize=7)
            ax.invert_yaxis()
            # coverage differs between axes (KEGG annotates fewer orthogroups than Pfam, and the
            # effector call covers a subset), so both counts travel with the axis label
            ax.set_xlabel(f"% of orthogroups (n = {len(a):,} / {len(b):,})", fontsize=6.5)
            ax.tick_params(axis="x", labelsize=6.5)
            ax.set_title(f"{letters[ri*3+ci]}  {title} — {sp_label}", fontsize=8.5,
                         loc="left", fontweight="bold", fontfamily="sans-serif")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="x", color="#EAEAEA", lw=0.6)
            ax.set_axisbelow(True)
            if ri == 0 and ci == 0:
                handles, labs = ax.get_legend_handles_labels()
                fig.legend(handles, labs, fontsize=7.5, frameon=False,
                           loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.015))

    fig.tight_layout(h_pad=1.6, w_pad=2.2, rect=(0, 0, 1, 0.96))
    pdf = os.path.join(OUTDIR, BASE + ".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, BASE + ".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(pdf, ROOT)} ({os.path.getsize(pdf)/1024:.1f} KB)")
    print(f"wrote {os.path.relpath(pdf[:-4] + '.png', ROOT)}")


if __name__ == "__main__":
    main()
