#!/usr/bin/env python3
"""Supplementary Figure S41 - effector and resistance loci across clinical and non-clinical isolates.

Four panels: one per species (F. solani, F. falciformis) and one per axis (effector orthogroups,
antifungal-resistance orthogroups). Rows are the orthogroups of that axis whose presence differs
most between clinical and non-clinical isolates (Mann-Whitney on per-genome presence, BH within
the axis and species; capped at the 50 smallest q). Columns are the isolates, clinical first.
Cells give the number of GENE COPIES the isolate carries in that orthogroup, i.e. the gene-level
quantity the manuscript reports; 0 is white and the scale saturates at five copies. The column
strip marks isolation source with the same three classes as Figure 3E and 3G (clinical, non-clinical,
unspecified source); the test that ranks the rows treats clinical against everything else, which is
the binary label the random-forest screen itself uses.

Writes 00.manuscript/supplementary/Supplementary-Fig-S41-effector-resistance-heatmap.{pdf,png}
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "00.manuscript/supplementary")
BASE = "Supplementary-Fig-S41-effector-resistance-heatmap"

HOG_COUNT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
PER_GENOME = os.path.join(ROOT, "09.rf_clinical_genes/rf_genes.per_genome.tsv")
RES_ASSOC = os.path.join(ROOT, "08.pan_020_genes/resistance_genes_clinical_vs_env.tsv")
EFF_ASSOC = os.path.join(ROOT, "08.pan_020_genes/effector_source_assoc.tsv")

SPECIES = [("solani", "F. solani"), ("falciformis", "F. falciformis")]
C_CLIN, C_NON, C_UNK = "#C92626", "#4D8AC1", "#9E9E9E"
SRC_CODE = {"env": 0, "unknown": 1, "med": 2}
SRC_CMAP = matplotlib.colors.ListedColormap([C_NON, C_UNK, C_CLIN])
TOPN = 50
VMAX = 5


def main():
    pg = pd.read_csv(PER_GENOME, sep="\t", usecols=["SampleID", "Source", "Species"])
    counts = pd.read_csv(HOG_COUNT, sep="\t")
    res = pd.read_csv(RES_ASSOC, sep="\t")
    res = res[res.q_bh < 0.05]
    eff = pd.read_csv(EFF_ASSOC, sep="\t")
    eff = eff[eff.q_bh < 0.05]

    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                         "axes.linewidth": 0.6})
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.5))
    letters = "abcd"
    ims = []

    for ri, (sp, sp_label) in enumerate(SPECIES):
        sub = pg[pg.Species == sp]
        order = list(sub[sub.Source == "med"].SampleID) + list(sub[sub.Source != "med"].SampleID)
        for ci, (tab, title) in enumerate(((eff, "Effector orthogroups"),
                                           (res, "Antifungal-resistance orthogroups"))):
            ax = axes[ri][ci]
            t = tab[tab.species == sp].nsmallest(TOPN, "q_bh")
            if not len(t):
                ax.axis("off")
                continue
            wide = (counts[counts.HOG_ID.isin(t.HOG_ID)]
                    .groupby(["HOG_ID", "SampleID"]).GeneCount.sum().unstack()
                    .reindex(columns=order).fillna(0))
            wide = wide.reindex(t.HOG_ID)          # rows ordered by q
            mat = wide.values
            mat = np.where(mat > VMAX, VMAX, mat)
            im = ax.imshow(mat, aspect="auto", cmap="Reds", vmin=0, vmax=VMAX,
                           interpolation="nearest")
            ims.append(im)
            divider = make_axes_locatable(ax)
            # column annotation strip on top: isolation source
            cax = divider.append_axes("top", size="2.5%", pad=0.06)
            srcs = [sub.set_index("SampleID").Source[i] for i in order]
            bar = np.array([[SRC_CODE.get(s, 1) for s in srcs]], dtype=float)
            cax.imshow(bar, aspect="auto", cmap=SRC_CMAP,
                       interpolation="nearest", vmin=-0.5, vmax=2.5)
            cax.set_xticks([]); cax.set_yticks([])
            # row annotation strip on the left: direction of the association
            rax = divider.append_axes("left", size="2.5%", pad=0.06)
            dirs = np.array([[1.0 if d > 0 else 0.0 for d in
                              (t.clinical_rate - t.nonclinical_rate).values]]).T
            rax.imshow(dirs, aspect="auto", cmap=matplotlib.colors.ListedColormap([C_NON, C_CLIN]),
                       interpolation="nearest", vmin=-0.5, vmax=1.5)
            rax.set_xticks([]); rax.set_yticks([])
            # clinical / non-clinical separator
            nclin = int((sub.Source == "med").sum())
            ax.axvline(nclin - 0.5, color="black", lw=1.0)
            ax.set_yticks([])
            ax.set_xticks([])
            n = len(t)
            cax.set_title(f"{letters[ri*2+ci]}  {sp_label}: {title} (n = {n})",
                          fontsize=9, loc="left", fontweight="bold")
            ax.set_ylabel("orthogroups, ordered by q", fontsize=7.5)

    cbar_ax = fig.add_axes([0.905, 0.11, 0.015, 0.76])
    cb = fig.colorbar(ims[0], cax=cbar_ax)
    cb.set_label("gene copies\nin the orthogroup", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    fig.text(0.52, 0.015, "columns are isolates; top strip marks isolation source "
             "(red, clinical; blue, non-clinical; grey, unspecified source); left strip marks "
             "whether the orthogroup is more prevalent in clinical (red) or non-clinical (blue) "
             "isolates",
             fontsize=7.5, ha="center")
    fig.subplots_adjust(left=0.045, right=0.895, top=0.975, bottom=0.055,
                         wspace=0.18, hspace=0.42)
    pdf = os.path.join(OUTDIR, BASE + ".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, BASE + ".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(pdf, ROOT)} ({os.path.getsize(pdf)/1024:.1f} KB)")
    print(f"panels: eff solani {len(eff[eff.species=='solani'])}, res solani "
          f"{len(res[res.species=='solani'])}, eff falciformis {len(eff[eff.species=='falciformis'])}, "
          f"res falciformis {len(res[res.species=='falciformis'])}")


if __name__ == "__main__":
    main()
