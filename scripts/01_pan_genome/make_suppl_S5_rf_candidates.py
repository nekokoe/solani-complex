#!/usr/bin/env python3
"""Supplementary Figure S5 - the random-forest candidate orthogroups: how many, how prevalent,
and what they encode.

Four panels, replacing a figure that had no surviving generator on disk and whose functional
panel classified orthogroups with a "first matching keyword over the union of a HOG's domains"
rule. That rule tracks orthogroup size, so it put 50-64% of the candidates in a single
"Transporter" bin; labelled by the majority class among each orthogroup's own annotated member
genes, transporters are 5.5% of the candidates.

  a      candidates as a fraction of the 29,080 orthogroups
  b, c   for each candidate, the fraction of clinical isolates carrying it against the fraction
         of non-clinical isolates carrying it, within its own species
  d      the candidates split by whether they are more prevalent in clinical or in
         non-clinical isolates, and counted in five mutually exclusive categories
         (priority: effector > kinase > transporter > regulatory > other)

Writes 00.manuscript/supplementary/Supplementary-Fig-S5-rf-clinical-HOGs.{pdf,png}
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR = os.path.join(ROOT, "00.manuscript/supplementary")
BASE = "Supplementary-Fig-S5-rf-clinical-HOGs"

PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
GROUPS = os.path.join(ROOT, "05.pan-genome/accessory_genes/bmu_ref_groups.shorter.txt")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
LABELS = os.path.join(ROOT, "02.annotation/hog_function_labels.tsv")
STATS = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun/rf_pav_hog_stats.tsv")

SPECIES = [("falciformis", "F. falciformis"), ("solani", "F. solani")]
N_HOG = 29080
C_CLIN = "#C92626"
C_NON = "#4D8AC1"
CATS = ["Effector", "Kinase", "Transporter", "Regulatory", "Other"]


def main():
    grp = pd.read_csv(GROUPS, sep="\t")
    grp = grp[grp.MedEnv.isin(["med", "env"])]
    src = grp.set_index("BMU_ID").MedEnv.map({"med": "clinical", "env": "non-clinical"})
    sp_of = grp.set_index("BMU_ID").Species

    pav = pd.read_csv(PAV, sep="\t")
    lab = pd.read_csv(LABELS, sep="\t").set_index("HOG_ID").go_mf_class
    st = pd.read_csv(STATS, sep="\t").set_index("HOG_ID")

    plt.rcParams.update({"font.family": "sans-serif", "pdf.fonttype": 42,
                         "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(9.6, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.25], hspace=0.75, wspace=0.42)

    # ---- a: fraction of the pan-genome ----
    ax = fig.add_subplot(gs[0, 0])
    cand = {sp: pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"), sep="\t")
            for sp, _ in SPECIES}
    frac = [len(cand[sp]) / N_HOG * 100 for sp, _ in SPECIES]
    ax.bar([s for _, s in SPECIES], frac, color=[C_NON, C_CLIN], width=0.55)
    for i, (sp, _) in enumerate(SPECIES):
        ax.text(i, frac[i] + 0.03, f"{len(cand[sp])} / {N_HOG:,}\n{frac[i]:.2f}%",
                ha="center", va="bottom", fontsize=7)
    ax.set_ylim(0, max(frac) * 1.45)
    ax.set_ylabel("% of pan-genome", fontsize=7.5)
    ax.set_title("a  Candidates as a fraction\nof the pan-genome", fontsize=8.5, loc="left",
                 fontweight="bold")
    ax.tick_params(labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    # ---- b, c: prevalence in clinical vs non-clinical isolates ----
    preval = {}
    for sp, sp_label in SPECIES:
        hogs = set(cand[sp].HOG_ID)
        genomes = sp_of[sp_of == sp].index
        sub = pav[pav.SampleID.isin(genomes) & pav.HOG_ID.isin(hogs)]
        sub = sub.assign(source=sub.SampleID.map(src))
        prev = (sub.groupby(["HOG_ID", "source"]).Present.mean().unstack())
        preval[sp] = prev
    for ci, (sp, sp_label) in enumerate(SPECIES):
        ax = fig.add_subplot(gs[0, ci + 1])
        prev = preval[sp]
        data = [prev["clinical"].dropna().values, prev["non-clinical"].dropna().values]
        bp = ax.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=0.9))
        for patch, col in zip(bp["boxes"], (C_CLIN, C_NON)):
            patch.set_facecolor(col)
            patch.set_alpha(0.75)
            patch.set_edgecolor("#333333")
        ax.set_xticklabels(["clinical", "non-clinical"], fontsize=7.5)
        ax.set_ylabel("Prevalence among isolates\nof the species", fontsize=7)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(f"{'bc'[ci]}  Prevalence in {sp_label}\n"
                     f"({len(prev)} candidate orthogroups)", fontsize=8.5, loc="left",
                     fontweight="bold")
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    # ---- d: composition by clinical / non-clinical dominance ----
    def category(hog, row):
        if bool(row.IsEffector):
            return "Effector"
        if bool(row.IsKinase):
            return "Kinase"
        c = lab.get(hog, "")
        if c == "Translocase / transporter activity (EC 7)":
            return "Transporter"
        if c == "Transcription regulator activity":
            return "Regulatory"
        return "Other"

    rows = []
    for sp, sp_label in SPECIES:
        for _, r in cand[sp].iterrows():
            rows.append(dict(species=sp_label,
                             group="clinical-dominant" if r.rate_delta > 0
                             else "non-clinical-dominant",
                             cat=category(r.HOG_ID, st.loc[r.HOG_ID])))
    comp = pd.DataFrame(rows)
    ax = fig.add_subplot(gs[1, :])
    width = 0.19
    groups = [f"{lbl} — {g}" for _, lbl in SPECIES
              for g in ("clinical-dominant", "non-clinical-dominant")]
    xs = np.arange(len(groups))
    for k, cat in enumerate(CATS):
        vals = [int(((comp.species == lbl) & (comp.group == g) & (comp.cat == cat)).sum())
                for _, lbl in SPECIES
                for g in ("clinical-dominant", "non-clinical-dominant")]
        ax.bar(xs + (k - 2) * width, vals, width=width, label=cat)
        for x, v in zip(xs + (k - 2) * width, vals):
            if v:
                ax.text(x, v + 1.5, str(v), ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(xs)
    ax.set_xticklabels([g.replace(" — ", "\n") for g in groups], fontsize=7)
    ax.set_ylabel("Number of candidate orthogroups", fontsize=7.5)
    ax.set_title("d  Composition of the candidates, split by whether they are more prevalent in "
                 "clinical or in non-clinical isolates", fontsize=8.5, loc="left",
                 fontweight="bold")
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(fontsize=7, frameon=False, ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.18))
    ax.spines[["top", "right"]].set_visible(False)

    pdf = os.path.join(OUTDIR, BASE + ".pdf")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(os.path.join(OUTDIR, BASE + ".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {os.path.relpath(pdf, ROOT)} ({os.path.getsize(pdf)/1024:.1f} KB)")
    print("\npanel d counts:")
    print(pd.crosstab([comp.species, comp.group], comp.cat).reindex(columns=CATS).to_string())
    print("\ncandidate category shares overall:",
          (comp.cat.value_counts(normalize=True) * 100).round(1).to_dict())


if __name__ == "__main__":
    main()
