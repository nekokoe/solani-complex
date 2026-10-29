#!/usr/bin/env python3
"""Effector genes in the FSSC pan-genome and in the clinical-versus-environmental comparison.

The effector axis is reported GENE-LEVEL throughout, in the same three questions as the
antifungal-resistance axis and with the same cohort:

  1. effector genes as a share of ALL genes in each presence interval (gene level; the HOG-level
     "contains at least one effector" share inflates with orthogroup size and is not used);
  2. the number of effector genes per genome in clinical versus non-clinical isolates;
  3. the effector-gene share inside the random-forest candidate orthogroups versus elsewhere.

Effector calls come from the EffectorP-3 consensus (08.pan_020_genes/effector_master.gene_level.tsv;
note that EffectorP ran on the SignalP-positive secretome, so "effector" and "non-effector" are
explicit calls only for those proteins). Four isolates have no EffectorP output at all
(BMU12683-BMU12686) and are excluded from the per-genome comparison, because their zero is a
missing annotation rather than an absence of effectors. Source and species labels are those of the
random-forest cohort (09.rf_clinical_genes/rf_genes.per_genome.tsv): 60 clinical / 26 non-clinical
for F. solani and 48 / 18 for F. falciformis.

Writes 08.pan_020_genes/effector_genes_by_interval.tsv
       08.pan_020_genes/effector_genes_per_genome.tsv
       08.pan_020_genes/effector_source_assoc.tsv  (per-orthogroup presence test, used only for
        the Supplementary Fig. S41 heatmap; the manuscript's effector statements stay gene-level)
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENE_EFF = os.path.join(ROOT, "08.pan_020_genes/effector_master.gene_level.tsv")
IV = os.path.join(ROOT, "05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv")
GENE_HOG = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv")
PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
HOG_EFF = os.path.join(ROOT, "08.pan_020_genes/effector_master.hog_level.tsv")
PER_GENOME = os.path.join(ROOT, "09.rf_clinical_genes/rf_genes.per_genome.tsv")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
OUT_DIR = os.path.join(ROOT, "08.pan_020_genes")
BINS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]
SPECIES = ["solani", "falciformis"]

def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, p[order[i]] * n / (i + 1))
        q[order[i]] = prev
    return np.clip(q, 0, 1)



def main():
    eff = pd.read_csv(GENE_EFF, sep="\t", usecols=["GeneID", "SampleID", "IsEffector"])
    eff["IsEffector"] = eff.IsEffector.astype(str).str.lower().eq("true")
    eff = eff.drop_duplicates("GeneID")

    gvh = pd.read_csv(GENE_HOG, sep="\t", usecols=["GeneID", "HOG_ID"])
    iv = pd.read_csv(IV, sep="\t")
    iv.columns = ["HOG_ID", "rate", "bin"]
    m = gvh.merge(iv, on="HOG_ID").merge(eff, on="GeneID", how="left")
    m["IsEffector"] = m.IsEffector.fillna(False)

    # ---- 1. effector genes / ALL genes per presence interval ----
    rows = []
    for b in BINS:
        s = m[m.bin == b]
        rows.append(dict(interval=b, n_genes=len(s), n_effector_genes=int(s.IsEffector.sum()),
                         pct=round(s.IsEffector.mean() * 100, 2)))
    by_iv = pd.DataFrame(rows)
    by_iv["fold_vs_0-20"] = (by_iv.pct / by_iv.pct.iloc[0]).round(2)
    by_iv.to_csv(os.path.join(OUT_DIR, "effector_genes_by_interval.tsv"), sep="\t", index=False)
    print("effector genes per presence interval (gene level):")
    print(by_iv.to_string(index=False))

    # ---- 2. effector genes per genome, clinical versus non-clinical ----
    per = eff.groupby("SampleID").IsEffector.sum().rename("n_effector_genes")
    pg = pd.read_csv(PER_GENOME, sep="\t", usecols=["SampleID", "Source", "Species"])
    pg = pg[pg.Species.isin(SPECIES)]
    gm = pg.merge(per, on="SampleID", how="left")
    no_ann = set(gm[gm.n_effector_genes.isna()].SampleID)
    print(f"\nno EffectorP output (excluded): {sorted(no_ann)}")
    gm = gm[~gm.SampleID.isin(no_ann)]
    rows = []
    for sp in SPECIES:
        sub = gm[gm.Species == sp]
        a = sub[sub.Source == "med"].n_effector_genes
        b = sub[sub.Source != "med"].n_effector_genes
        p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
        rows.append(dict(species=sp, n_clinical=len(a), n_nonclinical=len(b),
                         median_clinical=float(a.median()), median_nonclinical=float(b.median()),
                         median_delta=float(a.median() - b.median()), p_mannwhitney=p))
        print(f"  {sp}: clinical median {a.median():.0f} vs non-clinical {b.median():.0f} "
              f"(delta {a.median()-b.median():+.0f})  P = {p:.3g}")
    per_sp = pd.DataFrame(rows)
    per_sp.to_csv(os.path.join(OUT_DIR, "effector_genes_per_genome.tsv"), sep="\t", index=False)

    # ---- 3. per-orthogroup presence test, clinical versus non-clinical (heatmap display only) ----
    hog_eff = pd.read_csv(HOG_EFF, sep="\t", usecols=["HOG_ID", "EffectorHOG"])
    eff_hogs = set(hog_eff[hog_eff.EffectorHOG == True].HOG_ID)
    pav = pd.read_csv(PAV, sep="\t")
    rows_h = []
    for sp in SPECIES:
        sub = pg[pg.Species == sp]
        clinical = sub[sub.Source == "med"].SampleID
        nonclin = sub[sub.Source != "med"].SampleID
        wide = (pav[pav.SampleID.isin(sub.SampleID) & pav.HOG_ID.isin(eff_hogs)]
                .groupby(["HOG_ID", "SampleID"]).Present.mean().unstack())
        for hog, r in wide.iterrows():
            a, b = r[clinical].dropna(), r[nonclin].dropna()
            if len(a) < 2 or len(b) < 2:
                continue
            try:
                p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
            rows_h.append(dict(species=sp, HOG_ID=hog,
                               clinical_rate=round(float(a.mean()), 4),
                               nonclinical_rate=round(float(b.mean()), 4),
                               p_mannwhitney=p))
    res_h = pd.DataFrame(rows_h)
    res_h["q_bh"] = np.nan
    for sp, grp in res_h.groupby("species"):
        ok = grp.p_mannwhitney.notna()
        res_h.loc[grp.index[ok], "q_bh"] = bh(grp.p_mannwhitney[ok])
    res_h.sort_values(["species", "q_bh"]).to_csv(
        os.path.join(OUT_DIR, "effector_source_assoc.tsv"), sep="\t", index=False)
    print(f"\neffector per-orthogroup presence test (heatmap only): "
          f"{res_h[res_h.q_bh < 0.05].groupby('species').size().to_dict()}")

    # ---- 4. effector-gene share inside the candidates versus elsewhere ----
    cand = set()
    for sp in SPECIES:
        cand |= set(pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"),
                                sep="\t").HOG_ID)
    m["inCand"] = m.HOG_ID.isin(cand)
    for k, lab in ((True, "candidate orthogroups"), (False, "all other orthogroups")):
        s = m[m.inCand == k]
        print(f"  {lab}: {int(s.IsEffector.sum()):,}/{len(s):,} = {s.IsEffector.mean()*100:.2f}%")


if __name__ == "__main__":
    main()
