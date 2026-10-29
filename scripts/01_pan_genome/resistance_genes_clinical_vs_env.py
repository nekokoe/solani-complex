#!/usr/bin/env python3
"""Antifungal-resistance genes in the FSSC pan-genome and in the clinical-versus-environmental
comparison.

Resistance is one of the expected features of clinical isolates, so the pan-genome section has to
report it alongside the virulence (effector) axis. This project has no CARD/RGI run — and CARD is
bacterial-AMR-centric, so a negative RGI result would say nothing about fungi — therefore the
catalogue below is curated from the eggNOG preferred names, grouped by the drug class the gene
bears on:

  azole core      ERG11 (lanosterol 14-alpha demethylase, the azole target) and UPC2 (the sterol
                  regulator that controls ERG expression)
  ergosterol      the biosynthetic pathway around that target: ERG1-ERG13, ERG20, ERG24-ERG27,
                  HMG1, IDI1, MVD1
  efflux          CDR1, SNQ2 (ABC), MDR1, ATR1 (MFS), PDR6, PDR16
  other targets   FKS1 and RHO1 (echinocandin target and its regulator), FUR1 (5-fluorocytosine)

Three questions are answered:
  1. where these genes sit on the presence gradient;
  2. whether they are over-represented among the random-forest clinical candidates;
  3. whether individual resistance orthogroups differ in prevalence or copy number between
     clinical and non-clinical isolates of the same species.

Writes 08.pan_020_genes/resistance_genes_clinical_vs_env.tsv
       08.pan_020_genes/resistance_genes_by_interval.tsv
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EGGNOG = os.path.join(ROOT, "02.annotation/all_isolates/all.fusarium.gene_eggnog.tsv")
GENE_HOG = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv")
PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
IV = os.path.join(ROOT, "05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv")
# the source/species table the random-forest screen itself uses, so that the resistance
# comparison is run on exactly the same cohort (solani 60/26, falciformis 48/18)
PER_GENOME = os.path.join(ROOT, "09.rf_clinical_genes/rf_genes.per_genome.tsv")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
OUT_DIR = os.path.join(ROOT, "08.pan_020_genes")

AZOLE = {"ERG11", "UPC2"}
ERGOSTEROL = {f"ERG{i}" for i in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 20,
                                  24, 25, 26, 27)} | {"HMG1", "IDI1", "MVD1"}
EFFLUX = {"CDR1", "SNQ2", "MDR1", "ATR1", "PDR6", "PDR16"}
OTHER = {"FKS1", "RHO1", "FUR1"}
CATALOGUE = {"azole core": AZOLE, "ergosterol pathway": ERGOSTEROL,
             "efflux": EFFLUX, "other drug targets": OTHER}
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
    eg = pd.read_csv(EGGNOG, sep="\t", header=None, usecols=[0, 1],
                     names=["GeneID", "name"], dtype=str, low_memory=False)
    eg = eg.drop_duplicates("GeneID")
    eg["category"] = np.nan
    for cat, names in CATALOGUE.items():
        eg.loc[eg.name.isin(names) & eg.category.isna(), "category"] = cat
    res = eg[eg.category.notna()].copy()

    gvh = pd.read_csv(GENE_HOG, sep="\t", usecols=["GeneID", "HOG_ID"])
    res_m = res.merge(gvh, on="GeneID", how="inner")
    res_hogs = res_m.groupby("HOG_ID").category.apply(set).to_dict()
    print(f"catalogue: {len(res):,} annotated genes in {len(res_hogs):,} orthogroups "
          f"({len(res_hogs) / 29080 * 100:.1f}% of 29,080)")
    for cat in CATALOGUE:
        n = sum(1 for v in res_hogs.values() if cat in v)
        print(f"   {cat:22s} {n:5d} orthogroups")

    # ---- 1. resistance genes / ALL genes per presence interval (gene level) ----
    # The HOG-level share ("orthogroups carrying a resistance gene") inflates with orthogroup
    # size in exactly the same way the old transporter figure did, so the interval sheet counts
    # GENE COPIES: a resistance gene is one eggNOG-named gene, and the denominator is every gene
    # of the interval.
    iv = pd.read_csv(IV, sep="\t")
    iv.columns = ["HOG_ID", "rate", "bin"]
    res_genes = set(res.GeneID)
    gm = gvh.merge(iv, on="HOG_ID")
    gm["is_res"] = gm.GeneID.isin(res_genes)
    rows = []
    for b in BINS:
        sub = gm[gm.bin == b]
        rows.append(dict(interval=b, n_genes=len(sub), n_resistance_genes=int(sub.is_res.sum()),
                         pct=round(sub.is_res.mean() * 100, 2)))
    by_iv = pd.DataFrame(rows)
    by_iv["fold_vs_0-20"] = (by_iv.pct / by_iv.pct.iloc[0]).round(2)
    by_iv.to_csv(os.path.join(OUT_DIR, "resistance_genes_by_interval.tsv"), sep="\t", index=False)
    print("\nresistance genes on the presence gradient (gene level):")
    print(by_iv.to_string(index=False))

    # ---- 2. resistance genes among the candidate orthogroups, at gene level ----
    cand = set()
    for sp in SPECIES:
        cand |= set(pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"),
                                sep="\t").HOG_ID)
    gm["inCand"] = gm.HOG_ID.isin(cand)
    for k, lab in ((True, "candidate orthogroups"), (False, "all other orthogroups")):
        sub = gm[gm.inCand == k]
        print(f"  resistance genes in {lab}: {int(sub.is_res.sum()):,}/{len(sub):,} = "
              f"{sub.is_res.mean()*100:.2f}%")

    # ---- 3. clinical versus non-clinical, per species ----
    pav = pd.read_csv(PAV, sep="\t")
    pg = pd.read_csv(PER_GENOME, sep="\t", usecols=["SampleID", "Source", "Species"])
    pg["Species"] = pg.Species.fillna("")
    pg = pg[pg.Species.isin(SPECIES)]
    rows = []
    for sp in SPECIES:
        sub_pg = pg[pg.Species == sp]
        genomes = list(sub_pg.SampleID)
        clinical = list(sub_pg[sub_pg.Source == "med"].SampleID)
        nonclin = list(sub_pg[sub_pg.Source != "med"].SampleID)
        sub = pav[pav.SampleID.isin(genomes) & pav.HOG_ID.isin(res_hogs)]
        wide = sub.groupby(["HOG_ID", "SampleID"]).Present.mean().unstack()
        print(f"   {sp}: {len(genomes)} genomes "
              f"({len(clinical)} clinical / {len(nonclin)} non-clinical)")
        for hog, r in wide.iterrows():
            a, b = r[clinical].dropna(), r[nonclin].dropna()
            if not len(a) or not len(b):
                continue
            try:
                p = stats.mannwhitneyu(a, b, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
            rows.append(dict(species=sp, HOG_ID=hog, categories=";".join(sorted(res_hogs[hog])),
                             clinical_rate=a.mean(), nonclinical_rate=b.mean(),
                             rate_delta=a.mean() - b.mean(), p_mannwhitney=p))
    res_tab = pd.DataFrame(rows)
    if len(res_tab):
        res_tab["q_bh"] = np.nan
        for sp, grp in res_tab.groupby("species"):
            ok = grp.p_mannwhitney.notna()
            res_tab.loc[grp.index[ok], "q_bh"] = bh(grp.p_mannwhitney[ok])
        res_tab.sort_values(["species", "p_mannwhitney"]).to_csv(
            os.path.join(OUT_DIR, "resistance_genes_clinical_vs_env.tsv"), sep="\t", index=False)
        print("\nclinical vs non-clinical prevalence of resistance orthogroups:")
        for sp, grp in res_tab.groupby("species"):
            print(f"   {sp}: tested {len(grp)}, nominal p<0.05 "
                  f"{int((grp.p_mannwhitney < 0.05).sum())}, BH q<0.05 "
                  f"{int((grp.q_bh < 0.05).sum())}")
        top = res_tab.nsmallest(6, "p_mannwhitney")
        print(top[["species", "HOG_ID", "categories", "clinical_rate", "nonclinical_rate",
                   "p_mannwhitney", "q_bh"]].to_string(index=False))


if __name__ == "__main__":
    main()
