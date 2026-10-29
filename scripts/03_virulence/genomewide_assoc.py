#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 3: genome-wide association of HOG copy number with virulence phenotypes.

SUPERSEDED as the manuscript's universe on 2026-10-22: the phenotypes are now
1 - AUC at 10^3 and 10^4 with their two-dose means and difference, produced by
`assoc_auc_no10e5`-style runs into `hog_phenotype_assoc_auc_no10e5.tsv`. This
mortality version is kept because the archived supplementary figures S24/S36 were
drawn from it, and because its `tested` flag documents the ceiling-saturation filter.

Approach (small-n robust): per-condition Spearman of HOG copy number vs
mortality; permutation-based FDR; also association with derived traits
(temperature effect, concentration slope)."""

import os
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
OUT = "07.lethal_groups/genomewide_assoc"
os.makedirs(OUT, exist_ok=True)

pheno = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)
# read the strain list from the phenotype table instead of hard-coding it, so the
# two can never drift apart
STRAINS = pheno.index.tolist()
print(f"strains: {len(STRAINS)}")
conds = ["28Ce3","28Ce4","28Ce5","37Ce3","37Ce4","37Ce5"]
traits_all = [f"mort_{c}" for c in conds] + ["mort_28_mean","mort_37_mean",
              "temp_effect"]
# Drop ceiling-saturated traits: at 10^5 conidia nearly every strain is dead, so
# mort_28Ce5 has 2 distinct values and mort_37Ce5 4, both with IQR 0. Spearman
# against such a vector is driven by ties and produced 4,482 "hits" (= every HOG)
# in the first 21-strain run. Require a non-zero interquartile range.
def saturated(col):
    v = pheno[col].dropna()
    return (v.quantile(0.75) - v.quantile(0.25)) == 0

dropped = [t for t in traits_all if saturated(t)]
traits = [t for t in traits_all if t not in dropped]
print(f"traits tested: {len(traits)}; dropped as ceiling-saturated (IQR = 0): {dropped}")
for t in traits_all:
    print(f"   {t:16s} distinct {pheno[t].nunique():3d}  IQR {pheno[t].quantile(.75)-pheno[t].quantile(.25):.3f}")

# load HOG copy-number matrix for the 15 strains
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv",
                  sep="\t")
print("hog_count cols:", hog.columns.tolist()[:5], "...")
# long format: SampleID / HOG_ID / GeneCount
missing = [s for s in STRAINS if s not in set(hog.iloc[:, 0])]
print("missing strains in hog_count:", missing)
sub = hog[hog.iloc[:, 0].isin(STRAINS)]
mat = sub.pivot(index="HOG_ID", columns="SampleID", values="GeneCount")
mat = mat.reindex(columns=STRAINS)
print("HOG matrix:", mat.shape, "| zero fraction:", float((mat == 0).mean().mean()).__round__(3))

# filter: variable HOGs (present in >=3 strains, count varies, and absent from at
# least 2 strains). The upper bound follows n -- it used to be the literal 13 for
# 15 strains (i.e. n-2) and must move with the cohort.
n_pres = (mat > 0).sum(axis=1)
var = mat[STRAINS].nunique(axis=1) > 1
UPPER = len(STRAINS) - 2
keep = mat[(n_pres >= 3) & (n_pres <= UPPER) & var]
print(f"variable HOGs kept: {len(keep)}  (presence between 3 and {UPPER} of {len(STRAINS)})")

# association: Spearman copy vs trait
# vectorised: rank each row (ties averaged, as scipy.stats.rankdata does) and a
# fixed trait vector, then the t-approximation p -- scipy's own default for n>10.
# Needed because the relaxed presence filter makes the test set far larger than
# the old 3,630, and the permutation loop below would otherwise be intractable.
from scipy.stats import rankdata

def rank_rows(M):
    """Row-wise average ranks - the HOG rank matrix, which never changes."""
    return np.apply_along_axis(rankdata, 1, np.asarray(M, dtype=float))


def spearman_from_ranks(R, y, Rc=None):
    """Spearman rho + two-sided p for every row of precomputed ranks.

    Factoring the ranking out of the permutation loop is what makes 200
    permutations tractable: only the trait vector is re-ranked per draw.
    """
    R = np.asarray(R, dtype=float)
    n = len(y)
    if Rc is None:
        Rc = R - R.mean(axis=1, keepdims=True)
    ry = rankdata(y)
    Rc = R - R.mean(axis=1, keepdims=True)
    yc = ry - ry.mean()
    den = np.sqrt((Rc ** 2).sum(axis=1) * (yc ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (Rc * yc).sum(axis=1) / den
    rho = np.clip(rho, -1.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = rho * np.sqrt((n - 2) / (1 - rho ** 2))
        p = 2 * stats.t.sf(np.abs(tt), df=n - 2)
    return rho, np.nan_to_num(p, nan=1.0)

KEEP = keep[STRAINS].to_numpy()
RANKS = rank_rows(KEEP)                      # computed once, reused by every draw
RC = RANKS - RANKS.mean(axis=1, keepdims=True)
_rho, _p = spearman_from_ranks(rank_rows(KEEP[:300]),
                               pheno[traits[0]].reindex(STRAINS).values)

_chk_rho = [_rho[i] - stats.spearmanr(KEEP[i], pheno[traits[0]].reindex(STRAINS).values).statistic
            for i in range(50)]
_chk_p = [_p[i] - stats.spearmanr(KEEP[i], pheno[traits[0]].reindex(STRAINS).values).pvalue
          for i in range(50)]
print(f"vectorised Spearman vs scipy (50 rows): max |d rho| = {max(map(abs,_chk_rho)):.2e}, "
      f"max |d p| = {max(map(abs,_chk_p)):.2e}")
assert max(map(abs, _chk_rho)) < 1e-10 and max(map(abs, _chk_p)) < 1e-10

rows = []
y_cache = {t: pheno[t].reindex(STRAINS) for t in traits_all}
for t in traits_all:
    y = y_cache[t].values
    rho, pv = spearman_from_ranks(RANKS, y, RC)
    rows.append(pd.DataFrame({"trait": t, "rho": rho, "p": pv},
                             index=keep.index))
assoc = pd.concat(rows)
assoc.index.name = "HOG_ID"
assoc = assoc.reset_index()

# permutation FDR per trait (200 permutations of trait labels)
rng = np.random.default_rng(42)
N_PERM = 200
perm_null = {}
for t in traits:
    y = y_cache[t].values
    pv = np.empty((N_PERM, len(keep)))
    for i in range(N_PERM):
        _, r = spearman_from_ranks(RANKS, rng.permutation(y), RC)  # all HOGs, no cap
        pv[i] = r
    perm_null[t] = pv
    print(f"{t}: permutation min-p 5% = {np.quantile(pv.min(axis=1), 0.05):.4f}, "
          f"min = {pv.min():.2e}")

# Per-test FDR is plain Benjamini-Hochberg. The previous version multiplied the BH
# q by (observed min-p / permutation 5% min-p), but with thousands of tests the
# permutation minimum collapses to its 1e-8 floor, the multiplier explodes and
# every q saturates to 1 - i.e. the old "q < 0.1: 0" was an artefact of the guard,
# not evidence of no association. The permutation evidence is now reported as a
# per-trait FDR curve instead of a multiplier.
def fdr_by_trait(g):
    p = g.p.values
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(1, ranked)
    return q

# rho/p are reported for every trait (the figure scripts display all of them), but
# the FDR decision covers only the traits that passed the variance screen; the
# ceiling-saturated ones keep their rho/p and get q = NaN plus tested = False.
assoc["tested"] = assoc.trait.isin(traits)
_q = (assoc[assoc.tested].groupby("trait", group_keys=False)
      .apply(lambda g: pd.Series(fdr_by_trait(g), index=g.index)))
assoc["q"] = _q
assoc = assoc.sort_values(["tested", "q"], ascending=[False, True])
assoc.to_csv(f"{OUT}/hog_phenotype_assoc.tsv", sep="\t", index=False)
sig = assoc[assoc.q < 0.1]
print(f"\ntotal tests: {len(assoc)}; q<0.1: {len(sig)}")
print(sig.groupby("trait").size())

# permutation-based FDR of each trait's hit set (expected null hits / observed)
fdr_rows = []
for t in traits:
    p_t = assoc[assoc.trait == t].p.values
    for alpha in (0.001, 0.005, 0.01, 0.05):
        obs = int((p_t <= alpha).sum())
        null_mean = float((perm_null[t] <= alpha).sum(axis=1).mean())
        fdr_rows.append({"trait": t, "alpha": alpha, "observed_hits": obs,
                         "expected_null_hits": null_mean,
                         "perm_FDR": (null_mean / obs) if obs else np.nan})
fdf = pd.DataFrame(fdr_rows)
fdf.to_csv(f"{OUT}/permutation_fdr.tsv", sep="\t", index=False)
print("\n=== permutation FDR of each trait's hit set ===")
print(fdf.round(4).to_string(index=False))
# top hits per trait
print("\ntop hit per trait:")
print(assoc.sort_values("p").groupby("trait").head(1)[["HOG_ID","trait","rho","p","q"]]
      .round(4).to_string(index=False))

# annotate top hits
annot = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv",
                    sep="\t")
if "HOG_ID" in annot.columns:
    sig2 = assoc[assoc.q < 0.1].merge(
        annot[["HOG_ID","Gene_Name","NOG_Description","Pfam_Domain","KEGG_KO","PHI_Category"]],
        on="HOG_ID", how="left")
    sig2 = sig2.sort_values("q")
    sig2.to_csv(f"{OUT}/hog_phenotype_assoc.sig.tsv", sep="\t", index=False)
    print("\nsignificant hits annotated (top 15):")
    print(sig2.head(15)[["HOG_ID","trait","rho","q","Gene_Name","NOG_Description","Pfam_Domain"]]
          .to_string(index=False))
print("\ndone ->", OUT)
