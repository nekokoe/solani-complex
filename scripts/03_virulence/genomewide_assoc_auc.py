#!/usr/bin/env python3
"""Genome-wide HOG association screen with survival-AUC phenotypes (parallel run).

Mirror of `genomewide_assoc.py` with one difference: the phenotype is the area under the
daily survival curve rather than the day-10 mortality. Everything else — HOG matrix,
variable-HOG filter, vectorised Spearman, 200 permutations, per-trait Benjamini-Hochberg —
is identical, so the two screens are directly comparable.

    AUC   = trapz(survival, days 0..10) / 10      (0 = all dead, 1 = maximal survival)
    vir   = 1 - AUC                               (1 = maximal lethality, sign matches mortality)

Traits (9): the six conditions, the two temperature means, and the temperature difference.
None of them is degenerate: each has 18-21 distinct values, so all nine enter the FDR
instead of seven.

Outputs (all suffixed _auc so nothing overwrites the mortality screen)
    07.lethal_groups/genomewide_assoc/hog_phenotype_assoc_auc.tsv
    07.lethal_groups/genomewide_assoc/permutation_fdr_auc.tsv
"""

import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "07.lethal_groups/genomewide_assoc")
PA = os.path.join(BASE, "07.lethal_groups/phenotype_analysis")
HOGCOUNT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")

CONDS = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
N_PERM = 200


def build_traits():
    auc = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    t = pd.DataFrame(index=auc.index)
    for c in CONDS:
        t[f"vir_{c}"] = 1.0 - auc[f"auc_{c}"]
    t["vir_28_mean"] = t[[f"vir_{c}" for c in CONDS[:3]]].mean(axis=1)
    t["vir_37_mean"] = t[[f"vir_{c}" for c in CONDS[3:]]].mean(axis=1)
    t["vir_temp_diff"] = t["vir_28_mean"] - t["vir_37_mean"]   # positive = 28 °C more lethal
    return t


def rank_rows(M):
    return np.apply_along_axis(rankdata, 1, np.asarray(M, dtype=float))


def spearman_from_ranks(R, y):
    """Spearman rho + two-sided p for every row of precomputed ranks."""
    R = np.asarray(R, dtype=float)
    n = len(y)
    Rc = R - R.mean(axis=1, keepdims=True)
    ry = rankdata(y)
    yc = ry - ry.mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (Rc * yc).sum(axis=1) / np.sqrt((Rc ** 2).sum(axis=1) * (yc ** 2).sum())
    rho = np.clip(rho, -1.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = rho * np.sqrt((n - 2) / (1 - rho ** 2))
        p = 2 * stats.t.sf(np.abs(tt), df=n - 2)
    return rho, np.nan_to_num(p, nan=1.0)


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(1, ranked)
    return q


def main():
    traits_df = build_traits()
    STRAINS = traits_df.index.tolist()
    print(f"strains: {len(STRAINS)}")
    traits = list(traits_df.columns)
    for t in traits:
        v = traits_df[t]
        print(f"   {t:16s} distinct {v.nunique():3d}  IQR {v.quantile(.75) - v.quantile(.25):.3f}"
              f"  range {v.min():.2f}-{v.max():.2f}")

    hog = pd.read_csv(HOGCOUNT, sep="\t")
    sub = hog[hog.iloc[:, 0].isin(STRAINS)]
    mat = sub.pivot(index="HOG_ID", columns="SampleID", values="GeneCount").reindex(columns=STRAINS)
    n_pres = (mat > 0).sum(axis=1)
    var = mat[STRAINS].nunique(axis=1) > 1
    UPPER = len(STRAINS) - 2
    keep = mat[(n_pres >= 3) & (n_pres <= UPPER) & var]
    print(f"variable HOGs kept: {len(keep)} (presence {3}..{UPPER} of {len(STRAINS)})")

    KEEP = keep[STRAINS].to_numpy()
    RANKS = rank_rows(KEEP)

    # reproduce scipy exactly on a sample before trusting the vectorised path
    t0 = traits[0]
    y0 = traits_df[t0].reindex(STRAINS).values
    r_v, p_v = spearman_from_ranks(RANKS, y0)
    dr = [r_v[i] - stats.spearmanr(KEEP[i], y0).statistic for i in range(50)]
    dp = [p_v[i] - stats.spearmanr(KEEP[i], y0).pvalue for i in range(50)]
    print(f"vectorised vs scipy (50 rows): max |d rho| = {max(map(abs, dr)):.1e}, "
          f"max |d p| = {max(map(abs, dp)):.1e}")
    assert max(map(abs, dr)) < 1e-10 and max(map(abs, dp)) < 1e-10

    rows, perm_null = [], {}
    rng = np.random.default_rng(42)
    for t in traits:
        y = traits_df[t].reindex(STRAINS).values
        rho, pv = spearman_from_ranks(RANKS, y)
        rows.append(pd.DataFrame({"trait": t, "rho": rho, "p": pv}, index=keep.index))
        pvs = np.empty((N_PERM, len(keep)))
        for i in range(N_PERM):
            _, r = spearman_from_ranks(RANKS, rng.permutation(y))
            pvs[i] = r
        perm_null[t] = pvs
        print(f"{t}: permutation min-p 5% = {np.quantile(pvs.min(axis=1), 0.05):.4f}, "
              f"min = {pvs.min():.2e}")

    assoc = pd.concat(rows)
    assoc.index.name = "HOG_ID"
    assoc = assoc.reset_index()
    assoc["tested"] = True
    assoc["q"] = (assoc.groupby("trait", group_keys=False)
                  .apply(lambda g: pd.Series(bh(g.p.values), index=g.index))
                  .reset_index(level=0, drop=True).reindex(assoc.index))
    assoc = assoc.sort_values(["q", "p"])
    path = f"{OUT}/hog_phenotype_assoc_auc.tsv"
    assoc.to_csv(path, sep="\t", index=False)
    print(f"\ntotal rows: {len(assoc)}; tests: {len(assoc)}; q<0.1: {int((assoc.q < 0.1).sum())}")
    print(f"min q = {assoc.q.min():.4f}")

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
    fdf.to_csv(f"{OUT}/permutation_fdr_auc.tsv", sep="\t", index=False)
    print("\n=== permutation FDR of each AUC trait's hit set ===")
    print(fdf.round(4).to_string(index=False))

    print("\ntop hit per trait:")
    print(assoc.sort_values("p").groupby("trait").head(1)[["HOG_ID", "trait", "rho", "p", "q"]]
          .round(5).to_string(index=False))
    print("\nwrote", path)


if __name__ == "__main__":
    main()
