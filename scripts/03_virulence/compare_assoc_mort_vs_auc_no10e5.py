#!/usr/bin/env python3
"""Same 7-trait, 10^5-free association universe, but built on survival AUC.

The mortality version uses day-10 rates at 10^3 and 10^4, their two-dose means and the
difference of those means. This builds the exact AUC analogue — 1 - AUC at the same two
doses, the same two means and the same difference — and runs the identical screen, so the
two can be compared on equal footing.

Rationale for the AUC arm: at 10^3 and 10^4 the endpoint is not on the ceiling (10-12
distinct values, IQR 0.23-0.39), so the trajectory adds kinetic information — how fast a
strain kills — that the day-10 rate throws away. At 10^5 that kinetic information is an
artefact of an overwhelming dose and is excluded from both arms.

Writes
    07.lethal_groups/phenotype_analysis/phenotype_metrics_auc_no10e5.tsv
    07.lethal_groups/genomewide_assoc/hog_phenotype_assoc_auc_no10e5.tsv
"""

import os
from collections import Counter
from math import comb

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PA = os.path.join(BASE, "07.lethal_groups/phenotype_analysis")
GW = os.path.join(BASE, "07.lethal_groups/genomewide_assoc")
HOGCOUNT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
ANNOT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
RHO_MIN = 0.50
N_PERM = 200


def auc_traits():
    a = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    out = pd.DataFrame(index=a.index)
    for t in ("28", "37"):
        out[f"vir_{t}Ce3"] = 1.0 - a[f"auc_{t}Ce3"]
        out[f"vir_{t}Ce4"] = 1.0 - a[f"auc_{t}Ce4"]
        out[f"vir_{t}_mean2"] = (out[f"vir_{t}Ce3"] + out[f"vir_{t}Ce4"]) / 2
    out["vir_temp_effect2"] = out["vir_28_mean2"] - out["vir_37_mean2"]
    out.to_csv(f"{PA}/phenotype_metrics_auc_no10e5.tsv", sep="\t")
    return out


def rank_rows(M):
    return np.apply_along_axis(rankdata, 1, np.asarray(M, dtype=float))


def rho_p(R, y):
    n = len(y)
    Rc = R - R.mean(axis=1, keepdims=True)
    ry = rankdata(y)
    yc = ry - ry.mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (Rc * yc).sum(axis=1) / np.sqrt((Rc ** 2).sum(axis=1) * (yc ** 2).sum())
    rho = np.clip(rho, -1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = rho * np.sqrt((n - 2) / (1 - rho ** 2))
        return rho, np.nan_to_num(2 * stats.t.sf(np.abs(tt), df=n - 2), nan=1.0)


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.empty(n)
    q[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def hyper_sf(k, m, N, K):
    if k <= 0:
        return 1.0
    return sum(comb(m, x) * comb(N - m, K - x) for x in range(k, min(m, K) + 1)) / comb(N, K)


def screen(traits_df, keep_index, R, path, strains):
    rows = []
    for t in traits_df.columns:
        rho, p = rho_p(R, traits_df[t].reindex(strains).values)
        rows.append(pd.DataFrame({"trait": t, "rho": rho, "p": p}, index=keep_index))
    assoc = pd.concat(rows).reset_index().rename(columns={"index": "HOG_ID"})
    assoc["q"] = (assoc.groupby("trait", group_keys=False)
                  .apply(lambda g: pd.Series(bh(g.p.values), index=g.index))
                  .reset_index(level=0, drop=True).reindex(assoc.index))
    assoc.to_csv(path, sep="\t", index=False)
    return assoc


def modules(assoc, universe, tag_traits):
    mx = (assoc[assoc.trait.isin(tag_traits)].groupby("HOG_ID")["rho"]
          .apply(lambda x: x.abs().max()))
    return sorted(set(mx[mx >= RHO_MIN].index) & universe)


def main():
    auc = auc_traits()
    mort = pd.read_csv(f"{PA}/phenotype_metrics_no10e5.tsv", sep="\t", index_col=0)
    print(f"AUC universe: {auc.shape[1]} traits · mortality universe: {mort.shape[1]} traits")
    print(f"temperature-preference trait: Spearman = "
          f"{stats.spearmanr(auc.vir_temp_effect2, mort.temp_effect2).statistic:.3f}")

    strains = list(auc.index)
    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = (hog[hog.SampleID.isin(strains)].pivot(index="HOG_ID", columns="SampleID",
                                                 values="GeneCount")
           .reindex(columns=strains).fillna(0))
    n_pres = (mat > 0).sum(axis=1)
    var = mat[strains].nunique(axis=1) > 1
    keep = mat[(n_pres >= 3) & (n_pres <= len(strains) - 2) & var]
    R = rank_rows(keep[strains].to_numpy())
    universe = set(keep.index)
    print(f"variable HOGs: {len(keep)}")

    a_assoc = screen(auc, keep.index, R, f"{GW}/hog_phenotype_assoc_auc_no10e5.tsv", strains)
    m_assoc = pd.read_csv(f"{GW}/hog_phenotype_assoc_no10e5.tsv", sep="\t")

    print("\n" + "=" * 74)
    print("SCREEN COMPARISON (same 7 traits, same 4,484 HOGs, 10^5 excluded in both)")
    print("=" * 74)
    for name, d in (("mortality", m_assoc), ("AUC", a_assoc)):
        print(f"  {name:10s} tests {len(d)}  min q {d.q.min():.4f}  q<0.05 "
              f"{int((d.q < 0.05).sum()):2d}  p<=1e-3 {int((d.p <= 1e-3).sum()):2d}")

    print("\nper matched trait: how much do the two phenotypes agree about the HOGs?")
    pairs = [(f"mort_{t}Ce{d}", f"vir_{t}Ce{d}") for t in ("28", "37") for d in (3, 4)] + \
            [(f"mort_{t}_mean2", f"vir_{t}_mean2") for t in ("28", "37")] + \
            [("temp_effect2", "vir_temp_effect2")]
    for mt, at in pairs:
        x = m_assoc[m_assoc.trait == mt].set_index("HOG_ID").rho
        y = a_assoc[a_assoc.trait == at].set_index("HOG_ID").rho
        j = pd.concat([x, y], axis=1, keys=["m", "a"]).dropna()
        r = stats.spearmanr(j.m, j.a).statistic
        print(f"  {mt:16s} vs {at:18s} rho-of-rhos = {r:+.3f}  top20 overlap "
              f"{len(set(j.nsmallest(20,'m').index) & set(j.nsmallest(20,'a').index)):2d}/20")

    print("\npermutation null for the strongest AUC trait")
    best = a_assoc.groupby("trait").p.min().idxmin()
    rng = np.random.default_rng(42)
    obs = int((a_assoc[a_assoc.trait == best].p <= 1e-3).sum())
    null = np.mean([int((rho_p(R, rng.permutation(auc[best].reindex(strains).values))[1]
                         <= 1e-3).sum()) for _ in range(N_PERM)])
    print(f"  {best}: {obs} hits at p<=1e-3 vs null {null:.1f} -> perm FDR "
          f"{null / obs if obs else float('nan'):.2f}")

    print("\n" + "=" * 74)
    print("S7 MODULES under the AUC traits")
    print("=" * 74)
    annot = pd.read_csv(ANNOT, sep="\t", dtype=str).drop_duplicates("HOG_ID").set_index("HOG_ID")
    def doms(h):
        v = annot.at[h, "Pfam_Domain"] if h in annot.index else ""
        return [d for d in str(v).split(";") if d and d.lower() not in ("nan", "-")]
    bg = Counter(d for h in keep.index for d in set(doms(h)))
    for tag, tl, ref in (("28C", ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"], 176),
                         ("37C", ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"], 155)):
        mod = modules(a_assoc, universe, tl)
        sc = Counter(d for h in mod for d in set(doms(h)))
        prow = [(d, k, bg.get(d, 0), hyper_sf(k, bg.get(d, 0), len(keep), len(mod)))
                for d, k in sc.items()]
        qs = bh([x[3] for x in prow])
        best = min(zip(prow, qs), key=lambda z: z[1])
        print(f"  {tag}: {len(mod)} HOGs (mortality arm: {ref}) · best domain {best[0][0]} "
              f"{best[0][1]}/{len(mod)} · q = {best[1]:.3f}")
    print(f"\nwrote {GW}/hog_phenotype_assoc_auc_no10e5.tsv")


if __name__ == "__main__":
    main()
