#!/usr/bin/env python3
"""Re-run the association and the S7 modules with the 10^5 dose removed everywhere.

The 10^5 condition is on the ceiling (16 of 21 isolates reach >= 95% mortality at both
temperatures), so it was already excluded from the FDR-tested GWAS traits. Two places
still carried it indirectly: the two temperature MEANS, which average all three doses,
and the S7 module definition, whose trait list included the 10^5 rates. This script
rebuilds both without 10^5 and prints an old-vs-new comparison so the manuscript
numbers can be updated from data rather than from memory.

Reduced phenotype set (7 traits, same count as before -> the test count stays 31,388)
    mort_{28,37}Ce3        day-10 mortality at 10^3
    mort_{28,37}Ce4        day-10 mortality at 10^4
    mort_{28,37}_mean2     mean of the two doses above
    temp_effect2           mort_28_mean2 - mort_37_mean2

Writes
    07.lethal_groups/phenotype_analysis/phenotype_metrics_no10e5.tsv
    07.lethal_groups/genomewide_assoc/hog_phenotype_assoc_no10e5.tsv
    07.lethal_groups/genomewide_assoc/S7_modules_no10e5.tsv
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


def reduced_phenotypes():
    m = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
    out = pd.DataFrame(index=m.index)
    for t in ("28", "37"):
        out[f"mort_{t}Ce3"] = m[f"mort_{t}Ce3"]
        out[f"mort_{t}Ce4"] = m[f"mort_{t}Ce4"]
        out[f"mort_{t}_mean2"] = (m[f"mort_{t}Ce3"] + m[f"mort_{t}Ce4"]) / 2
    out["temp_effect2"] = out["mort_28_mean2"] - out["mort_37_mean2"]
    out.to_csv(f"{PA}/phenotype_metrics_no10e5.tsv", sep="\t")
    return out, m


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


def main():
    red, full = reduced_phenotypes()
    strains = list(red.index)
    print(f"reduced phenotype table: {red.shape[0]} strains x {red.shape[1]} traits")
    print(f"  temp_effect2 vs temp_effect Spearman = "
          f"{stats.spearmanr(red.temp_effect2, full.temp_effect).statistic:.3f}")
    print(f"  mean2 vs mean(3 doses): 28 °C "
          f"{stats.spearmanr(red.mort_28_mean2, full.mort_28_mean).statistic:.3f}, 37 °C "
          f"{stats.spearmanr(red.mort_37_mean2, full.mort_37_mean).statistic:.3f}")

    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = (hog[hog.SampleID.isin(strains)].pivot(index="HOG_ID", columns="SampleID",
                                                 values="GeneCount")
           .reindex(columns=strains).fillna(0))
    n_pres = (mat > 0).sum(axis=1)
    var = mat[strains].nunique(axis=1) > 1
    keep = mat[(n_pres >= 3) & (n_pres <= len(strains) - 2) & var]
    R = rank_rows(keep[strains].to_numpy())
    print(f"\nvariable HOGs: {len(keep)}")

    traits = list(red.columns)
    rows = []
    for t in traits:
        rho, p = rho_p(R, red[t].values)
        rows.append(pd.DataFrame({"trait": t, "rho": rho, "p": p}, index=keep.index))
    assoc = pd.concat(rows).reset_index()
    assoc["q"] = (assoc.groupby("trait", group_keys=False)
                  .apply(lambda g: pd.Series(bh(g.p.values), index=g.index))
                  .reset_index(level=0, drop=True).reindex(assoc.index))
    assoc.to_csv(f"{GW}/hog_phenotype_assoc_no10e5.tsv", sep="\t", index=False)
    print(f"tests: {len(assoc)}  min q = {assoc.q.min():.4f}  q<0.05 = "
          f"{int((assoc.q < 0.05).sum())}  p<=1e-3 = {int((assoc.p <= 1e-3).sum())}")

    # permutation null for the strongest trait
    best = assoc.groupby("trait").p.min().idxmin()
    rng = np.random.default_rng(42)
    obs = int((assoc[assoc.trait == best].p <= 1e-3).sum())
    null = np.mean([int((rho_p(R, rng.permutation(red[best].values))[1] <= 1e-3).sum())
                    for _ in range(N_PERM)])
    print(f"strongest trait {best}: {obs} hits at p<=1e-3, permutation null "
          f"{null:.1f} -> perm FDR {null / obs if obs else float('nan'):.2f}")

    print("\n" + "=" * 74)
    print("S7 MODULES without 10^5  (|rho| >= 0.50 against that temperature's traits)")
    print("=" * 74)
    annot = pd.read_csv(ANNOT, sep="\t", dtype=str).drop_duplicates("HOG_ID").set_index("HOG_ID")

    def doms(h):
        v = annot.at[h, "Pfam_Domain"] if h in annot.index else ""
        return [d for d in str(v).split(";") if d and d.lower() not in ("nan", "-")]

    universe = set(keep.index)
    bg = Counter(d for h in keep.index for d in set(doms(h)))
    out = []
    for tag, tl in (("28C", ["mort_28Ce3", "mort_28Ce4", "mort_28_mean2"]),
                    ("37C", ["mort_37Ce3", "mort_37Ce4", "mort_37_mean2"])):
        mx = (assoc[assoc.trait.isin(tl)].groupby("HOG_ID")["rho"]
              .apply(lambda x: x.abs().max()))
        mod = sorted(set(mx[mx >= RHO_MIN].index) & universe)
        sc = Counter(d for h in mod for d in set(doms(h)))
        prow = []
        for d, k in sc.items():
            m = bg.get(d, 0)
            prow.append((d, k, m, hyper_sf(k, m, len(keep), len(mod))))
        qs = bh([x[3] for x in prow])
        prow = sorted([(d, k, m, p_, q) for (d, k, m, p_), q in zip(prow, qs)],
                      key=lambda x: x[3])
        print(f"  {tag}: {len(mod)} HOGs (was {'304' if tag == '28C' else '228'}) · "
              f"{len(prow)} Pfam domains · best q = {min(x[4] for x in prow):.3f}")
        for d, k, m, p_, q in prow[:3]:
            print(f"      {d:16s} {k:3d}/{len(mod)} · bg {m:4d}/{len(keep)} · "
                  f"P = {p_:.4f} · q = {q:.3f}")
        out += [dict(module=tag, HOG_ID=h) for h in mod]
    pd.DataFrame(out).to_csv(f"{GW}/S7_modules_no10e5.tsv", sep="\t", index=False)
    print(f"\nwrote {GW}/S7_modules_no10e5.tsv")


if __name__ == "__main__":
    main()
