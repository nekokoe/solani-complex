#!/usr/bin/env python3
"""Site-level SNP association screen — the missing generator, now with an AUC arm.

The published `snp_phenotype_assoc_summary.tsv` and `snp_vs_hog_crossvalidation.tsv` had
no writer anywhere in the repository (only the cached arrays survived in /tmp). This
script reconstructs that screen and adds a parallel run on survival-AUC phenotypes.

Site filter (reproduces the published 345,988 sites exactly):
    21 experimental strains complete (no missing genotype) and MAF >= 0.10

Phenotypes
    mortality : the seven non-degenerate traits (the two 10^5 mortality rates are excluded
                because their IQR is 0) -> 7 x 345,988 = 2,421,916 tests, as published
    auc       : nine traits from 1 - AUC (all non-degenerate) -> 9 x 345,988 tests

Outputs
    07.lethal_groups/genomewide_assoc/snp_phenotype_assoc_mortality_rebuilt.tsv
    07.lethal_groups/genomewide_assoc/snp_phenotype_assoc_auc.tsv

Cache
    The genotype matrix is read from `snp_cache/` inside the project so the run no longer
    depends on /tmp; `--cache-dir` overrides it.
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "07.lethal_groups/genomewide_assoc")
PA = os.path.join(BASE, "07.lethal_groups/phenotype_analysis")
CONDS = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
N_PERM = 300
MAF = 0.10


def load_sites(cache_dir, strains):
    geno = np.load(os.path.join(cache_dir, "snp21.npy")).astype(np.float32)
    labels = [str(x) for x in np.load(os.path.join(cache_dir, "snpdist_labels.npy"),
                                      allow_pickle=True)]
    assert len(labels) == geno.shape[1], (len(labels), geno.shape)
    idx = [labels.index(s) for s in strains]
    geno = geno[:, idx]

    complete = (geno >= 0).all(axis=1)
    p = geno.mean(axis=1) / 2.0
    maf = np.minimum(p, 1 - p)
    keep = complete & (maf >= MAF)
    print(f"sites: {geno.shape[0]} total -> {int(complete.sum())} complete -> "
          f"{int(keep.sum())} with MAF >= {MAF}")
    return geno[keep], keep


def build_traits(kind, strains):
    if kind == "mortality":
        m = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
        traits = {f"mort_{c}": m[f"mort_{c}"] for c in CONDS}
        traits.update({k: m[k] for k in ("mort_28_mean", "mort_37_mean", "temp_effect")})
        traits = {k: v for k, v in traits.items()
                  if (v.quantile(.75) - v.quantile(.25)) > 0}
        dropped = [k for k in ("mort_28Ce5", "mort_37Ce5") if k not in traits]
        print(f"mortality traits: {len(traits)} (dropped as degenerate: {dropped})")
    else:
        a = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
        traits = {f"vir_{c}": 1.0 - a[f"auc_{c}"] for c in CONDS}
        traits["vir_28_mean"] = sum(traits[f"vir_{c}"] for c in CONDS[:3]) / 3
        traits["vir_37_mean"] = sum(traits[f"vir_{c}"] for c in CONDS[3:]) / 3
        traits["vir_temp_diff"] = traits["vir_28_mean"] - traits["vir_37_mean"]
        print(f"AUC traits: {len(traits)}")
    return {k: v.reindex(strains).values for k, v in traits.items()}


def rho_p(ranks, rcen, y):
    n = len(y)
    ry = rankdata(y)
    yc = ry - ry.mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (rcen * yc).sum(axis=1) / np.sqrt((rcen ** 2).sum(axis=1) * (yc ** 2).sum())
    rho = np.clip(rho, -1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = rho * np.sqrt((n - 2) / (1 - rho ** 2))
        p = 2 * stats.t.sf(np.abs(tt), df=n - 2)
    return rho, np.nan_to_num(p, nan=1.0)


def bh(p):
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.empty(n)
    q[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", default=os.path.join(OUT, "snp_cache"))
    ap.add_argument("--phenotype", choices=["mortality", "auc", "both"], default="both")
    ap.add_argument("--permutations", type=int, default=N_PERM)
    args = ap.parse_args()

    ph = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
    strains = list(ph.index)
    geno, _ = load_sites(args.cache_dir, strains)
    ranks = rankdata(geno, axis=1).astype(np.float32)
    rcen = ranks - ranks.mean(axis=1, keepdims=True)
    print(f"rank matrix: {ranks.shape}")

    kinds = ["mortality", "auc"] if args.phenotype == "both" else [args.phenotype]
    rng = np.random.default_rng(42)
    for kind in kinds:
        traits = build_traits(kind, strains)
        rows = []
        for name, y in traits.items():
            rho, p = rho_p(ranks, rcen, y)
            q = bh(p)
            null_min = np.empty(args.permutations)
            for i in range(args.permutations):
                _, pv = rho_p(ranks, rcen, rng.permutation(y))
                null_min[i] = pv.min()
            thr = float(np.quantile(null_min, 0.05))
            rows.append(dict(trait=name, n_snp=int(ranks.shape[0]),
                             min_p=float(p.min()),
                             n_p1e4=int((p < 1e-4).sum()),
                             n_p1e5=int((p < 1e-5).sum()),
                             n_q05=int((q < 0.05).sum()),
                             perm_thr=thr,
                             n_pass_perm=int((p < thr).sum()),
                             null_min_median=float(np.median(null_min))))
            print(f"  {name:16s} min_p {p.min():.3e}  n_p<1e-4 {(p < 1e-4).sum():3d}  "
                  f"q<0.05 {(q < 0.05).sum():3d}  perm_thr {thr:.2e}  pass {rows[-1]['n_pass_perm']}")
        out = os.path.join(OUT, f"snp_phenotype_assoc_{kind}_rebuilt.tsv" if kind == "mortality"
                           else "snp_phenotype_assoc_auc.tsv")
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        print(f"wrote {out}\n")


if __name__ == "__main__":
    main()
