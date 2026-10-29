#!/usr/bin/env python3
# DEPRECATED PHENOTYPE (2026-10-22): reads the mortality-era association table
# (nine phenotypes, 10^5 included). The live universe is 1 - AUC at 10^3 and 10^4
# with the two-dose means and their difference -> hog_phenotype_assoc_auc_no10e5.tsv.
# Kept only to regenerate the archived artefacts; do not use for new figures.
# -*- coding: utf-8 -*-
"""gwa_guard_diagnostic.py — does the empirical min-p guard in
scripts/genomewide_assoc.py suppress every hit?

That script applies, per trait,

    thr = max(perm_5pct_min_p, 1e-8)
    q   = min(1, bh_q * (observed_min_p / thr))

With thousands of HOGs the permutation null's minimum p collapses towards 0, so
`thr` hits its 1e-8 floor and the multiplier becomes huge for any observed
minimum above it -> every q saturates to 1 and the reported result is "0 hits"
regardless of the data. This script recomputes BH per trait WITHOUT the guard
from the produced p-values and reports what the data actually supports, plus a
permutation-free sanity view.

Usage: python3 scripts/gwa_guard_diagnostic.py
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "07.lethal_groups", "genomewide_assoc")
assoc = pd.read_csv(f"{OUT}/hog_phenotype_assoc.tsv", sep="\t")
print(f"tests: {len(assoc)}  traits: {assoc.trait.nunique()}")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(1, ranked)
    return q


print("\n=== per trait: observed vs guarded ===")
rows = []
for t, g in assoc.groupby("trait"):
    q_bh = bh(g.p.values)
    n_guarded = int((g.q < 0.1).sum())
    n_bh10 = int((q_bh < 0.1).sum())
    n_bh25 = int((q_bh < 0.25).sum())
    # Benjamini-Yekutieli style expectation from the p-values themselves:
    # how many would we expect by chance at the same alpha?
    exp10 = 0.1 * len(g)
    rows.append({"trait": t, "n_HOGs": len(g), "min_p": g.p.min(),
                 "n_q_lt_0.1_guarded": n_guarded, "n_q_lt_0.1_BH": n_bh10,
                 "n_q_lt_0.25_BH": n_bh25, "expected_by_chance_at_0.1": exp10,
                 "max_abs_rho": g.rho.abs().max()})
r = pd.DataFrame(rows)
print(r.round(5).to_string(index=False))
r.to_csv(f"{OUT}/guard_diagnostic.tsv", sep="\t", index=False)

tot_bh1 = int((r["n_q_lt_0.1_BH"] > 0).sum())
print(f"\ntraits with any BH q<0.1 once the guard is removed: {tot_bh1}/{len(r)}")
print(f"total HOG-trait pairs at BH q<0.1: {int(r["n_q_lt_0.1_BH"].sum())} "
      f"(guarded: {int(r["n_q_lt_0.1_guarded"].sum())})")

top = assoc.reindex(assoc.p.sort_values().index).head(15)
print("\n=== top 15 by raw p (unguarded view) ===")
print(top[["HOG_ID", "trait", "rho", "p", "q"]].round(6).to_string(index=False))
