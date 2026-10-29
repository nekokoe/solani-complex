#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 plink .raw 基因型构建简单匹配距离矩阵 (1 - matching proportion)，供 pegas amova。
输入: tmp/amova/<sp>_geno.raw (plink --recode A)
输出: tmp/amova/<sp>_dist.npy + <sp>_groups.tsv
"""
import sys
import numpy as np
import pandas as pd

sp = sys.argv[1]
raw_path = f"tmp/amova/{sp}_geno.raw"

# read header only for IID
with open(raw_path) as f:
    header = f.readline().split()
# sample IDs
iids = []
for line in open(raw_path):
    parts = line.split()
    iids.append(parts[1])
# read genotype block as numpy (columns from index 6 onward)
geno = np.loadtxt(raw_path, dtype=np.int8, skiprows=1, usecols=range(6, len(header)))
# missing (-9) -> NaN
geno[geno < 0] = -1
n, m = geno.shape
print(f"{sp}: {n} samples x {m} SNPs")

# pairwise simple-matching distance over non-missing loci
D = np.zeros((n, n))
for i in range(n):
    gi = geno[i]
    for j in range(i + 1, n):
        gj = geno[j]
        both = (gi >= 0) & (gj >= 0)
        nb = both.sum()
        if nb == 0:
            D[i, j] = D[j, i] = 0.5
        else:
            D[i, j] = D[j, i] = 1.0 - (gi[both] == gj[both]).sum() / nb
np.save(f"tmp/amova/{sp}_dist.npy", D)
pd.DataFrame({"IID": iids}).to_csv(f"tmp/amova/{sp}_samples.tsv", sep="\t", index=False)
print(f"saved {sp} distance matrix ({n}x{n})")
