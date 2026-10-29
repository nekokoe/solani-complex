#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 2: literature pathogenicity genes vs phenotype.

Align the literature pathogenicity-gene proteins to the 21 experimental
strains' proteomes with mmseqs easy-search (or diamond fallback), then
correlate per-gene copy number / identity with mortality across conditions.

Outputs:
  - gene_strain_copy.tsv       21 strains x N genes (best-hit copies)
  - gene_strain_identity.tsv   21 strains x N genes (best identity)
  - gene_phenotype_assoc.tsv   per-gene x condition Spearman assoc + FDR
  - figure_literature_genes.pdf
"""

import os, subprocess, glob, sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
# read the strain list from the phenotype table (21 strains since 2026-10-12) so
# the two cannot drift apart
STRAINS = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                      sep="\t", index_col=0).index.tolist()
GENE_DIR = "05.pan-genome/pathogenicity_genes"
OUT = "07.lethal_groups/literature_genes"
os.makedirs(OUT, exist_ok=True)
tmp = f"{OUT}/tmp"
os.makedirs(tmp, exist_ok=True)

genes = sorted(f[:-6] for f in os.listdir(GENE_DIR) if f.endswith(".fasta"))

# merged proteome with strain-prefixed headers
prot_all = f"{tmp}/all_strains.proteins.fa"
with open(prot_all, "w") as w:
    for s in STRAINS:
        p = f"02.annotation/per_isolates/{s}/funannotate/predict_results/fusarium.proteins.fa"
        for line in open(p):
            if line.startswith(">"):
                w.write(f">{s}|{line[1:]}")
            else:
                w.write(line)
print(f"merged proteome -> {prot_all}")

# mmseqs easy-search (query=literature genes, target=proteome)
q = f"{tmp}/lit_genes.fa"
with open(q, "w") as w:
    for g in genes:
        w.write(open(f"{GENE_DIR}/{g}.fasta").read())
db = f"{tmp}/prot_db"
subprocess.run(["mmseqs", "createdb", prot_all, db], check=True)
subprocess.run(["mmseqs", "createindex", db, f"{tmp}/idx", "--remove-tmp-files", "1"],
               check=True)
res = f"{OUT}/lit_vs_strains.tsv"
subprocess.run(["mmseqs", "easy-search", q, db, res, f"{tmp}/search_tmp",
                "--format-output", "query,target,fident,alnlen,mismatch,qcov,tcov,evalue,bits",
                "--threads", "8"], check=True)
hits = pd.read_csv(res, sep="\t", header=None,
                    names=["query", "target", "fident", "alnlen", "mismatch",
                           "qcov", "tcov", "evalue", "bits"])
# accession -> gene name (normalise sp|ACC|NAME to ACC)
acc2gene = {}
for g in genes:
    for line in open(f"{GENE_DIR}/{g}.fasta"):
        if line.startswith(">"):
            hdr = line[1:].split()[0]
            parts = hdr.split("|")
            acc = parts[1] if len(parts) >= 3 and parts[0] == "sp" else parts[0]
            acc2gene[acc] = g
hits["gene"] = hits["query"].map(acc2gene)
hits["strain"] = hits.target.str.split("|").str[0]
print("hits:", len(hits), "| genes hit:", hits.gene.nunique())

# copy number = hits per gene x strain (fident>=30, qcov>=40)
hits = hits[(hits.fident >= 0.3) & (hits.qcov >= 0.4)]
copy = hits.pivot_table(index="strain", columns="gene", values="target",
                        aggfunc="count", fill_value=0)
copy = copy.reindex(index=STRAINS, columns=genes, fill_value=0)
ident = hits.groupby(["strain", "gene"]).fident.max().unstack(fill_value=0.0)
ident = ident.reindex(index=STRAINS, columns=genes, fill_value=0.0)
copy.to_csv(f"{OUT}/gene_strain_copy.tsv", sep="\t")
ident.to_csv(f"{OUT}/gene_strain_identity.tsv", sep="\t")
print("\ncopy number matrix (nonzero fraction):")
print((copy > 0).mean().round(2))

# ---- association with phenotype ----
pheno = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)
# Same universe as the orthogroup screen: 1 - AUC at 10^3 and 10^4 plus the
# two-dose means and their difference. The 10^5 dose is excluded because 16 of
# the 21 isolates sit on the mortality ceiling at both temperatures.
AUC = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics_auc_no10e5.tsv",
                  sep="\t", index_col=0)
conds = ["vir_28Ce3", "vir_28Ce4", "vir_37Ce3", "vir_37Ce4",
         "vir_28_mean2", "vir_37_mean2", "vir_temp_effect2"]
rows = []
for g in genes:
    for c in conds:
        y = AUC[c].reindex(copy.index)
        x = copy[g]
        if x.nunique() <= 1:
            continue
        rho, p = stats.spearmanr(x, y)
        rows.append({"gene": g, "cond": c, "n_strains_with": int((x > 0).sum()),
                     "spearman_rho": rho, "p": p})
assoc = pd.DataFrame(rows)
# FDR across all tests
def bh_fdr(pvals):
    p = np.asarray(pvals)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum(1, p[order] * n / (np.arange(n) + 1))
    q[order] = np.minimum.accumulate(q[order][::-1])[::-1]
    return q
assoc["q"] = bh_fdr(assoc["p"].values)
assoc = assoc.sort_values("p")
assoc.to_csv(f"{OUT}/gene_phenotype_assoc.tsv", sep="\t", index=False)
sig = assoc[assoc.q < 0.1]
print(f"\nassociations (Spearman) tested: {len(assoc)}; q<0.1: {len(sig)}")
print(sig[["gene", "cond", "n_strains_with", "spearman_rho", "p", "q"]]
      .round(3).to_string(index=False))

# ---- figure ----
fig, axes = plt.subplots(2, 2, figsize=(10, 6), dpi=300)
# a: heatmap copy number
ax = axes[0, 0]
im = ax.imshow(copy.T.values, cmap="YlOrRd", aspect="auto")
ax.set_xticks(range(len(copy))); ax.set_xticklabels(copy.index, rotation=90, fontsize=5)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
ax.set_title("a  Copy number of literature pathogenicity genes")
fig.colorbar(im, ax=ax, fraction=0.046)
# b: top association scatter example
ax = axes[0, 1]
if len(sig):
    top = sig.iloc[0]
    g, c = top.gene, top.cond
    x = copy[g]; y = AUC[c].reindex(copy.index)
    ax.scatter(x, y, s=30, c="#3B7DD8")
    for s_ in copy.index:
        ax.annotate(s_.replace("BMU", ""), (x[s_], y[s_]), fontsize=5)
    ax.set_xlabel(f"{g} copies"); ax.set_ylabel(f"1 - AUC ({c})")
    ax.set_title(f"b  {g} vs {c} (rho={top.spearman_rho:.2f}, p={top.p:.3g})")
# c: |rho| heatmap gene x condition
ax = axes[1, 0]
hm = assoc.pivot_table(index="gene", columns="cond", values="spearman_rho")
hm = hm.reindex(genes)
im = ax.imshow(hm.values, cmap="RdBu_r", aspect="auto", vmin=-1, vmax=1)
ax.set_xticks(range(len(hm.columns))); ax.set_xticklabels(hm.columns, fontsize=6, rotation=45, ha="right")
ax.set_yticks(range(len(hm.index))); ax.set_yticklabels(hm.index, fontsize=6)
ax.set_title("c  Spearman rho: gene copy vs mortality per condition")
fig.colorbar(im, ax=ax, fraction=0.046)
# d: identity heatmap
ax = axes[1, 1]
im = ax.imshow(ident.T.values * 100, cmap="viridis", aspect="auto", vmin=0, vmax=100)
ax.set_xticks(range(len(ident))); ax.set_xticklabels(ident.index, rotation=90, fontsize=5)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
ax.set_title("d  Best-hit identity (%)")
fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_literature_genes.pdf")
fig.savefig(f"{OUT}/figure_literature_genes.png")
print("\ndone ->", OUT)
