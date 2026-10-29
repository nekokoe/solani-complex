#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 2 (AUC arm): literature pathogenicity genes vs survival-AUC phenotype.

Copy of `literature_genes_phenotype.py` with exactly one change: the phenotype
source. Instead of the day-10 mortality rates in `phenotype_metrics.tsv`, this
run uses the survival-curve AUC-derived traits

    vir_<cond> = 1 - auc_<cond>        (cond in 28Ce3 ... 37Ce5)

read from `07.lethal_groups/phenotype_analysis/phenotype_auc.tsv`. Higher
`vir` = more lethal, i.e. the same direction as mortality, so rho signs stay
comparable with the mortality arm.

Everything else (gene set, 21 strains, copy-number definition, Spearman,
BH-FDR) is unchanged, so the two arms are directly comparable.

The per-gene copy-number matrix is reused from `gene_strain_copy.tsv` when it
already exists (it does, from the mortality run); the mmseqs search is only
re-run — into `*_auc` paths, never over the original outputs — if that matrix
is missing. No figure is produced here; this script only writes the table.

Outputs:
  - 07.lethal_groups/literature_genes/gene_phenotype_assoc_auc.tsv
"""

import os, subprocess
import numpy as np
import pandas as pd
from scipy import stats

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

copy_file = f"{OUT}/gene_strain_copy.tsv"
ident_file = f"{OUT}/gene_strain_identity.tsv"

if os.path.exists(copy_file):
    print(f"reusing copy-number matrix {copy_file}")
    copy = pd.read_csv(copy_file, sep="\t", index_col=0)
    copy = copy.reindex(index=STRAINS, columns=genes, fill_value=0)
    ident = (pd.read_csv(ident_file, sep="\t", index_col=0)
             if os.path.exists(ident_file) else None)
else:
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
    res = f"{OUT}/lit_vs_strains_auc.tsv"
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
    copy.to_csv(f"{OUT}/gene_strain_copy_auc.tsv", sep="\t")
    ident.to_csv(f"{OUT}/gene_strain_identity_auc.tsv", sep="\t")

print("copy number matrix (nonzero fraction):")
print((copy > 0).mean().round(2))

# ---- association with the AUC-derived phenotype (THE ONLY CHANGE) ----
auc = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_auc.tsv",
                  sep="\t", index_col=0)
conds = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
pheno = pd.DataFrame({c: 1.0 - auc[f"auc_{c}"] for c in conds})  # vir_<cond>
pheno = pheno.reindex(copy.index)                                # 1 - AUC, 0..1
print("\nAUC-derived phenotypes (vir = 1 - AUC):")
print(pheno.round(3).to_string())
rows = []
for g in genes:
    for c in conds:
        y = pheno[c].reindex(copy.index)
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
out_tsv = f"{OUT}/gene_phenotype_assoc_auc.tsv"
assoc.to_csv(out_tsv, sep="\t", index=False)
sig = assoc[assoc.q < 0.1]
print(f"\nAUC-arm associations (Spearman) tested: {len(assoc)}; q<0.1: {len(sig)}")
print(f"min q: {assoc['q'].min():.4g}")
print(assoc.head(5)[["gene", "cond", "n_strains_with", "spearman_rho", "p", "q"]]
      .round(4).to_string(index=False))
print("\ndone ->", out_tsv)
