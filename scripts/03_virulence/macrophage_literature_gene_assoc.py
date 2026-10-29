#!/usr/bin/env python3
"""Macrophage readouts x literature virulence genes, Spearman + Benjamini-Hochberg.

The 21-strain macrophage data (THP-1 and PBMC, four qPCR markers and two ELISA
proteins, each at two concentrations) give 24 readout conditions. The 16
literature-reported virulence genes whose copy number varies across the 21 strains
are correlated against each readout, giving gene x readout tests; q values are
Benjamini-Hochberg across all tests, as in the Galleria literature-gene screen.

Writes 07.lethal_groups/cell_assay/cell_literature_gene_assoc.tsv
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CELL = os.path.join(ROOT, "07.lethal_groups/cell_assay/cell_assay_summary.tsv")
GENES = os.path.join(ROOT, "07.lethal_groups/literature_genes/gene_strain_copy.tsv")
OUT = os.path.join(ROOT, "07.lethal_groups/cell_assay/cell_literature_gene_assoc.tsv")


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
    cell = pd.read_csv(CELL, sep="\t")
    genes = pd.read_csv(GENES, sep="\t")
    gene_cols = [c for c in genes.columns if c != "strain"]
    # keep only genes whose copy number varies across the 21 strains
    variable = [g for g in gene_cols if genes[g].std() > 0]
    rows = []
    for (model, cyto, assay, conc), sub in cell.groupby(["model", "cytokine", "assay", "concentration"]):
        sub = sub.sort_values("strain")
        y = sub.set_index("strain")["fold_vs_ctrl"]
        for g in variable:
            x = genes.set_index("strain")[g]
            both = x.index.intersection(y.index)
            if len(both) < 10:
                continue
            rho, p = stats.spearmanr(x[both], y[both])
            rows.append(dict(model=model, cytokine=cyto, assay=assay, concentration=int(conc),
                             gene=g, n=len(both), spearman_rho=round(float(rho), 4),
                             p=float(p)))
    res = pd.DataFrame(rows)
    res["q_bh"] = bh(res.p)
    res = res.sort_values(["q_bh", "p"])
    res.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {os.path.relpath(OUT, ROOT)} ({len(res)} tests, {len(variable)} variable genes)")
    print(f"best q = {res.q_bh.min():.4f} | tests with nominal p < 0.05: {(res.p < 0.05).sum()}")
    print(res.head(8)[["model", "cytokine", "assay", "concentration", "gene", "spearman_rho", "p", "q_bh"]].to_string(index=False))


if __name__ == "__main__":
    main()
