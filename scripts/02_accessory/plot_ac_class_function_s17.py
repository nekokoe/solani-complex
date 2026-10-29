#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate S17 figure (ac_class_function) with new naming:
transfer classes = bulk (T1) / scatter (T2)."""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_BULK = "#C44E52"   # red = bulk (Type-1)
C_SCAT = "#3B7DD8"   # blue = scatter (Type-2)

f = pd.read_csv("06.accessory-chrs/gene_transfer_211/suppl_table.ac_class_function.tsv", sep="\t")
# only significant rows, ordered by test_frac desc (matches original layout)
f = f[f.FDR < 0.05].sort_values("test_frac", ascending=False)

cats = f.Category.tolist()
x = range(len(cats))
fig, ax = plt.subplots(figsize=(6.4, 2.9), dpi=300)

w = 0.36
ax.bar([i - w/2 for i in x], f.test_frac * 100, w, color=C_BULK, label="bulk AC genes")
ax.bar([i + w/2 for i in x], f.control_frac * 100, w, color=C_SCAT, alpha=0.85, label="scatter AC genes")

for i, r in f.iterrows():
    if r.FDR < 0.001:
        ax.text(i, r.test_frac * 100 + 0.15, "***", ha="center", fontsize=7)
    elif r.FDR < 0.01:
        ax.text(i, r.test_frac * 100 + 0.15, "**", ha="center", fontsize=7)
    else:
        ax.text(i, r.test_frac * 100 + 0.15, "*", ha="center", fontsize=7)

ax.set_xticks(list(x))
ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=6.5)
ax.set_ylabel("Fraction of transferred genes (%)", fontsize=7.5)
ax.set_title("Function of genes transferred by bulk- vs scatter-transfer ACs", fontsize=8.5)
ax.legend(fontsize=6.5, frameon=False, loc="upper right")
ax.tick_params(labelsize=6.5)
ax.set_ylim(0, max(f.test_frac) * 100 * 1.18)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("06.accessory-chrs/gene_transfer_211/figure_ac_class_function.pdf")
fig.savefig("06.accessory-chrs/gene_transfer_211/figure_ac_class_function.png")
print("S17 regenerated:", len(cats), "categories")
