#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the sidC phylogram with ETE (run under the braker python 3.9 env).
Output: 07.lethal_groups/literature_genes/tmp/tree_d.png
"""

import pandas as pd
from ete3 import Tree, TreeStyle, TextFace, CircleFace
from matplotlib import cm, colors

BASE = "/share/pasteur/gaosh/Fusarium/phase2/solaniPathogenicityPaper"
pheno = pd.read_csv(f"{BASE}/07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                    sep="\t", index_col=0)

t = Tree(f"{BASE}/07.lethal_groups/literature_genes/tmp/sidC_tree_midroot.nwk")

cmap = cm.get_cmap("YlOrRd")
norm = colors.Normalize(0, 1)

ts = TreeStyle()
ts.show_branch_length = True
ts.show_scale = True
ts.scale = 800            # px per branch-length unit
ts.branch_vertical_margin = 14
ts.title.add_face(TextFace("", fsize=1), column=0)  # title rendered by matplotlib
ts.margin_left = 4
ts.margin_right = 4
ts.margin_top = 4
ts.margin_bottom = 4

# branch-length labels on every branch
for node in t.traverse():
    if not node.is_root():
        node.add_face(TextFace(f" {node.dist:.4f}", fsize=9, fgcolor="grey40"),
                      column=0, position="branch-top")

for leaf in t.iter_leaves():
    name = leaf.name
    mort = pheno.loc[name, "mort_28Ce4"] if name in pheno.index else float("nan")
    col = colors.to_hex(cmap(norm(mort)))
    leaf.add_face(CircleFace(9, col), column=0, position="branch-right")
    leaf.add_face(TextFace(f"  {name[3:]}  {mort:.2f}", fsize=16),
                  column=1, position="branch-right")

# remove default leaf name display (we add our own)
ts.show_leaf_name = False

t.render(f"{BASE}/07.lethal_groups/literature_genes/tmp/tree_d.png",
         w=1000, h=1200, dpi=300, tree_style=ts)
print("tree_d.png rendered")
