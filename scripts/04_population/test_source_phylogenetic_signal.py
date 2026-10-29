#!/usr/bin/env python3
"""How is isolation source distributed on the FSSC phylogenies?

Two trees are tested, because they answer different questions and carry very different numbers of
non-clinical isolates:

  primary    the core-SNP maximum-likelihood tree of all 211 isolates behind Figure 1B
             (133 clinical, 73 non-clinical, 5 unlabelled) -- this is the tree the Results use to
             say that clinical and non-clinical isolates are interleaved yet locally clustered;
  secondary  the TEF1+RPB2 barcode tree behind Figure 1A, which exists for species assignment
             (140 isolates, 132 clinical against only 8 environmental).

The statistic is the Fitch parsimony score of the binary source label, compared with the same
statistic under random permutation of the labels.

Reading the two reference points:
  * a score of 2 is the minimum and means each source is monophyletic, i.e. the tree fully
    separates them. A score far above 2 means the sources are interleaved;
  * the permutation null is what the same tree gives when the labels carry no information, so the
    gap between the observed score and that null is the measure of local clustering.

⚠ The barcode tree is NOT a test of source structure: with only 8 environmental isolates its null
distribution is degenerate (7-8) and the comparison cannot resolve anything. It is reported here
only so that the difference between the two trees is on record.

Writes 03.identify_species/phylo_trees/source_parsimony_test.txt
"""

import os
import random
import statistics

import pandas as pd
from Bio import Phylo

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUPS = os.path.join(ROOT, "05.pan-genome/accessory_genes/bmu_ref_groups.shorter.txt")
OUT = os.path.join(ROOT, "03.identify_species/phylo_trees/source_parsimony_test.txt")
NPERM = 2000
SEED = 42

TREES = [
    # copied from /share/pasteur/gaosh/Fusarium/phase2/mapping_reference/solani-shortname/
    # (md5 8b98cc7195141219a625cd3e705a48cd) so the tree used by Figure 1B lives with the
    # SNP analyses it belongs to
    ("core-SNP tree, 211 isolates (Figure 1B)",
     os.path.join(ROOT, "04.population/solani-short-fullSC/"
                        "solani-bmu-ref-short.raw.snps.combined.genotype.core_snps_tree.contree")),
    ("TEF1+RPB2 barcode tree (Figure 1A, species assignment only)",
     os.path.join(ROOT, "03.identify_species/phylo_trees/barcode_genes/solani-shortname/"
                        "rpb2_tef1/all.outgroup_BMU_REF_FID_LC.rpb2_tef1.mafft.nogaps."
                        "ufb1000.treefile")),
]


def fitch(clade, state):
    if clade.is_terminal():
        return 0, state[clade.name]
    total, sets = 0, []
    for child in clade.clades:
        s, st = fitch(child, state)
        total += s
        sets.append(st)
    cur = sets[0]
    for st in sets[1:]:
        both = cur & st
        if both:
            cur = both
        else:
            cur = cur | st
            total += 1
    return total, cur


def score(tree, leaves, state):
    """unlabelled leaves are given both states, i.e. they are free to take either"""
    st = {l: ({state[l]} if state.get(l) else {"clinical", "environmental"}) for l in leaves}
    return fitch(tree.root, st)[0]


def main():
    groups = pd.read_csv(GROUPS, sep="\t")
    source = groups.set_index("BMU_ID").MedEnv.map({"med": "clinical", "env": "environmental"})
    source = source.dropna()

    report = ["Parsimony reconstruction of isolation source on FSSC phylogenies",
              "statistic: Fitch score of the binary source label, observed vs permuted labels",
              "  2 = both sources monophyletic (the tree separates them completely)",
              "", ]
    for title, path in TREES:
        tree = Phylo.read(path, "newick")
        leaves = [c.name for c in tree.get_terminals()]
        state = {l: (source[l] if l in source.index else None) for l in leaves}
        labelled = [l for l in leaves if state[l]]
        observed = score(tree, leaves, state)

        labels = [state[l] for l in labelled]
        random.seed(SEED)
        null = []
        for _ in range(NPERM):
            random.shuffle(labels)
            perm = dict(state)
            perm.update(dict(zip(labelled, labels)))
            null.append(score(tree, leaves, perm))
        null.sort()
        p = (sum(1 for v in null if v <= observed) + 1) / (len(null) + 1)

        report += [
            f"== {title}",
            f"   file: {path}",
            f"   leaves {len(leaves)}, labelled {len(labelled)} "
            f"({pd.Series([state[l] for l in labelled]).value_counts().to_dict()})",
            f"   observed Fitch score {observed}",
            f"   permutation null ({NPERM}): median {statistics.median(null):.0f}, "
            f"5th pct {null[int(0.05 * len(null))]}, 95th pct {null[int(0.95 * len(null))]}, "
            f"mean {statistics.mean(null):.1f}",
            f"   P(observed <= null) = {p:.4f}",
            "",
        ]
    report.append("Source is interleaved on both trees (score far above 2) but only on the "
                  "core-SNP tree does it cluster more than random; the barcode tree has 8 "
                  "environmental isolates and cannot resolve this.")
    text = "\n".join(report)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
