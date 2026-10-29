#!/usr/bin/env python3
"""Presence interval of every orthogroup computed WITHIN a single species.

The complex-wide interval table (`all/N0.HOG.hog_count.interval.tsv`) scores presence over all
211 genomes, so a gene that is fixed within a species that is represented by a handful of
genomes still falls in a low bin. Once the pan-genome section moves to a within-species analysis,
each species needs its own interval assignment, on its own genome set.

Do not substitute `05.pan-genome/accessory_genes/HOG_list.<sp>.interval_*.txt`: against the
recomputed values they agree on only about three quarters of orthogroups (881 orthogroups in the
*F. solani* 80-100% list have a within-species presence below 20%) and no script that produced
them survives on disk.

Writes 05.stats_pan_by_species_shorter/per_species/{solani,falciformis,keratoplastica}.HOG.interval.tsv
"""

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "05.stats_pan_by_species_shorter/per_species")

PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
SPECIES = os.path.join(ROOT, "05.stats_pan_by_species_shorter/bmu_ref_groups.species.txt")
GENOME_LISTS = os.path.join(ROOT, "05.pan-genome/accessory_genes")

# Presence_Rate in the reference table is a percentage and the bins are left-closed
# ([0,20), [20,40), ... [80,inf)). The upper edge must be open: 4,321 orthogroups sit at
# exactly 100% and a closed [80,100) bin would drop every one of them.
EDGES = [0, 20, 40, 60, 80, np.inf]
LABELS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]


def main():
    os.makedirs(OUT, exist_ok=True)
    pav = pd.read_csv(PAV, sep="\t")
    species = pd.read_csv(SPECIES, sep="\t", header=None, names=["SampleID", "species"])

    for sp in ("solani", "falciformis", "keratoplastica"):
        listfile = os.path.join(GENOME_LISTS, f"all_{sp}.id.shorter.txt")
        genomes = [x.strip() for x in open(listfile) if x.strip()]
        n_declared = int((species.species == sp).sum())
        assert len(genomes) == n_declared, (
            f"{sp}: {len(genomes)} genomes in the list vs {n_declared} in the species file")

        sub = pav[pav.SampleID.isin(genomes)]
        rate = (sub.groupby("HOG_ID").Present.mean() * 100).round(6)
        tab = pd.DataFrame({"HOG_ID": rate.index, "Presence_Rate": rate.values})
        tab["Presence_Interval"] = pd.cut(tab.Presence_Rate, EDGES, right=False,
                                          labels=LABELS).astype(str)

        path = os.path.join(OUT, f"{sp}.HOG.interval.tsv")
        tab.to_csv(path, sep="\t", index=False)
        counts = tab.Presence_Interval.value_counts().reindex(LABELS)
        print(f"{sp:16s} genomes={len(genomes):3d}  HOGs={len(tab):,}  "
              f"bins={counts.to_dict()}")
        assert len(tab) == pav.HOG_ID.nunique(), f"{sp}: HOG count mismatch"


if __name__ == "__main__":
    main()
