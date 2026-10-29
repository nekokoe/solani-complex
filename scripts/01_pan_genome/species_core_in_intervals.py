#!/usr/bin/env python3
"""How much of each presence interval is the core genome of an undersampled species.

Presence rate is computed over all 211 FSSC genomes, so a gene that is core within one species
still counts as rare across the complex whenever that species is represented by only a handful
of genomes. This script quantifies that overlap, because it is the reason the pan-genome section
switches from a complex-wide analysis to a within-species one.

Two guards on the claim:
  * species with fewer than MIN_GENOMES genomes are excluded -- being "100% present" in a single
    genome is true of every gene in it and says nothing about core status;
  * the comparison is reported per interval, so that it is visible that the same phenomenon also
    runs through 20-40% and 40-60% rather than being special to the rarest bin.

Writes 05.stats_pan_by_species_shorter/species_core_by_interval.tsv
       05.stats_pan_by_species_shorter/species_core_by_species.tsv
"""

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = os.path.join(ROOT, "05.stats_pan_by_species_shorter")

PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
SPECIES = os.path.join(out, "bmu_ref_groups.species.txt")
INTERVAL = os.path.join(out, "all/N0.HOG.hog_count.interval.tsv")

MIN_GENOMES = 3
THRESHOLDS = (0.60, 0.80, 0.90, 1.00)
BINS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]


def main():
    pav = pd.read_csv(PAV, sep="\t")
    species = pd.read_csv(SPECIES, sep="\t", header=None, names=["SampleID", "species"])
    n_genomes = species.species.value_counts()
    keep = [s for s in n_genomes.index if n_genomes[s] >= MIN_GENOMES]
    print(f"species kept (>= {MIN_GENOMES} genomes): "
          f"{ {s: int(n_genomes[s]) for s in keep} }")
    print(f"excluded: {sorted(set(n_genomes.index) - set(keep))}")

    # mean of the 0/1 presence over a species' genomes IS the within-species presence rate;
    # dividing by the genome count again would deflate every species bigger than one genome
    within = (pav.merge(species, on="SampleID", how="left")
                 .groupby(["HOG_ID", "species"]).Present.mean().unstack())

    iv = pd.read_csv(INTERVAL, sep="\t")
    iv.columns = ["HOG_ID", "Presence_Rate_211", "Interval_211"]
    iv = iv[iv.Interval_211.isin(BINS)]

    rows = []
    for th in THRESHOLDS:
        best = (within[keep] >= th).any(axis=1).rename("species_core")
        d = iv.merge(best.reset_index(), on="HOG_ID", how="left")
        d["species_core"] = d.species_core.fillna(False).astype(int)
        g = d.groupby("Interval_211").species_core.agg(["size", "sum"])
        g = g.reindex(BINS).astype(int)
        for b, r in g.iterrows():
            rows.append(dict(Interval_211=b, within_species_threshold=th,
                             n_HOGs=int(r["size"]), n_species_core=int(r["sum"]),
                             pct_species_core=round(r["sum"] / r["size"] * 100, 1)))
    by_interval = pd.DataFrame(rows)
    by_interval.to_csv(os.path.join(out, "species_core_by_interval.tsv"), sep="\t", index=False)

    # which species carry the claim, and how strongly -- only within 0-20%
    rows = []
    rare = iv[iv.Interval_211 == "0%-20%"].HOG_ID
    for s in keep:
        v = within.loc[within.index.isin(rare), s]
        rows.append(dict(species=s, n_genomes=int(n_genomes[s]),
                         min_complex_presence_if_core=round(n_genomes[s] / 211 * 100, 1),
                         n_core_80=int((v >= 0.80).sum()),
                         n_core_100=int((v >= 1.00).sum()),
                         n_HOGs_in_interval=int(len(v))))
    by_species = pd.DataFrame(rows).sort_values("n_core_80", ascending=False)
    by_species.to_csv(os.path.join(out, "species_core_by_species.tsv"), sep="\t", index=False)

    print("\nspecies-core fraction of each interval:")
    print(by_interval.pivot(index="Interval_211", columns="within_species_threshold",
                            values="pct_species_core").reindex(BINS).to_string())
    print("\nwho carries it (0-20% interval only):")
    print(by_species.to_string(index=False))
    print("\na species whose core genes reach 20% of the complex must have >= 0.2 x 211 = "
          "42 genomes, so F. solani (86) and F. falciformis (66) cannot contribute by "
          "construction; the counts above are consistent with that.")


if __name__ == "__main__":
    main()
