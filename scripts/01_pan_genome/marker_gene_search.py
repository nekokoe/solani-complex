#!/usr/bin/env python3
"""Is there any orthogroup that marks a clinical isolate? A genome-wide search.

The PCoA shows that accessory PAV separates clinical from non-clinical isolates only as a fuzzy
gradient. A gradient of that kind could still be carried by one gene that every clinical isolate
has and every non-clinical isolate lacks; this script tests that directly, over all 29,080
orthogroups and within each of the two analysable species.

Cohort and source labels are those of the random-forest screen
(09.rf_clinical_genes/rf_genes.per_genome.tsv): 60 clinical / 26 non-clinical for F. solani and
48 / 18 for F. falciformis, with the one unlabelled isolate counted as non-clinical, exactly as
the screen does.

For each orthogroup the presence rate is computed in the clinical and in the non-clinical
isolates of its species, and the orthogroups are then counted against a ladder of criteria,
from "present in every clinical isolate" down to "present in >= 80% of clinical and <= 20% of
non-clinical isolates".

Writes 08.pan_020_genes/marker_gene_search.tsv (the count ladder) and
        08.pan_020_genes/marker_gene_search_80_20.tsv (the orthogroups meeting the relaxed
        80% clinical / 20% non-clinical criterion, with their gene content; these are not
        conclusions but the natural shortlist for experimental follow-up).
"""

import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
PER_GENOME = os.path.join(ROOT, "09.rf_clinical_genes/rf_genes.per_genome.tsv")
OUT = os.path.join(ROOT, "08.pan_020_genes/marker_gene_search.tsv")

CRITERIA = [("100/0", 1.0, 0.0), ("90/10", 0.9, 0.1), ("80/20", 0.8, 0.2)]


def main():
    pg = pd.read_csv(PER_GENOME, sep="\t", usecols=["SampleID", "Source", "Species"])
    pg = pg[pg.Species.isin(["solani", "falciformis"])]
    pav = pd.read_csv(PAV, sep="\t")

    rows = []
    for sp in ("solani", "falciformis"):
        sub = pg[pg.Species == sp]
        clinical = sub[sub.Source == "med"].SampleID
        nonclin = sub[sub.Source != "med"].SampleID
        wide = (pav[pav.SampleID.isin(sub.SampleID)]
                .groupby(["HOG_ID", "SampleID"]).Present.mean().unstack())
        c = wide[clinical].mean(axis=1)
        e = wide[nonclin].mean(axis=1)
        rows.append(dict(species=sp, criterion="all orthogroups",
                         n_HOGs=len(wide),
                         note=f"clinical {len(clinical)} / non-clinical {len(nonclin)}"))
        rows.append(dict(species=sp, criterion="present in 100% of clinical isolates",
                         n_HOGs=int((c >= 0.999).sum())))
        for label, lo_c, hi_e in CRITERIA:
            rows.append(dict(species=sp,
                             criterion=f"present in >= {int(lo_c*100)}% of clinical AND "
                                       f"<= {int(hi_e*100)}% of non-clinical isolates",
                             n_HOGs=int(((c >= lo_c) & (e <= hi_e)).sum())))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT, sep="\t", index=False)
    print(tab.to_string(index=False))

    # ---- the relaxed 80/20 criterion, as a gene-level shortlist ----
    eggnog = pd.read_csv(os.path.join(ROOT, "02.annotation/all_isolates/all.fusarium.gene_eggnog.tsv"),
                         sep="\t", header=None, usecols=[0, 1],
                         names=["GeneID", "name"], dtype=str, low_memory=False)
    eggnog = eggnog.drop_duplicates("GeneID")
    gvh = pd.read_csv(os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv"),
                      sep="\t", usecols=["GeneID", "HOG_ID"])
    short = []
    for sp in ("solani", "falciformis"):
        sub = pg[pg.Species == sp]
        clinical = sub[sub.Source == "med"].SampleID
        nonclin = sub[sub.Source != "med"].SampleID
        wide = (pav[pav.SampleID.isin(sub.SampleID)]
                .groupby(["HOG_ID", "SampleID"]).Present.mean().unstack())
        c = wide[clinical].mean(axis=1)
        e = wide[nonclin].mean(axis=1)
        keep = (c >= 0.8) & (e <= 0.2)
        genes = (eggnog.merge(gvh, on="GeneID", how="inner")
                 .groupby("HOG_ID").name.agg(lambda v: ";".join(sorted({x for x in v if x != "-"}))))
        for hog in wide.index[keep]:
            short.append(dict(species=sp, HOG_ID=hog,
                              clinical_rate=round(float(c[hog]), 3),
                              nonclinical_rate=round(float(e[hog]), 3),
                              named_genes=genes.get(hog, "")))
    short_tab = pd.DataFrame(short).sort_values(["species", "clinical_rate", "nonclinical_rate"],
                                                ascending=[True, False, True])
    short_out = os.path.join(ROOT, "08.pan_020_genes/marker_gene_search_80_20.tsv")
    short_tab.to_csv(short_out, sep="\t", index=False)
    print(f"\n80/20 shortlist: {len(short_tab)} orthogroups "
          f"({short_tab.species.value_counts().to_dict()})")
    print(f"wrote {os.path.relpath(short_out, ROOT)}")


if __name__ == "__main__":
    main()
