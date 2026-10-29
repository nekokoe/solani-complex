#!/usr/bin/env python3
"""Rebuild the S27 plot intermediates from the current toolkit HOG sets.

Why this script exists
    `07.lethal_groups/genomewide_assoc/tmp/s27_{long,meta,ann}.tsv` are read by
    `plot_hog_bias_s27v4.R` but were not written by any script in the repo (they were
    produced ad hoc on 2026-10-12 for the 15/1 toolkit). Changing the toolkit therefore
    silently leaves S27 stale. This is the missing writer.

Rules are the ones `derive_virulence_modes.py` documents and `plot_hog_bias_s27v2.py`
implements: a HOG is chitin-related if its Pfam domains or gene name match CHIT, and
infection-relevant if they match INF; `tag` is the first gene name, else the first Pfam
domain.

Output (overwrites, with one-off backups on first run)
    tmp/s27_long.tsv   HOG, strain, copy
    tmp/s27_meta.tsv   HOG, block (37C/28C), chit, inf, tag
    tmp/s27_ann.tsv    HOG, Gene_Name, Pfam_Domain
"""

import os
import shutil

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GW = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc")
TMP = os.path.join(GW, "tmp")
PHENO = os.path.join(ROOT, "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv")
HOGCOUNT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
ANNOT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")

CHIT = ["Glyco_hydro_18", "Glyco_hydro_19", "Glyco_hydro_20", "Glyco_hydro_12",
        "Chitin_synth", "chitin", "Glyco_hydro_28"]
INF = ["Hce2", "CFEM", "LysM", "HET", "Peptidase_S", "Effector"]


def hog_ids(path):
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path, sep="\t", usecols=[0])
    return [h for h in df.iloc[:, 0] if isinstance(h, str) and h.startswith("N0.")]


def get(annot, h, col):
    row = annot[annot.HOG_ID == h]
    if row.empty:
        return ""
    v = str(row[col].iloc[0])
    return "" if v.lower() in ("nan", "-", "none") else v


def main():
    pheno = pd.read_csv(PHENO, sep="\t", index_col=0)
    strains = list(pheno.index)

    sp37 = hog_ids(os.path.join(GW, "hog_sp37.tsv"))
    sp28 = hog_ids(os.path.join(GW, "hog_sp28.tsv"))
    print(f"sp37 = {len(sp37)} HOGs · sp28 = {len(sp28)} HOGs")
    assert sp37 and sp28, "run scripts/build_toolkit_hog_sets.py first"

    annot = pd.read_csv(ANNOT, sep="\t", dtype=str)
    rows_meta, rows_ann = [], []
    for block, ids in (("37C", sp37), ("28C", sp28)):
        for h in ids:
            pf = get(annot, h, "Pfam_Domain")
            gn = get(annot, h, "Gene_Name")
            pfs = [p for p in pf.split(";") if p]
            chit = any(c in pf for c in CHIT) or "chitin" in gn.lower()
            inf = any(c in pf for c in INF)
            tag = (gn.split(";")[0] if gn else "") or (pfs[0] if pfs else "")
            rows_meta.append(dict(HOG=h, block=block, chit=chit, inf=inf, tag=tag))
            rows_ann.append(dict(HOG=h, Gene_Name=gn, Pfam_Domain=pf))

    meta = pd.DataFrame(rows_meta)
    # meta is read in row order as the display order: 37C block first, then 28C
    long_rows = []
    count = pd.read_csv(HOGCOUNT, sep="\t")
    mat = count[count.iloc[:, 0].isin(strains)].pivot(
        index="HOG_ID", columns="SampleID", values="GeneCount").fillna(0)
    for h in meta.HOG:
        for s in strains:
            long_rows.append(dict(HOG=h, strain=s,
                                  copy=int(mat.at[h, s]) if h in mat.index else 0))
    long = pd.DataFrame(long_rows)

    os.makedirs(TMP, exist_ok=True)
    for name, df in (("s27_long.tsv", long), ("s27_meta.tsv", meta),
                     ("s27_ann.tsv", pd.DataFrame(rows_ann))):
        path = os.path.join(TMP, name)
        if os.path.exists(path) and not os.path.exists(path + ".bak-before-cmh"):
            shutil.copy2(path, path + ".bak-before-cmh")
        df.to_csv(path, sep="\t", index=False)
        print(f"wrote {path}  ({len(df)} rows)")

    print(f"\nblocks: {meta.block.value_counts().to_dict()}")
    print(f"chitin HOGs: {int(meta.chit.sum())} · infection-relevant: {int(meta.inf.sum())}")


if __name__ == "__main__":
    main()
