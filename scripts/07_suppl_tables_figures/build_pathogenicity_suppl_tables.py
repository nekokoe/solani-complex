#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_pathogenicity_suppl_tables.py — the three supplementary tables the
21-strain re-analysis needs to be reproducible.

  S11  Galleria survival raw data  (21 isolates x 6 conditions x day 0-10)
  S12  cell-assay raw data         (2 models x [4 cytokines mRNA + 2 protein] x 2 doses x 3 reps)
  S13  the 71 publicly available genomes used, with their accessions and QC

Existing supplementary tables are numbered up to S10 in the text, so these take
S11-S13. Figures and tables run as separate series, so the S11 figure does not
collide with table S11.

Usage: python3 scripts/build_pathogenicity_suppl_tables.py
"""
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
SUPP = "00.manuscript/supplementary"
META = "05.pan-genome/bmu_ref_groups.211.txt"


def main():
    # ---------------- S11: Galleria survival raw ----------------
    raw = pd.read_csv("07.lethal_groups/caterpillar_37C_narrow.txt", sep="\t")
    phen = "07.lethal_groups/phenotype_analysis"
    mort = pd.read_csv(f"{phen}/phenotype_mortality.tsv", sep="\t", index_col=0)
    auc = pd.read_csv(f"{phen}/phenotype_auc.tsv", sep="\t", index_col=0)
    lt = pd.read_csv(f"{phen}/phenotype_lt50.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(META, sep="\t")
    strains = sorted(raw.Sample.unique())
    smeta = (meta[meta.BMU_ID.isin(strains)]
             .rename(columns={"BMU_ID": "Isolate"})
             [["Isolate", "MedEnv", "Host", "Species"]]
             .sort_values("Isolate"))
    with pd.ExcelWriter(f"{SUPP}/suppl_table.S11_galleria_survival.xlsx") as w:
        smeta.to_excel(w, sheet_name="isolates", index=False)
        raw.to_excel(w, sheet_name="survival_long", index=False)
        for sheet, df in (("day10_mortality", mort), ("survival_AUC", auc),
                          ("LT50_days", lt)):
            df.round(6).rename_axis("Isolate").to_excel(w, sheet_name=sheet)
    print(f"S11 galleria: {len(raw)} survival rows, {len(smeta)} isolates, "
          f"mortality/auc/LT50 {mort.shape}")

    # ---------------- S12: cell assay raw ----------------
    cell = "07.lethal_groups/cell_assay"
    long = pd.read_csv(f"{cell}/cell_assay_long.tsv", sep="\t")
    summ = pd.read_csv(f"{cell}/cell_assay_summary.tsv", sep="\t")
    cvenv = pd.read_csv(f"{cell}/cell_clinical_vs_env.tsv", sep="\t")
    cphen = pd.read_csv(f"{cell}/cell_vs_phenotype.tsv", sep="\t")
    with pd.ExcelWriter(f"{SUPP}/suppl_table.S12_cell_cytokine.xlsx") as w:
        long.to_excel(w, sheet_name="raw_long", index=False)
        summ.round(6).to_excel(w, sheet_name="fold_change", index=False)
        cvenv.round(6).to_excel(w, sheet_name="clinical_vs_environmental", index=False)
        cphen.round(6).to_excel(w, sheet_name="vs_galleria_phenotype", index=False)
    print(f"S12 cell: {len(long)} raw rows, {len(summ)} strain readouts, "
          f"{len(cvenv)} source tests, {len(cphen)} phenotype correlations")

    # ---------------- S13: the 71 public genomes ----------------
    src = f"{SUPP}/Fusarium_solani_species_complex.GCA_SRA.xlsx"
    qc = pd.read_excel(src, sheet_name="QC-final")
    mt = pd.read_excel(src, sheet_name="Meta")
    pub = qc.merge(mt, on="GCA", how="left")
    pub = pub.rename(columns={
        "ShortName": "Isolate",
        "Rawdata\nAvaliable?": "Raw data available",
        "Native BUSCO\nCompleteness": "Native BUSCO completeness",
        "Refreshed BUSCO\nCompleteness": "Re-assembled BUSCO completeness",
        "Re-identified \nSpecies": "Re-identified species"})
    keep = ["Isolate", "GCA", "Re-identified species", "Complex", "Raw data available",
            "Native BUSCO completeness", "Re-assembled BUSCO completeness",
            "Biosample", "SRA", "Host", "Location", "URL"]
    pub = pub[[c for c in keep if c in pub.columns]]
    with pd.ExcelWriter(f"{SUPP}/suppl_table.S13_public_genomes.xlsx") as w:
        pub.to_excel(w, sheet_name="public_genomes", index=False)
    n_meta = int(pub.Biosample.notna().sum()) if "Biosample" in pub else 0
    print(f"S13 public genomes: {len(pub)} rows "
          f"({n_meta} with Biosample/SRA metadata joined from the Meta sheet); "
          f"species: {pub['Re-identified species'].value_counts().to_dict()}")
    for f in ("S11_galleria_survival", "S12_cell_cytokine", "S13_public_genomes"):
        p = f"{SUPP}/suppl_table.{f}.xlsx"
        print(f"  wrote {os.path.relpath(p, ROOT)}  ({os.path.getsize(p)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
