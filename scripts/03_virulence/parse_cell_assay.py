#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""parse_cell_assay.py — 21-strain cell-assay xlsx -> a tidy long table.

The 2026-10-12 file (07.lethal_groups/21株-细胞原始数据.xlsx) has one sheet per
cytokine. Each sheet stacks blocks:

    <label row>  高孢子浓度 / 低孢子浓度 / ELISA-高浓度 / ELISA-低浓度
    <header row> Ctrl, LPS, then the 21 BMU strains
    3 rows        relative mRNA expression (or relative protein level)

So: 4 cytokines x mRNA + 2 cytokines x ELISA protein, each with 2 concentrations
x 3 replicates x (21 strains + Ctrl + LPS).

Two things to be careful about:
  * Ctrl and LPS carry byte-identical values in the high and low blocks - they are
    ONE shared control pair, not a separate control per concentration, so they
    must not be counted twice (and never as experimental strains).
  * The concentrations are 10^5 (high) and 10^3 (low) conidia.

Output: 07.lethal_groups/cell_assay_long.tsv and cell_assay_summary.tsv
"""
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
XLSX = "07.lethal_groups/21株-细胞原始数据.xlsx"
OUT = "07.lethal_groups/cell_assay"
CONC = {"高孢子浓度": 1e5, "低孢子浓度": 1e3, "ELISA-高浓度": 1e5, "ELISA-低浓度": 1e3}
ASSAY = {"高孢子浓度": "mRNA", "低孢子浓度": "mRNA",
         "ELISA-高浓度": "protein", "ELISA-低浓度": "protein"}


def parse_sheet(df, cytokine):
    rows = []
    labels = [i for i in range(len(df))
              if str(df.iat[i, 0]).strip() in CONC]
    for li, i in enumerate(labels):
        label = str(df.iat[i, 0]).strip()
        header = [str(x).strip() for x in df.iloc[i + 1].tolist()]
        end = labels[li + 1] if li + 1 < len(labels) else len(df)
        block = df.iloc[i + 2:end]
        block = block[block.iloc[:, 0].notna()]
        for _, r in block.iterrows():
            # NB the cell sheets have NO leading label column: column 0 of a data
            # row is already the first control ("Ctrl"), so start at 0. Starting at
            # 1 silently dropped the entire Ctrl reference.
            for c in range(len(header)):
                name = header[c]
                if not name or name == "nan":
                    continue
                v = r.iloc[c]
                if pd.isna(v):
                    continue
                rows.append({"cytokine": cytokine, "assay": ASSAY[label],
                             "concentration": CONC[label], "unit": name,
                             "value": float(v)})
    return rows


def main():
    xl = pd.ExcelFile(XLSX)
    allrows = []
    for sh in xl.sheet_names:
        d = pd.read_excel(XLSX, sheet_name=sh, header=None)
        allrows += parse_sheet(d, sh)
        print(f"{sh}: parsed")
    long = pd.DataFrame(allrows)
    long["is_control"] = long.unit.isin(["Ctrl", "LPS"])
    long["strain"] = np.where(long.is_control, NA := pd.NA, long.unit)
    os.makedirs(OUT, exist_ok=True)
    long.to_csv(f"{OUT}/cell_assay_long.tsv", sep="\t", index=False)

    # sanity: the control pair must be identical across concentrations
    ctrl = long[long.is_control]
    piv = (ctrl.groupby(["cytokine", "assay", "unit", "concentration"]).value.mean()
           .groupby(["cytokine", "assay", "unit"]).nunique())
    print(f"\nrows: {len(long)}  cytokines: {sorted(long.cytokine.unique())}")
    print(f"assays: {sorted(long.assay.unique())}  concentrations: "
          f"{sorted(long.concentration.unique())}")
    print(f"strains per cytokine/assay/conc: "
          f"{long[~long.is_control].groupby(['cytokine','assay','concentration']).strain.nunique().unique()}")
    print(f"control MEANS identical across the two concentrations: "
          f"{bool((piv == 1).all())}  (one shared control pair)")
    print(f"replicates per strain: "
          f"{sorted(long[~long.is_control].groupby(['cytokine','assay','concentration','strain']).size().unique())}")

    # fold change vs the Ctrl of the same cytokine/assay/concentration
    base = (ctrl[ctrl.unit == "Ctrl"]
            .groupby(["cytokine", "assay", "concentration"]).value.mean()
            .rename("ctrl_mean"))
    s = (long[~long.is_control]
         .groupby(["cytokine", "assay", "concentration", "strain"]).value
         .agg(["mean", "std", "size"]).reset_index()
         .merge(base.reset_index(), on=["cytokine", "assay", "concentration"]))
    s["fold_vs_ctrl"] = s["mean"] / s.ctrl_mean
    s.to_csv(f"{OUT}/cell_assay_summary.tsv", sep="\t", index=False)
    print(f"\nfold change vs Ctrl, by cytokine (10^5):")
    hi = s[s.concentration == 1e5]
    print(hi.groupby(["cytokine", "assay"]).fold_vs_ctrl.describe()[["mean", "min", "max"]]
          .round(2).to_string())
    print(f"\nwrote {OUT}/cell_assay_long.tsv and cell_assay_summary.tsv")


if __name__ == "__main__":
    main()
