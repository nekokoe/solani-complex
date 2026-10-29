#!/usr/bin/env python3
"""Supplementary Table S17 - population differentiation among and within species.

Collects the numbers behind the "Isolation source explains less variance than
subpopulation membership" subsection, none of which appears in a main figure except
the 11x11 FST matrix of Fig. 2G:

  amova               AMOVA variance partition between clinical and non-clinical
                      isolates (phi_ST, variance components, permutation P)
  source_fst          Weir-Cockerham FST between the two source groups, computed from
                      the vcftools 10-kb window output used in the main text
  subpopulation_fst   between-subpopulation FST, the comparison the text draws
  definitions         distance metric, permutations, software, per-species SNP counts

  python3 04.population/amova_suppl_table.py
"""

import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "00.manuscript", "supplementary",
                   "suppl_table.S17_population_differentiation.xlsx")
AMOVA = os.path.join(ROOT, "05.stats_pan_by_species_shorter", "amova_results.tsv")
MATRIX = os.path.join(HERE, "_fig2", "data", "panelJ_fst_matrix_v2.tsv")
SNP = {"solani": 454466, "falciformis": 828841, "keratoplastica": 566279}


def source_fst():
    rows = []
    for sp in ("falciformis", "solani", "keratoplastica"):
        path = glob.glob(os.path.join(HERE, f"solani-short-{sp}", f"fst.{sp}.med_vs_env.windowed.weir.fst"))
        if not path:
            continue
        df = pd.read_csv(path[0], sep="\t")
        df = df[(df.N_VARIANTS > 0)].dropna(subset=["MEAN_FST", "WEIGHTED_FST"])
        n = df.N_VARIANTS.to_numpy(float)
        rows.append(dict(
            species=sp,
            n_windows=int(len(df)),
            pct_windows_negative=round(float((df.MEAN_FST < 0).mean() * 100), 1),
            mean_of_window_MEAN_FST=round(float(df.MEAN_FST.mean()), 4),
            mean_of_window_WEIGHTED_FST=round(float(df.WEIGHTED_FST.mean()), 4),
            variant_weighted_FST=round(float(np.average(df.WEIGHTED_FST, weights=n)), 4),
        ))
    return pd.DataFrame(rows)


def subpopulation_fst(match):
    df = pd.read_csv(MATRIX, sep="\t")
    df = df[df.group_a.str.contains("-pop") & df.group_b.str.contains("-pop")]
    df = df[df.group_a.str.split("-").str[0] == df.group_b.str.split("-").str[0]]
    out = df.groupby(df.group_a.str.split("-").str[0]).agg(
        n_pairs=("fst", "size"), fst_min=("fst", "min"), fst_max=("fst", "max"),
        p_max=("p_perm", "max")).reset_index().rename(columns={"group_a": "species"})
    for c in ("fst_min", "fst_max", "p_max"):
        out[c] = out[c].round(4)
    return out


def main():
    amova = pd.read_csv(AMOVA, sep="\t")
    amova["phi_ST"] = amova.phi_ST.round(4)
    amova[["sigma2_source", "sigma2_error"]] = amova[["sigma2_source", "sigma2_error"]].round(5)
    amova["n_snps"] = amova.species.map(SNP)
    sf = source_fst()
    sub = subpopulation_fst(MATRIX)

    definitions = pd.DataFrame({
        "item": ["FST", "FST windows", "FST summary", "genes/SNPs", "AMOVA", "AMOVA input",
                 "permutations", "software", "why two statistics", "note"],
        "value": [
            "Weir & Cockerham (1984) estimator, biallelic core SNPs",
            "10-kb windows, 5-kb step (vcftools --fst-window-size 10000 --fst-window-step 5000)",
            "the main text quotes the mean over windows of the per-window MEAN_FST column, the same "
            "summary used for the between-subpopulation matrix of Fig. 2G",
            "; ".join(f"{k}: {v:,}" for k, v in SNP.items()),
            "analysis of molecular variance, phi-statistics (Excoffier et al. 1992)",
            "simple-matching distance over non-missing core-SNP genotypes (1 - matching proportion)",
            "999 permutations (AMOVA); the FST matrix of Fig. 2G uses individual-level permutations",
            "vcftools v0.1.16 for FST and pi; pegas (R) for AMOVA",
            "phi_ST is distance-based and not numerically comparable with FST, so the two are reported "
            "separately: FST for the magnitude of the source-versus-subpopulation partition, AMOVA for "
            "the significance of the source partition",
            "negative FST and phi_ST estimates arise when the true value is zero, because both are ratios "
            "of estimated variance components; they are read as no detectable differentiation",
        ],
    })

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        amova.to_excel(xw, sheet_name="amova", index=False)
        sf.to_excel(xw, sheet_name="source_fst", index=False)
        sub.to_excel(xw, sheet_name="subpopulation_fst", index=False)
        definitions.to_excel(xw, sheet_name="definitions", index=False)

    print(amova.to_string(index=False))
    print()
    print(sf.to_string(index=False))
    print()
    print(sub.to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
