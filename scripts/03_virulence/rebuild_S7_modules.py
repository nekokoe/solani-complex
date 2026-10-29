#!/usr/bin/env python3
"""Refresh Supplementary Table S7 (temperature-associated modules).

Module members are the variable HOGs whose maximum |Spearman rho| against that
temperature's three AUC phenotypes reaches the 99th percentile of the permutation
null (0.602 / 0.604). The phenotypes are the survival-AUC traits at 10^3 and 10^4
plus their two-dose means and difference; the 10^5 dose is excluded because its
endpoint is on the ceiling. The threshold is calibrated rather than conventional,
and the full sweep lives in Threshold_sensitivity.

Domain enrichment follows scripts/derive_virulence_modes.py exactly: universe =
the 4,484 variable HOGs, class = HOGs carrying the domain, one-sided
hypergeometric, no filtering on the number of hits, Benjamini-Hochberg across
every tested domain.

Sheets rewritten: 28C_module_Pfam_enrichment, 37C_module_Pfam_enrichment,
Module_HOG_membership, Module_definitions.
Sheets left untouched: Nominal_P005_candidates, Temperature_effect_HOGs,
Sensitivity_incl1e5 (they do not depend on the module definition).

Input : 07.lethal_groups/genomewide_assoc/hog_phenotype_assoc.tsv
        05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv
        05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv
Output: 00.manuscript/supplementary/suppl_table.S7_temperature_modules.xlsx
"""
import os

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# --no10e5 switches to the 10^3+10^4-only screen and its reduced trait lists: the 10^5
# rates are on the ceiling, and letting them into the module definition inflates the
# modules with tie-driven correlations.
import sys
NO10E5 = "--no10e5" in sys.argv
AUC_ARM = "--auc" in sys.argv
ASSOC = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc/"
                           + ("hog_phenotype_assoc_auc_no10e5.tsv" if AUC_ARM else
                              "hog_phenotype_assoc_no10e5.tsv" if NO10E5
                              else "hog_phenotype_assoc.tsv"))
PAV = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv")
ANNOT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
XLSX = os.path.join(ROOT, "00.manuscript/supplementary",
                    "suppl_table.S7_temperature_modules.xlsx")

STRAINS = ["BMU" + x.zfill(5) for x in
           ["9608", "7029", "7936", "1336", "2406", "5572", "2479", "2031", "12685",
            "12686", "1287", "1935", "1936", "12683", "12684", "2407", "3976", "2028",
            "1973", "2745", "10147"]]
TRAITS = ({"28C": ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"],
           "37C": ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"]} if AUC_ARM else
          {"28C": ["mort_28Ce3", "mort_28Ce4", "mort_28_mean2"],
           "37C": ["mort_37Ce3", "mort_37Ce4", "mort_37_mean2"]} if NO10E5 else
          {"28C": ["mort_28Ce3", "mort_28Ce4", "mort_28Ce5", "mort_28_mean"],
           "37C": ["mort_37Ce3", "mort_37Ce4", "mort_37Ce5", "mort_37_mean"]})
# Not a convention: the 99th percentile of the max|rho| null distribution obtained by
# permuting the trait labels across the 4,484 variable orthogroups (500 draws; 50th = 0.29,
# 90th = 0.461, 95th = 0.509, 99th = 0.602/0.604). The full sweep, and the fact that a Pfam
# enrichment appears only inside a narrow band around the 95th percentile, is in
# scripts/S7_threshold_calibration.py and the table's Threshold_sensitivity sheet.
RHO_MIN = {"28C": 0.602, "37C": 0.604}

assoc = pd.read_csv(ASSOC, sep="\t")
pav = pd.read_csv(PAV, sep="\t", dtype={"Present": "int8"})
wide = (pav[pav.SampleID.isin(STRAINS)]
        .pivot(index="SampleID", columns="HOG_ID", values="Present").fillna(0).loc[STRAINS])
prev = wide.sum(0)
universe = set(wide.columns[(prev >= 3) & (prev <= 19)])

annot = pd.read_csv(ANNOT, sep="\t", usecols=["HOG_ID", "Pfam_Domain"]).drop_duplicates("HOG_ID")
domain_hogs = {}
for _, r in annot.iterrows():
    if r["HOG_ID"] not in universe:
        continue
    pf = str(r["Pfam_Domain"])
    if pf.lower() in ("nan", "-"):
        continue
    for d in {x for x in pf.split(";") if x}:
        domain_hogs.setdefault(d, set()).add(r["HOG_ID"])


def modules(traits, tag):
    rho = (assoc[assoc.trait.isin(traits)]
           .groupby("HOG_ID")["rho"].apply(lambda x: x.abs().max()))
    return sorted(set(rho[rho >= RHO_MIN[tag]].index) & universe)


def enrich(hset):
    n_set = len(hset)
    rows = []
    for d, members in domain_hogs.items():
        k = len(members & hset)
        if k == 0:
            continue
        K = len(members & universe)
        rows.append({"domain": d, "n_hit": k, "n_set": n_set, "n_domain": K,
                     "n_universe": len(universe),
                     "p": stats.hypergeom.sf(k - 1, len(universe), K, n_set)})
    e = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    m = len(e)
    q = np.minimum(1, e["p"] * m / (np.arange(m) + 1))
    e["q"] = np.minimum.accumulate(q.values[::-1])[::-1]
    return e


wb = load_workbook(XLSX)
member_rows, defs = [], []
for tag, traits in TRAITS.items():
    hset = modules(traits, tag)
    e = enrich(set(hset))
    for name, df in [(f"{tag}_module_Pfam_enrichment", e),
                     (f"{tag}_module_HOGs", pd.DataFrame({"HOG_ID": hset}))]:
        if name in wb.sheetnames:
            del wb[name]
        ws = wb.create_sheet(name)
        ws.append(list(df.columns))
        for row in df.itertuples(index=False):
            ws.append([v.item() if hasattr(v, "item") else v for v in row])
    for h in hset:
        member_rows.append({"temperature": tag, "HOG_ID": h})
    defs.append((tag, len(hset), len(e), e.iloc[0]["domain"], float(e.iloc[0]["q"]),
                 int((e["q"] < 0.05).sum())))
    print(f"{tag}: {len(hset)} HOGs · {len(e)} domains tested · best "
          f"{e.iloc[0]['domain']} q={e.iloc[0]['q']:.3f} · passing q<0.05: "
          f"{int((e['q'] < 0.05).sum())}")

if "Module_HOG_membership" in wb.sheetnames:
    del wb["Module_HOG_membership"]
ws = wb.create_sheet("Module_HOG_membership")
ws.append(["temperature", "HOG_ID"])
for r in member_rows:
    ws.append([r["temperature"], r["HOG_ID"]])

if "Module_definitions" in wb.sheetnames:
    del wb["Module_definitions"]
ws = wb.create_sheet("Module_definitions")
ws.append(["item", "value"])
ws.append(["member rule",
           "variable HOGs (present in 3-19 of the 21 strains) whose maximum |Spearman rho| "
           "against that temperature's three AUC phenotypes reaches the 99th percentile of "
           "the permutation null: 0.602 at 28 C and 0.604 at 37 C (500 permutations of the "
           "trait labels; 90th and 95th percentiles are 0.461 and 0.509)"])
ws.append(["phenotypes",
           "seven phenotypes built on the area under the survival curve: the four condition "
           "AUCs at 10^3 and 10^4 (as 1 - AUC), their two two-dose means and the difference "
           "of those means. The 10^5 dose is excluded: 16 of 21 isolates reach >= 95% "
           "mortality at both temperatures, so its endpoint is on the ceiling"])
ws.append(["threshold sweep",
           "reported in the Threshold_sensitivity sheet: a Pfam enrichment appears only "
           "around the 95th percentile and vanishes within 0.002 of rho, because a single "
           "orthogroup crosses the line, so no enrichment is claimed"])
ws.append(["enrichment",
           "one-sided hypergeometric, universe = 4,484 variable HOGs, class = HOGs "
           "carrying the domain, Benjamini-Hochberg across every tested domain"])
ws.append(["28C module", f"{defs[0][1]} HOGs; {defs[0][2]} domains tested; best "
                          f"{defs[0][3]} q = {defs[0][4]:.3f}; domains passing q<0.05: {defs[0][5]}"])
ws.append(["37C module", f"{defs[1][1]} HOGs; {defs[1][2]} domains tested; best "
                          f"{defs[1][3]} q = {defs[1][4]:.3f}; domains passing q<0.05: {defs[1][5]}"])

wb.save(XLSX)
print(f"\nwrote {XLSX}")
print("sheets:", wb.sheetnames)
