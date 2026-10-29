#!/usr/bin/env python3
"""Supplementary Table S6 — candidate accessory HOGs of the PAV screen.

Source of truth: the PAV rerun of 2026-10-21 (`09.rf_clinical_genes/rf_pav_rerun.py`,
outputs in `09.rf_clinical_genes/pav_rerun/`). Candidate HOGs are selected on
presence/absence (PAV), not on gene copy number — the count-based run gave 319 / 196,
the PAV run gives 390 / 277, and the manuscript cites the PAV numbers.

Anchors reproduced here (they are the numbers the Results text quotes):
    candidates                 390 (F. solani) · 277 (F. falciformis) · 12 shared · 655 distinct
    higher presence in clinical 272 of 390 · 176 of 277
    Bonferroni-significant      48 of 390 · 110 of 277
                                (p < 0.05 / n_variable_HOGs = 3.89e-6 and 3.97e-6;
                                 correcting over the candidate set instead would give
                                 112 and 155, so the denominator matters)
    composition over the union  44.6% carry an effector, 19.1% a protein kinase

Writes 00.manuscript/supplementary/suppl_table.S6_candidate_HOGs_PAV.xlsx
"""

import os

import numpy as np
import openpyxl
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAV = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
OUT = os.path.join(ROOT, "00.manuscript/supplementary",
                   "suppl_table.S6_candidate_HOGs_PAV.xlsx")
N_PAN = 29080          # pan-genome HOG count, used for the % of pan-genome


def main():
    stats = pd.read_csv(f"{PAV}/rf_pav_hog_stats.tsv", sep="\t")
    summary = pd.read_csv(f"{PAV}/rf_pav_summary.tsv", sep="\t").set_index("species")

    frames, anchors = {}, []
    for sp, label in (("solani", "F. solani"), ("falciformis", "F. falciformis")):
        feat = pd.read_csv(f"{PAV}/rf_pav_{sp}_features.tsv", sep="\t")
        ann = stats[stats.Species.isin([sp, "both"])][
            ["HOG_ID", "IsEffector", "EffectorClass", "IsKinase", "PfamClass",
             "Pfam_Domain", "KEGG_Category", "Bonferroni_p"]]
        d = feat.merge(ann, on="HOG_ID", how="left").sort_values("importance", ascending=False)
        d.insert(0, "species", label)
        frames[sp] = d

        n_var = int(summary.at[sp, "n_variable_hogs"])
        thr = 0.05 / n_var
        anchors.append(dict(
            species=label, candidates=len(d),
            pct_of_pangenome=round(len(d) / N_PAN * 100, 2),
            n_variable_hogs=n_var, bonferroni_threshold=f"{thr:.2e}",
            bonferroni_significant=int((d.p_pav < thr).sum()),
            higher_presence_in_clinical=int((d.rate_delta > 0).sum()),
            auc_cv=summary.at[sp, "auc_cv"], auc_cv_sd=summary.at[sp, "auc_cv_sd"],
            clinical_rate_mean=round(d.clinical_rate.mean(), 3),
            nonclinical_rate_mean=round(d.nonclinical_rate.mean(), 3),
            n_effector=int(d.IsEffector.sum()), n_kinase=int(d.IsKinase.sum())))
        print(f"{label}: {len(d)} candidates, Bonferroni p<{thr:.2e} -> "
              f"{anchors[-1]['bonferroni_significant']}, higher in clinical "
              f"{anchors[-1]['higher_presence_in_clinical']}")

    union = stats.copy()
    print(f"\nunion: {len(union)} distinct HOGs; effector {union.IsEffector.mean():.1%}, "
          f"kinase {union.IsKinase.mean():.1%}")
    shared = union[union.Species == "both"].copy()

    definitions = pd.DataFrame([
        ("candidate rule",
         "HOGs whose presence/absence separates clinical from non-clinical isolates as "
         "called by the random-forest classifier, i.e. non-zero feature importance with a "
         "significant PAV difference; selection is on presence/absence, not copy number"),
        ("cohorts",
         "F. solani 86 isolates (60 clinical, 26 non-clinical); F. falciformis 66 "
         "(48 clinical, 18 non-clinical); F. keratoplastica not modelled (4 non-clinical)"),
        ("variable HOGs in the test",
         f"12,845 (F. solani) and 12,579 (F. falciformis) of {N_PAN:,} orthogroups"),
        ("Bonferroni denominator",
         "the number of variable HOGs tested (12,845 / 12,579), giving thresholds "
         "3.89e-6 and 3.97e-6; correcting over the candidate set instead would yield 112 "
         "and 155 significant HOGs rather than the 48 and 110 reported"),
        ("composition percentages",
         "44.6% effector / 19.1% protein kinase are computed over the 655 distinct HOGs "
         "of the union (12 HOGs are candidates in both species)"),
        ("columns",
         "clinical_rate and nonclinical_rate are the fraction of isolates of that class "
         "carrying the HOG; rate_delta = clinical_rate - nonclinical_rate; p_pav is the "
         "presence/absence test, q_pav its Benjamini-Hochberg value, Bonferroni_p the "
         "multiplied p; presence_interval is the pan-genome presence bin"),
    ], columns=["item", "value"])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with pd.ExcelWriter(OUT) as w:
        pd.DataFrame(anchors).to_excel(w, sheet_name="anchors", index=False)
        definitions.to_excel(w, sheet_name="definitions", index=False)
        frames["solani"].to_excel(w, sheet_name="candidates_solani", index=False)
        frames["falciformis"].to_excel(w, sheet_name="candidates_falciformis", index=False)
        shared[["HOG_ID", "Importance", "P_value", "Bonferroni_p", "Clinical_rate",
                "NonClinical_rate", "Rate_delta", "IsEffector", "IsKinase", "Pfam_Domain",
                "KEGG_Category"]].to_excel(w, sheet_name="shared_both_species", index=False)
    print(f"\nwrote {os.path.relpath(OUT, ROOT)}")
    wb = openpyxl.load_workbook(OUT)
    for ws in wb:
        print(f"  {ws.title:24s} {ws.max_row-1:4d} rows")


if __name__ == "__main__":
    main()
