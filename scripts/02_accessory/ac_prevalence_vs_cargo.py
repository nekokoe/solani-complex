#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C-3: is an accessory chromosome's PREVALENCE across the 211 FSSC isolates
related to its CARGO FUNCTION?

Motivation
----------
Puebla Planas et al. (bioRxiv 2026.10.04.756499) argue dispensable chromosomes
are maintained because they confer a fitness benefit, i.e. maintenance is the
price of their cargo. Our manuscript already describes the AC cargo as
metabolism / transport / xenobiotic-degradation genes and the ACs as clinically
enriched. This script turns "why are they kept?" into a testable statement on
our own data: does the breadth of sharing of an AC (how many of the 211
isolates carry it) track how metabolic its cargo is?

Inputs (all in-repo)
    prevalence : 06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv
                 78 ACs x 211 isolates, k-mer sharing proportion.
                 Carrier convention follows plot_kmer_analysis_figures.py:
                 an isolate carries the AC when sharing >= --share-threshold (0.4).
    genes      : 02.annotation/per_isolates/<ISO>/funannotate/predict_results/fusarium.gff3
    functions  : 02.annotation/per_isolates/<ISO>/<ISO>.emapper.emapper.annotations
                 (eggNOG-mapper, COG_category column)

Cargo classes (COG functional categories; the transport+metabolism set)
    metab_transport = C,E,F,G,H,I,P,Q
    secondary_Q     = Q  (secondary metabolites biosynthesis, transport and catabolism)
    unknown         = S plus unannotated ("-")
    other           = everything else
Rationale: C/E/F/G/H/I/P/Q are exactly the COG categories whose names contain
"transport and metabolism" (plus C = energy production and Q = secondary
metabolites), i.e. the functional space our manuscript attributes to AC cargo.

GATE: the eggNOG query IDs must join to the GFF3 mRNA IDs at a high rate. A
silent 0% join is the classic failure mode of this pipeline and would make every
downstream fraction meaningless, so the join rate is checked before reporting.

Outputs (06.accessory-chrs/starship_ac/, prefix ac_cargo_)
    ac_cargo_per_ac.tsv        per AC: carrier count, prevalence, cargo fractions
    ac_cargo_correlation.tsv   Spearman tests of prevalence vs each cargo fraction
    ac_cargo_summary.md        human-readable summary
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

METAB = list("CEFGHIPQ")
CARGO_CLASSES = ["metab_transport", "secondary_Q", "unknown", "other"]
JOIN_RATE_FLOOR = 0.50   # below this, the annotation join is broken -> abort


def parse_args():
    p = argparse.ArgumentParser(description="AC prevalence vs cargo function")
    p.add_argument("--annotation-root",
                   default=os.path.join(REPO, "02.annotation", "per_isolates"))
    p.add_argument("--matrix",
                   default=os.path.join(REPO, "06.accessory-chrs", "kmer_analysis",
                                        "kmer_matrix.ACx211.tsv"))
    p.add_argument("--out-dir",
                   default=os.path.join(REPO, "06.accessory-chrs", "starship_ac"))
    p.add_argument("--share-threshold", type=float, default=0.4)
    p.add_argument("--skip-gate", action="store_true")
    return p.parse_args()


def read_gene_contig(gff):
    d = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "mRNA":
                continue
            for a in p[8].split(";"):
                if a.startswith("ID="):
                    d[a[3:]] = p[0]
                    break
    return d


def read_cog(emapper):
    """query -> COG_category string (may be multi-letter / '-' / absent)."""
    d = {}
    with open(emapper) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            d[p[0]] = p[6]
    return d


def cargo_class(cog):
    if cog in ("-", ""):
        return "unknown"
    if "Q" in cog:
        return "secondary_Q"
    if any(c in cog for c in METAB):
        return "metab_transport"
    if "S" in cog:
        return "unknown"
    return "other"


def main():
    args = parse_args()
    if not os.path.exists(args.matrix):
        raise SystemExit(f"missing matrix: {args.matrix}")

    m = pd.read_csv(args.matrix, sep="\t", index_col=0)
    m.index.name = "AC"
    isolates = list(m.columns)
    print(f"matrix: {m.shape[0]} ACs x {m.shape[1]} isolates")

    carriers = (m >= args.share_threshold).sum(axis=1)
    prevalence = carriers / m.shape[1]

    # ------------------------------------------------------------ gather genes
    per_iso = defaultdict(dict)      # iso -> {ac_contig: [gene ids]}
    joined = total_genes = 0
    for ac_id in m.index:
        iso = ac_id.split("_")[0]
        per_iso[iso].setdefault(ac_id, [])
    for iso in sorted(per_iso):
        ann = os.path.join(args.annotation_root, iso)
        gff = os.path.join(ann, "funannotate", "predict_results", "fusarium.gff3")
        emp = os.path.join(ann, f"{iso}.emapper.emapper.annotations")
        for p in (gff, emp):
            if not os.path.exists(p):
                raise SystemExit(f"missing input: {p}")
        gc = read_gene_contig(gff)
        cog = read_cog(emp)
        total_genes += len(gc)
        joined += sum(1 for t in gc if t in cog)
        for tid, ct in gc.items():
            if ct in per_iso[iso]:
                per_iso[iso][ct].append(tid)
        per_iso[iso]["__cog__"] = cog

    rate = joined / total_genes if total_genes else 0.0
    print(f"eggNOG join rate: {joined}/{total_genes} = {rate*100:.1f}%")
    if not args.skip_gate and rate < JOIN_RATE_FLOOR:
        print(f"\nGATE FAILED: join rate {rate*100:.1f}% < {JOIN_RATE_FLOOR*100:.0f}%. "
              "The GFF3 transcript IDs and the eggNOG query IDs are not the same "
              "namespace here; every cargo fraction below would be meaningless.")
        return 1
    if not args.skip_gate:
        print("GATE PASSED -- annotation join is sound\n")

    # ------------------------------------------------------ per-AC cargo table
    rows = []
    for ac_id in m.index:
        iso = ac_id.split("_")[0]
        genes = per_iso[iso].get(ac_id, [])
        cog = per_iso[iso]["__cog__"]
        counts = defaultdict(int)
        for t in genes:
            counts[cargo_class(cog.get(t, "-"))] += 1
        n = len(genes)
        rec = {"AC": ac_id, "isolate": iso, "n_genes": n,
               "carriers_of_211": int(carriers[ac_id]),
               "prevalence": round(float(prevalence[ac_id]), 4)}
        for c in ["metab_transport", "secondary_Q", "unknown", "other"]:
            rec[f"frac_{c}"] = round(counts[c] / n, 4) if n else np.nan
        rec["frac_metab_plus_secondary"] = round(
            (counts["metab_transport"] + counts["secondary_Q"]) / n, 4) if n else np.nan
        rows.append(rec)
    df = pd.DataFrame(rows)
    os.makedirs(args.out_dir, exist_ok=True)
    df.to_csv(os.path.join(args.out_dir, "ac_cargo_per_ac.tsv"), sep="\t", index=False)

    # --------------------------------------------------------- correlations
    corr_rows = []
    for col in ["frac_metab_transport", "frac_secondary_Q",
                "frac_metab_plus_secondary", "frac_unknown", "n_genes"]:
        sub = df[["prevalence", col]].dropna()
        if len(sub) < 5 or sub[col].nunique() < 2:
            corr_rows.append({"variable": col, "n": len(sub), "rho": np.nan, "p": np.nan})
            continue
        rho, pv = stats.spearmanr(sub["prevalence"], sub[col])
        corr_rows.append({"variable": col, "n": len(sub),
                          "rho": round(float(rho), 4), "p": round(float(pv), 4)})
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(os.path.join(args.out_dir, "ac_cargo_correlation.tsv"), sep="\t", index=False)

    # ------------------------------------------------------------- reporting
    print(f"\nper-AC table: {len(df)} ACs, median genes/AC = {int(df.n_genes.median())}")
    print(f"carriers at k-mer sharing >= {args.share_threshold}: "
          f"median {int(df.carriers_of_211.median())} of 211")
    print("\n" + corr.to_string(index=False))

    lines = ["# Accessory-chromosome prevalence vs cargo function", "",
             f"Prevalence = share of the 211 FSSC isolates carrying the AC at "
             f"k-mer sharing >= {args.share_threshold} "
             f"(the carrier convention of plot_kmer_analysis_figures.py).", "",
             f"Cargo from eggNOG COG_category per gene "
             f"(join rate {rate*100:.1f}%). "
             f"`metab_transport` = C,E,F,G,H,I,P,Q; `secondary_Q` = Q.", "",
             "| — | value |", "|---|---|",
             f"| ACs analysed | {len(df)} |",
             f"| median genes per AC | {int(df.n_genes.median())} |",
             f"| median carriers (of 211) | {int(df.carriers_of_211.median())} |", "",
             "## Spearman correlation with prevalence", "",
             "| variable | n | rho | p |", "|---|---|---|---|"]
    for _, r in corr.iterrows():
        lines.append(f"| {r.variable} | {r.n} | {r.rho} | {r.p} |")
    lines += ["", "## Reading", "",
              "A positive rho for `frac_metab_transport` would support the "
              "fitness-benefit / cargo-value view of AC maintenance. A null or "
              "negative result is a genuine negative and should be reported as "
              "such rather than explained away.", "",
              "Caveats: (i) prevalence at a single threshold -- the threshold "
              "sensitivity should be reported alongside; (ii) COG categories are "
              "coarse and 'unknown' (S) is a large residual class here; "
              "(iii) 78 ACs from 13 assemblies means the 13 isolates are not "
              "independent replicates, so the p-values are anti-conservative.", ""]
    with open(os.path.join(args.out_dir, "ac_cargo_summary.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    # threshold sensitivity for prevalence
    sens = []
    for th in [0.2, 0.3, 0.4, 0.5, 0.6]:
        car = (m >= th).sum(axis=1) / m.shape[1]
        sub = pd.DataFrame({"prev": car}).join(df.set_index("AC").frac_metab_plus_secondary)
        sub = sub.dropna()
        rho, pv = stats.spearmanr(sub.prev, sub.frac_metab_plus_secondary)
        sens.append({"share_threshold": th, "median_carriers": int((m >= th).sum(axis=1).median()),
                     "rho_vs_metab_plus_secondary": round(float(rho), 4),
                     "p": round(float(pv), 4)})
    pd.DataFrame(sens).to_csv(os.path.join(args.out_dir, "ac_cargo_threshold_sensitivity.tsv"),
                              sep="\t", index=False)
    print("\nthreshold sensitivity:")
    print(pd.DataFrame(sens).to_string(index=False))
    print(f"\nwrote ac_cargo_per_ac.tsv, ac_cargo_correlation.tsv, "
          f"ac_cargo_threshold_sensitivity.tsv, ac_cargo_summary.md in {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
