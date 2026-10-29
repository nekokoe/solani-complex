#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core-chromosome control for the repeat-flanking enrichment of accessory-derived genes.

Why this control
    Comparing transferred genes with the whole genome is confounded: transferred genes
    sit on accessory chromosomes, which are repeat-rich. Accessory-resident
    "non-transferred" genes are a poor control too, because with only 13 assembled
    genomes a gene called non-transferred may simply have been transferred in an
    un-sampled genome. The biologically meaningful control is therefore ACCESSORY-
    DERIVED GENES THAT NOW RESIDE ON CORE CHROMOSOMES (the core-side copies of C-A
    clusters) against other core-resident genes that are never accessory-associated.
    Both sides sit in the same repeat-poor core environment, so the comparison is
    unaffected by the 13-genome sampling limit.

Input
    suppl_table.gene_transfer_211.tsv from scripts/stats_gene_transfer_211.py.
      test    Contig_type == core AND Transfer_type == C-A
      control every OTHER core gene (single, no-acc and the unclustered ones):
              test + control then equals the whole core-gene set, which is how the
              archived table was built (1,981 + 195,618 = 197,599)

    The producing script for this table was missing from the repository (the archived
    table had no code), so this is a re-implementation; on the archived 100%-identity
    table it reproduces the archived numbers exactly (n_test 1,981, n_control 195,618,
    DNA within 1 kb 15.45% vs 2.55%, OR 6.986).

Usage
    python3 scripts/transfer_core_control_211.py \
        --per-gene <.../suppl_table.gene_transfer_211.tsv> \
        --out      <.../suppl_table.repeat_enrichment_211.core_control.tsv>
"""

import argparse
import os
from collections import Counter

import pandas as pd
from scipy import stats

FAMILIES = ["DNA", "LTR", "LINE", "RC", "Unknown", "Simple_repeat", "Low_complexity"]
TEST_TYPE = "C-A"


def frac_flanked(series, fam):
    return sum(1 for s in series if fam in str(s).split(";"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-gene", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.per_gene, sep="\t")
    core = df[df["Contig_type"] == "core"]
    test = core[core["Transfer_type"] == TEST_TYPE]
    ctrl = core[core["Transfer_type"] != TEST_TYPE]
    assert len(test) + len(ctrl) == len(core), "test and control must partition the core genes"
    print(f"test (C-A on core): {len(test)}   control (other core genes): {len(ctrl)}"
          f"   = {len(core)} core genes")

    rows = []
    for fam in FAMILIES:
        for win, col in (("flank_1kb", "Flank1kb_families"), ("flank_5kb", "Flank5kb_families")):
            k1 = frac_flanked(test[col], fam)
            k2 = frac_flanked(ctrl[col], fam)
            n1, n2 = len(test), len(ctrl)
            tbl = [[k1, n1 - k1], [k2, n2 - k2]]
            orr, p = stats.fisher_exact(tbl)
            rows.append({
                "Family": fam, "Metric": win,
                "n_CA_core": n1, "n_core_background": n2,
                "CA_core_frac": k1 / n1, "Core_back_frac": k2 / n2,
                "Odds_ratio": orr, "Fisher_p": p,
            })
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
