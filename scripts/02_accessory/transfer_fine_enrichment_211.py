#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fine repeat-family enrichment of transferred genes (the S13 panel c input).

For every RepeatMasker fine family, compare the fraction of TRANSFERRED genes whose
1 kb / 5 kb flank contains that family against the fraction among the whole-genome
background (all other genes of the 13 chromosome-level assemblies). Two-sided Fisher
exact test per family; families are then ranked by odds ratio at each window.

Provenance
    suppl_table.repeat_enrichment_211.fine_both_flanks.tsv had no producing script in
    the repository (the grouped "both_flanks" table is just a relabelled copy of
    suppl_table.repeat_enrichment_211.tsv, but the FINE one is a distinct product).
    This is a re-implementation; on the archived 100%-identity table it reproduces the
    archived file (DNA/hAT-Tag1 OR 21.2505 / 22.4574, n 109 / 329).

Usage
    python3 scripts/transfer_fine_enrichment_211.py \
        --per-gene <.../suppl_table.gene_transfer_211.tsv> \
        --out      <.../suppl_table.repeat_enrichment_211.fine_both_flanks.tsv>
"""

import argparse
import os
from collections import Counter

import pandas as pd
from scipy import stats


def explode(series):
    c = Counter()
    for s in series:
        s = str(s)
        if s and s != "nan":
            for f in s.split(";"):
                if f:
                    c[f] += 1
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-gene", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    df = pd.read_csv(a.per_gene, sep="\t")
    tr = df[df["AC_related_transfer"].astype(str).str.lower() == "true"]
    bg = df[df["AC_related_transfer"].astype(str).str.lower() != "true"]
    print(f"transfer genes: {len(tr)}   background: {len(bg)}")

    rows = []
    for win, col in (("1kb", "Flank1kb_fine"), ("5kb", "Flank5kb_fine")):
        rows.append((win, explode(tr[col]), explode(bg[col])))
    fams = sorted({f for _, t, b in rows for f in set(t) | set(b)})

    out = []
    for f in fams:
        rec = {"Fine_family": f}
        for win, tc, bc in rows:
            k1, k2 = tc.get(f, 0), bc.get(f, 0)
            n1, n2 = len(tr), len(bg)
            orr, p = stats.fisher_exact([[k1, n1 - k1], [k2, n2 - k2]])
            rec[f"n_transfer_{win}"] = k1
            rec[f"Transfer_{win}"] = k1 / n1
            rec[f"Genome_{win}"] = k2 / n2
            rec[f"OR_{win}"] = orr
            rec[f"p_{win}"] = p
        out.append(rec)
    o = pd.DataFrame(out)
    o["rank1"] = o["OR_1kb"].rank(ascending=False, method="min")
    o["rank5"] = o["OR_5kb"].rank(ascending=False, method="min")
    o["mean_rank"] = o[["rank1", "rank5"]].mean(axis=1)
    o = o.sort_values("OR_1kb", ascending=False)

    cols = ["Fine_family", "n_transfer_1kb", "Transfer_1kb", "Genome_1kb", "OR_1kb", "p_1kb",
            "n_transfer_5kb", "Transfer_5kb", "Genome_5kb", "OR_5kb", "p_5kb",
            "rank1", "rank5", "mean_rank"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    o[cols].to_csv(a.out, sep="\t", index=False)
    print(f"wrote {a.out}  ({len(o)} fine families)")
    print(o[["Fine_family", "n_transfer_1kb", "OR_1kb", "OR_5kb"]].head(6).to_string(index=False))


if __name__ == "__main__":
    main()
