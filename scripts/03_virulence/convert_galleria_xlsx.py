#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""convert_galleria_xlsx.py — 21-strain Galleria raw xlsx -> the long table the
pathogenicity pipeline expects.

The 2026-10-12 raw file (07.lethal_groups/21株-大蜡螟原始数据.xlsx) is a block
matrix: six stacked blocks of 2 temperatures x 3 doses, each block holding

    row A : temperature label + 21 strain ids written as bare last-5-digits
    row B : "Day"         + the dose label (10^3 / 10^4 / 10^5)
    rows  : day 0..10     + surviving-larva counts
    row   : "Probability of Survival" (absent in the two 10^5 blocks)

The pipeline instead reads the long format `Sample / Measure / Count / Percent`
with Measure like `37Ce3_day00` (scripts/phenotype_modeling.py:27 parses it with
(\\d+C[eE]\\d+)_day(\\d+)). This script does that conversion.

Two rules that matter and are easy to get wrong:

  * Percent is denormalised against EACH column's own day-0 count. The initial
    number of larvae genuinely varies (9-14) between strains and blocks, so a
    fixed denominator of 10 or 12 is wrong.
  * Missing values (a column simply stops once every larva is dead) become
    Count 0 / Percent 0. Leaving them as NaN corrupts both the day-10 mortality
    (phenotype_modeling.py:47 uses 1 - survival_at_day10) and the AUC
    (phenotype_modeling.py:37-41 drops NaN days instead of scoring them 0).
    Every NaN run in the raw file is preceded by an explicit 0, which this
    script asserts.

Usage:
  python3 scripts/convert_galleria_xlsx.py                    # all 21 strains
  python3 scripts/convert_galleria_xlsx.py --strains 07.lethal_groups/selected_sample.list \\
      --out /tmp/old15.txt                                    # 15-strain validation run
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_XLSX = os.path.join(ROOT, "07.lethal_groups", "21株-大蜡螟原始数据.xlsx")
# The converter writes the CANONICAL path in place rather than a new name: only
# two scripts read it (phenotype_modeling.py:21, plot_phenotype_s25.py:40) and both
# are name-coupled, so replacing in place keeps the whole chain working with zero
# reader edits. The 15-strain original is preserved at
# 07.lethal_groups/.bak-21strain/caterpillar_37C_narrow.txt.
DEFAULT_OUT = os.path.join(ROOT, "07.lethal_groups", "caterpillar_37C_narrow.txt")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", default=DEFAULT_XLSX)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--strains", default=None,
                   help="optional strain-list file; restrict the output to these "
                        "ids (used to re-run the converter on the old 15 strains)")
    return p.parse_args()


def norm_strain(v):
    """9608 -> BMU09608; pass through anything already carrying the prefix."""
    s = str(v).strip()
    if s.upper().startswith("BMU"):
        return s.upper()
    if s.endswith(".0"):
        s = s[:-2]
    return "BMU" + s.zfill(5)


def read_blocks(raw):
    """-> list of (cond, {strain: {day: count}}), in file order."""
    n_row, n_col = raw.shape
    starts = [i for i in range(n_row) if "℃" in str(raw.iat[i, 0])]
    if not starts:
        sys.exit("no temperature-labelled blocks found - is this the right file?")

    blocks = []
    for bi, s in enumerate(starts):
        tlabel = str(raw.iat[s, 0])
        temp = "37C" if "37" in tlabel else "28C" if "28" in tlabel else None
        if temp is None:
            sys.exit(f"cannot read a temperature from {tlabel!r} at row {s}")
        dose_lab = str(raw.iat[s + 1, 1]).strip()
        if not dose_lab.startswith("10^"):
            sys.exit(f"cannot read a dose from {dose_lab!r} at row {s + 1}")
        cond = f"{temp}e{dose_lab.split('^')[1]}"
        ids = {c: norm_strain(raw.iat[s, c]) for c in range(1, n_col)}

        end = starts[bi + 1] if bi + 1 < len(starts) else n_row
        series = {sid: {} for sid in ids.values()}
        for r in range(s + 2, end):
            try:
                day = int(float(raw.iat[r, 0]))
            except (TypeError, ValueError):
                continue                       # the 'Probability of Survival' row
            for c, sid in ids.items():
                series[sid][day] = raw.iat[r, c]
        blocks.append((cond, series))
    return blocks


def build(blocks, keep=None):
    rows = []
    for cond, series in blocks:
        for sid, days in series.items():
            if keep is not None and sid not in keep:
                continue
            if not days:
                continue
            counts = []
            for day in sorted(days):
                v = days[day]
                counts.append((day, v))
            n0 = counts[0][1]
            if pd.isna(n0) or float(n0) <= 0:
                sys.exit(f"day-0 count missing/zero for {sid} {cond}")
            saw_zero = False
            for i, (day, v) in enumerate(counts):
                if pd.isna(v):
                    if not saw_zero and i > 0:
                        sys.exit(f"NaN before any zero: {sid} {cond} day {day}")
                    v = 0.0
                if float(v) == 0:
                    saw_zero = True
                rows.append({"Sample": sid, "Measure": f"{cond}_day{day:02d}",
                             "Count": int(float(v)),
                             "Percent": f"{float(v) / float(n0):.9f}"})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    raw = pd.read_excel(args.xlsx, header=None)
    blocks = read_blocks(raw)
    print("blocks in file order:", [c for c, _ in blocks])

    keep = None
    if args.strains:
        keep = {norm_strain(x) for x in open(args.strains) if x.strip()}
        print(f"restricting to {len(keep)} strains from {args.strains}")

    out = build(blocks, keep)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    n0s = out[out.Measure.str.endswith("_day00")].Count
    print(f"strains: {out.Sample.nunique()}  conditions: "
          f"{out.Measure.str.rsplit('_day', n=1).str[0].nunique()}  rows: {len(out)}")
    print(f"initial N: min {n0s.min()} max {n0s.max()}  (distinct: {sorted(set(n0s))})")
    print(f"wrote {os.path.relpath(args.out, ROOT)}")


if __name__ == "__main__":
    main()
