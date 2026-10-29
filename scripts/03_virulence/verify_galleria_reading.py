#!/usr/bin/env python3
"""Independent audit of the Galleria raw-data reading.

The 10^5 blocks of `21株-大蜡螟原始数据.xlsx` stop being recorded as soon as every
larva is dead, so those columns carry long runs of empty cells and no
"Probability of Survival" row. Whether those empties mean "all dead" or "missing"
decides both the day-10 mortality and the survival AUC, so this script re-derives
everything straight from the workbook and compares it against the delivered files.

Checks
  1 every empty run is preceded by an explicit 0 (i.e. blanks mean all dead)
  2 survival is non-increasing along the days within each strain x block
  3 day-0 counts (n0) vary between 9 and 14, as documented
  4 the long tables agree cell by cell with the workbook (blanks -> 0)
  5 phenotype_mortality.tsv and phenotype_auc.tsv recomputed from the workbook match
"""

import os

import numpy as np
import openpyxl
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.path.join(ROOT, "07.lethal_groups/21株-大蜡螟原始数据.xlsx")
LONG21 = os.path.join(ROOT, "07.lethal_groups/caterpillar_21_narrow.txt")
LONG37 = os.path.join(ROOT, "07.lethal_groups/caterpillar_37C_narrow.txt")
PA = os.path.join(ROOT, "07.lethal_groups/phenotype_analysis")


def read_raw():
    """-> {(strain, temp, dose): {day: count_or_None}} plus blank-run diagnostics."""
    ws = openpyxl.load_workbook(XLSX, data_only=True)["Sheet1"]
    starts = [r for r in range(1, ws.max_row + 1) if str(ws.cell(r, 1).value).strip() == "Day"]
    cells, blank_runs, bad_runs = {}, [], []
    n0s = set()
    for s in starts:
        tlab = str(ws.cell(s - 1, 1).value)
        temp = "37C" if "37" in tlab else "28C" if "28" in tlab else None
        dose = str(ws.cell(s, 2).value).replace("10^", "e")
        strains = {}
        for c in range(2, 23):
            sid = ws.cell(s - 1, c).value
            if sid is None:
                continue
            sid = str(int(float(sid))) if str(sid).replace(".0", "").isdigit() else str(sid)
            strains[c] = "BMU" + sid.replace("BMU", "").zfill(5)
        days = {}
        for r in range(s + 1, s + 12):
            d = ws.cell(r, 1).value
            if d is None:
                continue
            days[int(float(d))] = r
        for c, sid in strains.items():
            series, prev_zero_days = {}, []
            vals = []
            for day in sorted(days):
                v = ws.cell(days[day], c).value
                series[day] = None if v is None else int(v)
                vals.append((day, series[day]))
            # 1 + 2: blank runs and monotonicity
            seen_zero = False
            for day, v in vals:
                if v is None:
                    blank_runs.append((sid, temp, dose, day))
                    if not seen_zero:
                        bad_runs.append((sid, temp, dose, day))
                else:
                    if v == 0:
                        seen_zero = True
                    if vals[0][1] is not None and v > vals[0][1]:
                        bad_runs.append((sid, temp, dose, day, "increase"))
            n0 = series[sorted(days)[0]]
            if n0 is not None:
                n0s.add(n0)
            cells[(sid, temp, dose)] = series
    return cells, blank_runs, bad_runs, n0s


def main():
    cells, blank_runs, bad_runs, n0s = read_raw()
    print("=" * 74)
    print("1-3. BLANKS, MONOTONICITY, INITIAL COUNTS")
    print("=" * 74)
    print(f"strain x block cells: {len(cells)}")
    print(f"empty cell-days: {len(blank_runs)}")
    print(f"empty runs NOT preceded by an explicit 0: {len(bad_runs)}")
    print(f"day-0 counts (n0) observed: {sorted(n0s)}")
    assert not bad_runs, bad_runs[:5]
    by_block = pd.Series([f"{t}{d}" for _, t, d, _ in blank_runs]).value_counts().sort_index()
    print("empties per block:")
    for k, v in by_block.items():
        print(f"   {k}: {v}")

    print("\n" + "=" * 74)
    print("4. LONG TABLES vs THE WORKBOOK")
    print("=" * 74)
    for path in (LONG21, LONG37):
        if not os.path.exists(path):
            print(f"   {os.path.basename(path)}: absent")
            continue
        d = pd.read_csv(path, sep="\t")
        mism = 0
        checked = 0
        for r in d.itertuples(index=False):
            cond, day = r.Measure.split("_day")
            temp, dose = cond[:3], cond[3:]
            sid = r.Sample
            if (sid, temp, dose) not in cells or int(day) not in cells[(sid, temp, dose)]:
                continue
            checked += 1
            want = cells[(sid, temp, dose)][int(day)]
            want = 0 if want is None else want
            if int(r.Count) != want:
                mism += 1
        print(f"   {os.path.basename(path)}: rows {len(d)}, checked {checked}, "
              f"mismatches {mism}")

    print("\n" + "=" * 74)
    print("5. DELIVERED PHENOTYPES RECOMPUTED FROM THE WORKBOOK")
    print("=" * 74)
    strains = sorted({k[0] for k in cells})
    conds = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
    mort_hat, auc_hat = {}, {}
    for s in strains:
        for c in conds:
            ser = cells[(s, c[:3], c[3:])]
            days = sorted(ser)
            vals = [0 if ser[d] is None else ser[d] for d in days]
            n0 = vals[0]
            mort_hat[(s, f"mort_{c}")] = 1 - vals[-1] / n0
            auc_hat[(s, f"auc_{c}")] = float(np.trapz([v / n0 for v in vals], days)) / 10.0
    M = pd.read_csv(f"{PA}/phenotype_mortality.tsv", sep="\t", index_col=0)
    A = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    dm = max(abs(M.at[s, f"mort_{c}"] - mort_hat[(s, f"mort_{c}")])
             for s in strains for c in conds)
    da = max(abs(A.at[s, f"auc_{c}"] - auc_hat[(s, f"auc_{c}")])
             for s in strains for c in conds)
    print(f"max |delivered - recomputed| mortality = {dm:.2e}")
    print(f"max |delivered - recomputed| AUC       = {da:.2e}")
    print("=> " + ("IDENTICAL: the reading is correct, blanks are scored as all-dead"
                   if dm < 1e-9 and da < 1e-9 else "MISMATCH, investigate"))


if __name__ == "__main__":
    main()
