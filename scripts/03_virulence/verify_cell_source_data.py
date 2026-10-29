#!/usr/bin/env python3
"""Audit the claims made about the 2026-10-21 cell source data.

`261020细胞蜡螟结果_and_like方法修改意见/Source data-cell and larvae.xlsx`, sheet `cell`,
holds every assay as 23 groups of THREE ADJACENT COLUMNS (Ctrl, LPS, then the 21
strains) and two rows per assay: `High concentration` and `Low concentration`. The
replicates are therefore columns here — unlike the 2026-10-12 workbook
(`21株-细胞原始数据.xlsx`), where the three replicates are rows. Reading the wrong
axis silently turns 3 replicates into 3 samples, so this script parses the sheet
itself and then recomputes every statistic the Results paragraph quotes.

Blocks present
    PBMC-derived : TNFa qPCR + ELISA, IL1b qPCR + ELISA, csf2 qPCR, cxcl1 qPCR
    THP1-derived : IL1b qPCR + ELISA, TNFa qPCR + ELISA, csf2 qPCR, cxcl1 qPCR
Only IL-1b and TNF-a have ELISA, which is what the text states.
The low-concentration rows carry no Ctrl/LPS (6 of 69 values); the 21 strains are complete.
"""

import os
import re

import numpy as np
import openpyxl
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.path.join(ROOT, "261020细胞蜡螟结果_and_like方法修改意见",
                    "Source data-cell and larvae.xlsx")
GROUPS = os.path.join(ROOT, "bmu_ref_groups.211.txt")
CACHE = os.path.join(ROOT, "07.lethal_groups/cell_assay/source_2026_10_21_long.tsv")


def parse():
    ws = openpyxl.load_workbook(XLSX, data_only=True)["cell"]
    model = None
    rows = []
    header, cyt, assay = {}, None, None
    for r in range(1, ws.max_row + 1):
        a = ws.cell(r, 1).value
        if a is None:
            continue
        a = str(a).strip()
        if a in ("PBMCs-derived macrophage", "THP1-derived macrophage"):
            model = "PBMC" if a.startswith("PBMC") else "THP1"
        elif a in ("TNFa", "IL1b", "csf2", "cxcl1"):
            cyt = a
        elif a in ("qPCR", "ELISA"):
            assay = a
            header = {}
            for c in range(2, ws.max_column + 1):
                v = ws.cell(r, c).value
                if v not in (None, ""):
                    header[c] = str(v).strip()
        elif a in ("High concentration", "Low concentration"):
            conc = "high" if a.startswith("High") else "low"
            for c, name in header.items():
                vals = [ws.cell(r, c + k).value for k in range(3)]
                vals = [float(v) for v in vals if isinstance(v, (int, float))]
                if not vals:
                    continue
                rows.append(dict(model=model, cytokine=cyt, assay=assay,
                                 concentration=conc, group=name,
                                 rep=len(vals), mean=float(np.mean(vals)),
                                 sd=float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan))
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    df.to_csv(CACHE, sep="\t", index=False)
    return df


def audit(df):
    print("=" * 78)
    print("0. STRUCTURE")
    print("=" * 78)
    n_rep = df.groupby(["model", "cytokine", "assay", "concentration", "group"]).rep.max()
    print(f"assay x concentration x group combinations: {len(n_rep)}")
    print(f"replicates per group (should be 3): {sorted(set(n_rep))}")
    miss = df[(df.concentration == "low") & (df.group.isin(["Ctrl", "LPS"]))]
    print(f"low-concentration control rows present: {len(miss)} (expected 0 -> Ctrl/LPS "
          f"exist only in the high block)")
    print(f"ELISA available for: {sorted(df[df.assay == 'ELISA'].cytokine.unique())}")
    print(f"qPCR available for : {sorted(df[df.assay == 'qPCR'].cytokine.unique())}")

    strains = sorted(g for g in df.group.unique() if g.startswith("BMU"))
    print(f"strains: {len(strains)}")

    print("\n" + "=" * 78)
    print("1. HIGH vs LOW, PER MODEL x MARKER  (text: Wilcoxon P = 9.5e-7..5.9e-5,")
    print("   median high/low ratio 1.18-1.67)")
    print("=" * 78)
    print(f"{'model':6s} {'marker':6s} {'assay':6s} {'P':>10s} {'median ratio':>13s}")
    ps, ratios = [], []
    for (m, c), g in df[df.assay == "qPCR"].groupby(["model", "cytokine"]):
        hi = g[g.concentration == "high"].set_index("group")["mean"]
        lo = g[g.concentration == "low"].set_index("group")["mean"]
        both = [s for s in strains if s in hi.index and s in lo.index]
        p = stats.wilcoxon(hi[both], lo[both]).pvalue
        ratio = float(np.median(hi[both] / lo[both]))
        ps.append(p)
        ratios.append(ratio)
        print(f"{m:6s} {c:6s} {'qPCR':6s} {p:10.2e} {ratio:13.2f}")
    print(f"  -> P range {min(ps):.1e}..{max(ps):.1e}; median-ratio range "
          f"{min(ratios):.2f}..{max(ratios):.2f}")

    print("\n" + "=" * 78)
    print("2. LPS FOLD INDUCTION OVER CTRL  (text: TNF-a 31x THP-1, 149x PBMC)")
    print("=" * 78)
    for m in ("THP1", "PBMC"):
        g = df[(df.model == m) & (df.cytokine == "TNFa") & (df.assay == "qPCR")
               & (df.concentration == "high")].set_index("group")["mean"]
        print(f"  {m:5s} TNFa qPCR: Ctrl {g['Ctrl']:.2f}, LPS {g['LPS']:.1f} "
              f"-> {g['LPS'] / g['Ctrl']:.1f}-fold")

    print("\n" + "=" * 78)
    print("3. MEDIAN ISOLATE FOLD OVER CTRL  (text: 8x THP-1, 41x PBMC)")
    print("=" * 78)
    for m in ("THP1", "PBMC"):
        g = df[(df.model == m) & (df.cytokine == "TNFa") & (df.assay == "qPCR")
               & (df.concentration == "high")].set_index("group")["mean"]
        fold = np.median([g[s] / g["Ctrl"] for s in strains])
        print(f"  {m:5s} median TNFa fold over Ctrl = {fold:.1f}")

    print("\n" + "=" * 78)
    print("4. TOP RESPONDERS  (text: THP-1 = 12684/12685/01936; PBMC = 12684/05572/12685)")
    print("=" * 78)
    for m in ("THP1", "PBMC"):
        for c in ("TNFa", "IL1b", "csf2", "cxcl1"):
            g = df[(df.model == m) & (df.cytokine == c) & (df.assay == "qPCR")
                   & (df.concentration == "high") & (df.group.isin(strains))]
            top = g.nlargest(3, "mean").group.tolist()
            print(f"  {m:5s} {c:6s} top3 = {', '.join(x.replace('BMU','') for x in top)}")

    print("\n" + "=" * 78)
    print("5. CLINICAL vs ENVIRONMENTAL  (text: Mann-Whitney P = 0.08-1.00, 12 combos)")
    print("=" * 78)
    src = pd.read_csv(GROUPS, sep="\t")
    env = dict(zip(src.BMU_ID.astype(str), src.MedEnv.astype(str)))
    med = [s for s in strains if env.get(s) == "med"]
    envs = [s for s in strains if env.get(s) == "env"]
    print(f"  clinical {len(med)} vs environmental {len(envs)}")
    out = []
    for (m, c, a), g in df[(df.concentration == "high") & df.group.isin(strains)].groupby(
            ["model", "cytokine", "assay"]):
        s = g.set_index("group")["mean"]
        p = stats.mannwhitneyu([s[x] for x in med], [s[x] for x in envs]).pvalue
        out.append((m, c, a, p))
    for m, c, a, p in sorted(out, key=lambda x: x[3]):
        print(f"  {m:5s} {c:6s} {a:6s} P = {p:.3f}")
    print(f"  -> {len(out)} combinations, P range {min(o[3] for o in out):.3f}.."
          f"{max(o[3] for o in out):.3f}")

    print("\n" + "=" * 78)
    print("6. CROSS-MODEL CONCORDANCE  (text: TNFa 0.40/0.07, IL1b 0.49/0.024,")
    print("   CXCL1 0.66/0.001, CSF2 0.03/0.91)")
    print("=" * 78)
    for c in ("TNFa", "IL1b", "cxcl1", "csf2"):
        a = df[(df.model == "THP1") & (df.cytokine == c) & (df.assay == "qPCR")
               & (df.concentration == "high")].set_index("group")["mean"].reindex(strains)
        b = df[(df.model == "PBMC") & (df.cytokine == c) & (df.assay == "qPCR")
               & (df.concentration == "high")].set_index("group")["mean"].reindex(strains)
        r = stats.spearmanr(a, b)
        print(f"  {c:6s} rho = {r.statistic:.2f}  P = {r.pvalue:.3f}")


if __name__ == "__main__":
    d = parse() if not os.path.exists(CACHE) else pd.read_csv(CACHE, sep="\t")
    audit(d)
    print(f"\ncache: {CACHE}")
