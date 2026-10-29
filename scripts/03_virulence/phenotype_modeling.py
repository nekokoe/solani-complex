#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 1: phenotype modelling from raw survival data.

Input: 07.lethal_groups/caterpillar_37C_narrow.txt
       21 strains x 6 conditions (28Ce3/e4/e5, 37Ce3/e4/e5) x day0-10 survival
       (regenerated 2026-10-12 from 21株-大蜡螟原始数据.xlsx by
        scripts/convert_galleria_xlsx.py; the 15-strain original is preserved at
        07.lethal_groups/.bak-21strain/caterpillar_37C_narrow.txt)
Outputs:
  - phenotype_mortality.tsv    21 strains x 6 conditions (day10 mortality)
  - phenotype_auc.tsv          21 strains x 6 conditions (survival AUC)
  - phenotype_lt50.tsv         21 strains x 6 conditions (LT50, days)
  - phenotype_metrics.tsv      9 traits (6 mort + mort_28/37_mean + temp_effect)
  - figure_phenotype.pdf       survival curves + heatmap + temperature split
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW = "07.lethal_groups/caterpillar_37C_narrow.txt"
OUT = "07.lethal_groups/phenotype_analysis"
import os
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(RAW, sep="\t")
d[["Condition", "Day"]] = d.Measure.str.extract(r"(\d+C[eE]\d+)_day(\d+)")
d["Day"] = d.Day.astype(int)
print("conditions:", sorted(d.Condition.unique()))
print("strains:", d.Sample.nunique())

# survival rate per strain x condition x day
surv = d.pivot_table(index=["Sample", "Condition"], columns="Day",
                     values="Percent", aggfunc="first")

# ---- phenotype metrics ----
def auc10(row):
    """Area under survival curve day0-10 (0=all dead, 10=max survival)."""
    xs = sorted(row.dropna().index)
    ys = [row[x] for x in xs]
    return float(np.trapz(ys, xs)) / 10.0

pheno = pd.DataFrame(index=d.Sample.unique())
for cond in ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]:
    sub = surv.loc[(slice(None), cond), :]
    sub = sub.droplevel("Condition")
    mortality = 1 - sub.get(10, np.nan)
    pheno[f"mort_{cond}"] = mortality.reindex(pheno.index)
auc_df = pd.DataFrame(index=d.Sample.unique())
for cond in ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]:
    sub = surv.loc[(slice(None), cond), :].droplevel("Condition")
    auc_df[f"auc_{cond}"] = sub.apply(auc10, axis=1).reindex(auc_df.index)

pheno.to_csv(f"{OUT}/phenotype_mortality.tsv", sep="\t")
auc_df.to_csv(f"{OUT}/phenotype_auc.tsv", sep="\t")

# ---- LT50 (days until 50% of the cohort is dead) ----
# Not previously generated anywhere in the repository, although the manuscript
# quotes per-strain LT50 values. Linear interpolation between the last day at or
# above 0.5 survival and the first day below it; a cohort that never drops below
# 0.5 is reported as >10 (censored) via NaN plus a flag column.
def lt50(row):
    xs = sorted(row.dropna().index)
    prev = None
    for x in xs:
        s = row[x]
        if s <= 0.5:
            if prev is None:
                return float(x)                      # already <=0.5 on day 0
            x0, s0 = prev
            if s0 == s:
                return float(x)
            return x0 + (s0 - 0.5) * (x - x0) / (s0 - s)
        prev = (x, s)
    return np.nan

lt = pd.DataFrame(index=d.Sample.unique())
for cond in ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]:
    sub = surv.loc[(slice(None), cond), :].droplevel("Condition")
    lt[f"lt50_{cond}"] = sub.apply(lt50, axis=1).reindex(lt.index)
lt.to_csv(f"{OUT}/phenotype_lt50.tsv", sep="\t")
lt_all = lt.mean(axis=1)
print(f"\nLT50 (days, mean over the 6 conditions): "
      f"{lt_all.min():.2f} ({lt_all.idxmin()}) .. {lt_all.max():.2f} ({lt_all.idxmax()}); "
      f"median {lt_all.median():.2f}; censored (never <50% dead) {int(lt_all.isna().sum())}")

# ---- effects ----
rows = []
# temperature effect per strain (mean over concentrations, 28 vs 37)
m28 = pheno[["mort_28Ce3", "mort_28Ce4", "mort_28Ce5"]].mean(axis=1)
m37 = pheno[["mort_37Ce3", "mort_37Ce4", "mort_37Ce5"]].mean(axis=1)
pheno["mort_28_mean"] = m28
pheno["mort_37_mean"] = m37
pheno["temp_effect"] = m28 - m37
# (concentration slope traits removed by decision 2026-08-13)

pheno.to_csv(f"{OUT}/phenotype_metrics.tsv", sep="\t")

w_tt = stats.wilcoxon(m28, m37)
print(f"\ntemperature effect: 28C mort {m28.mean():.2f} vs 37C {m37.mean():.2f}, "
      f"paired Wilcoxon p={w_tt.pvalue:.4f}")
for c in ["Ce3", "Ce4", "Ce5"]:
    a = pheno[f"mort_28{c}"].dropna(); b = pheno[f"mort_37{c}"].dropna()
    p = stats.wilcoxon(a, b).pvalue
    print(f"  {c}: 28C {a.mean():.2f} vs 37C {b.mean():.2f} (p={p:.4f})")
fr = stats.f_oneway(*[pheno[f"mort_28{c}"].dropna() for c in ["Ce3", "Ce4", "Ce5"]])
fr37 = stats.f_oneway(*[pheno[f"mort_37{c}"].dropna() for c in ["Ce3", "Ce4", "Ce5"]])
print(f"concentration effect: 28C one-way p={fr.pvalue:.2e}; 37C p={fr37.pvalue:.2e}")
# strain x temperature interaction (Friedman across conditions)
print("\nstrain x temperature interaction candidates (|temp effect| top):")
for s, te in pheno.temp_effect.sort_values(key=abs, ascending=False).items():
    print(f"  {s}: 28C {pheno.loc[s,'mort_28_mean']:.2f} vs 37C {pheno.loc[s,'mort_37_mean']:.2f} (diff {te:+.2f})")

# ---- figure ----
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6), dpi=300)
# a: survival curves (37C e4 example)
ax = axes[0]
for s in pheno.index:
    sub = surv.loc[(s, "37Ce4")]
    ax.plot(sorted(sub.index), sub.reindex(sorted(sub.index)), "o-", ms=2.5, lw=0.8,
            alpha=0.6, color="grey")
ax.set_xlabel("Day"); ax.set_ylabel("Survival")
ax.set_title(f"a  37C e4 survival curves ({len(pheno)} strains)")
# b: mortality heatmap strain x condition
ax = axes[1]
M = pheno[[f"mort_{c}" for c in ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]]]
im = ax.imshow(M.values, cmap="RdYlGn_r", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(6))
ax.set_xticklabels(["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"],
                   rotation=45, ha="right", fontsize=6)
ax.set_yticks(range(len(M)))
ax.set_yticklabels(M.index, fontsize=5.5)
ax.set_title("b  Mortality by condition")
fig.colorbar(im, ax=ax, fraction=0.046)
# c: temperature split
ax = axes[2]
for s in pheno.index:
    ax.plot([0, 1], [pheno.loc[s, "mort_28_mean"], pheno.loc[s, "mort_37_mean"]],
            "o-", color="grey", ms=3, lw=0.7, alpha=0.6)
ax.set_xticks([0, 1]); ax.set_xticklabels(["28C", "37C"])
ax.set_ylabel("Mean mortality")
ax.set_title("c  Temperature effect per strain")
for s, te in pheno.temp_effect.sort_values(key=abs, ascending=False).head(4).items():
    col = "#C44E52" if te > 0 else "#3B7DD8"
    ax.plot([0, 1], [pheno.loc[s, "mort_28_mean"], pheno.loc[s, "mort_37_mean"]],
            "o-", color=col, ms=4, lw=1.2, label=f"{s} ({te:+.2f})")
ax.legend(fontsize=5.5, frameon=False)
for s in ["top", "right"]:
    for axx in axes:
        axx.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/figure_phenotype.pdf")
fig.savefig(f"{OUT}/figure_phenotype.png")
print("\ndone ->", OUT)
