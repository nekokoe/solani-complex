#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phenotype_stats.py — cohort tests and the Table-1 statistics for the
21-strain Galleria experiment.

Nothing in the repository generated the published Table 1 or the
clinical-versus-environmental comparison, so both are reconstructed here from
the phenotype tables produced by scripts/phenotype_modeling.py.

Outputs (07.lethal_groups/phenotype_analysis/):
  cohort_tests.tsv        clinical (17) vs environmental (4), Mann-Whitney per trait
  table1_stats.tsv        dose / temperature / strain effects + dose-LT50 correlation
  (temperature_bias.tsv is written by cmh_temperature_preference.py, not here)

Cohort labels come from 05.pan-genome/bmu_ref_groups.211.txt (MedEnv/Host/Species).
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "07.lethal_groups", "phenotype_analysis")
META = os.path.join(ROOT, "05.pan-genome", "bmu_ref_groups.211.txt")
CONDS = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]

metrics = pd.read_csv(f"{OUT}/phenotype_metrics.tsv", sep="\t", index_col=0)
lt50 = pd.read_csv(f"{OUT}/phenotype_lt50.tsv", sep="\t", index_col=0)
auc = pd.read_csv(f"{OUT}/phenotype_auc.tsv", sep="\t", index_col=0)
mort = pd.read_csv(f"{OUT}/phenotype_mortality.tsv", sep="\t", index_col=0)

meta = pd.read_csv(META, sep="\t")
lab = meta.set_index("BMU_ID").reindex(metrics.index)
clinical = lab.index[lab.MedEnv == "med"]
envir = lab.index[lab.MedEnv != "med"]
print(f"cohort: {len(clinical)} clinical (Human) + {len(envir)} environmental "
      f"({sorted(set(lab.Host[lab.MedEnv != 'med']))})")

# ---------------------------------------------------------------- cohort tests
# Two 13-trait sets, because the published text says "13 traits" without listing
# them and both readings give 13. Reporting both also shows whether the
# conclusion depends on the choice.
def trait_set(name, use_auc=False, use_lt50=False):
    # NB: parameter names must not shadow the module-level frames (mort/auc/lt50)
    cols = {f"mort_{c}": mort[f"mort_{c}"] for c in CONDS}
    if use_auc:
        cols.update({f"auc_{c}": auc[f"auc_{c}"] for c in CONDS})
    if use_lt50:
        cols.update({f"lt50_{c}": lt50[f"lt50_{c}"] for c in CONDS})
    cols["temp_effect"] = metrics["temp_effect"]
    return name, pd.DataFrame(cols)


sets = [trait_set("mort+AUC+temp_effect", use_auc=True),
        trait_set("mort+LT50+temp_effect", use_lt50=True)]

rows, done = [], []
for setname, T in sets:
    for col in T.columns:
        a = T[col].reindex(clinical).dropna()
        b = T[col].reindex(envir).dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"trait_set": setname, "trait": col,
                     "n_clinical": len(a), "n_env": len(b),
                     "median_clinical": a.median(), "median_env": b.median(),
                     "mean_clinical": a.mean(), "mean_env": b.mean(),
                     "U": u, "p": p,
                     "clinical_higher": a.mean() > b.mean()})
    sub = [r for r in rows if r["trait_set"] == setname]
    ps = [r["p"] for r in sub]
    done.append((setname, len(sub), min(ps), max(ps),
                 sum(r["clinical_higher"] for r in sub)))
ct = pd.DataFrame(rows)
ct.to_csv(f"{OUT}/cohort_tests.tsv", sep="\t", index=False)

print("\n=== clinical vs environmental (Mann-Whitney), two 13-trait readings ===")
for setname, n, pmin, pmax, nhigher in done:
    print(f"  {setname:22s} traits={n}  p range {pmin:.3f}-{pmax:.3f}  "
          f"clinical higher in {nhigher}/{n}")

# ------------------------------------------------------------------- Table 1
# long form: one row per strain x condition
long = (mort.reset_index().melt(id_vars="index", var_name="condition", value_name="mort")
        .rename(columns={"index": "strain"}))
# NB the melted column names carry the 'mort_' prefix, so the temperature cannot be
# read off the leading characters - extract it with a pattern that skips the prefix.
long[["temp_n", "dose"]] = long.condition.str.extract(r"(28|37)C(e\d)")
long["temp"] = long.temp_n + "C"
assert long.temp.nunique() == 2, f"expected 2 temperatures, got {sorted(set(long.temp))}"
assert long.dose.nunique() == 3, f"expected 3 doses, got {sorted(set(long.dose))}"

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

model = ols("mort ~ C(dose) + C(temp) + C(dose):C(temp)", data=long).fit()
aov = sm.stats.anova_lm(model, typ=2)

t1 = []
for term, r in aov.iterrows():
    label = {"C(dose)": "Main effect of dose",
             "C(temp)": "Main effect of temperature",
             "C(dose):C(temp)": "Dose x temperature interaction",
             "Residual": "Residual"}.get(term, term)
    t1.append({"term": label, "df": r["df"], "F": r["F"], "p": r["PR(>F)"],
               "note": ""})

# strain main effect (one-way across strains) and across-dose Tukey
fr = stats.f_oneway(*[g.mort.values for _, g in long.groupby("strain")])
t1.append({"term": "Main effect of strain (one-way ANOVA)", "df": len(long.strain.unique()) - 1,
           "F": fr.statistic, "p": fr.pvalue, "note": ""})
tk_dose = pairwise_tukeyhsd(long.mort, long.dose)
t1.append({"term": "Multiple comparisons of dose (Tukey HSD)", "df": np.nan, "F": np.nan,
           "p": tk_dose.pvalues.min(),
           "note": f"max of {len(tk_dose.pvalues)} pairwise p; all significant={bool(tk_dose.reject.all())}"})
tk_temp = pairwise_tukeyhsd(long.mort, long.temp)
t1.append({"term": "Multiple comparisons of temperature (Tukey HSD)", "df": np.nan, "F": np.nan,
           "p": float(tk_temp.pvalues[0]),
           "note": f"28C vs 37C; significant={bool(tk_temp.reject.all())}"})

# dose vs LT50: one point per strain x dose (mean over the two temperatures)
lt_long = (lt50.reset_index().melt(id_vars="index", var_name="condition", value_name="lt50")
           .rename(columns={"index": "strain"}))
lt_long["dose"] = lt_long.condition.str.extract(r"(e\d)")[0]
lt_long["dose_n"] = lt_long.dose.map({"e3": 1e3, "e4": 1e4, "e5": 1e5})
per_dose = lt_long.dropna(subset=["lt50"])
for method, fn in [("Pearson", stats.pearsonr), ("Spearman", stats.spearmanr)]:
    res = fn(per_dose.dose_n, per_dose.lt50)
    t1.append({"term": f"Correlation: dose vs LT50 ({method})", "df": len(per_dose) - 2,
               "F": np.nan, "p": res.pvalue,
               "note": f"r = {res.statistic:.3f}, n = {len(per_dose)} (strain x dose, "
                       f"mean of the two temperatures); censored cells dropped"})

t1df = pd.DataFrame(t1)
t1df.to_csv(f"{OUT}/table1_stats.tsv", sep="\t", index=False)
print("\n=== reconstructed Table 1 ===")
print(t1df.to_string(index=False))
print("\nper-dose mean LT50:", per_dose.groupby("dose").lt50.mean().round(2).to_dict())

# ------------------------------------------------------------ temperature bias
te = metrics["temp_effect"].sort_values()
# The temperature-bias classification is NOT produced here. A cutoff on temp_effect
# (an absolute threshold on temp_effect) was retired: the preference is a comparison of
# binomial counts across three
# inoculum doses, so scripts/cmh_temperature_preference.py owns temperature_bias.tsv
# (Cochran-Mantel-Haenszel stratified by dose). This script only reports the effect size.
print(f"\n=== temperature effect (continuous, not a classification) ===")
print(f"  range: {te.min():+.3f} to {te.max():+.3f} (n = {len(te)})")
print("  classification: run scripts/cmh_temperature_preference.py")
print(f"\nwrote cohort_tests.tsv, table1_stats.tsv -> {os.path.relpath(OUT, ROOT)}")
