#!/usr/bin/env python3
"""How much does the phenotype definition change the answer? mortality vs survival AUC.

Three layers, each answered with numbers rather than opinion.

  1 phenotype   per condition: how well do the two phenotypes agree, and where does the
                day-10 endpoint lose information the trajectory keeps?
  2 association per matched trait pair: does the HOG screen pick the same loci? hits,
                q minima, effective independent signals after de-collinearisation,
                leave-one-strain-out stability of each screen's top hit.
  3 preference  does the strain-level temperature preference survive the change of
                phenotype, and which strains move?

Output  07.lethal_groups/genomewide_assoc/mortality_vs_auc_comparison.tsv
"""

import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GW = os.path.join(BASE, "07.lethal_groups/genomewide_assoc")
PA = os.path.join(BASE, "07.lethal_groups/phenotype_analysis")
HOGCOUNT = os.path.join(BASE, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
CONDS = ["28Ce3", "28Ce4", "28Ce5", "37Ce3", "37Ce4", "37Ce5"]
TIE_MORT, GAIN_AUC = 0.05, 0.10


def pheno_frames():
    m = pd.read_csv(f"{PA}/phenotype_mortality.tsv", sep="\t", index_col=0)
    a = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    vir = pd.DataFrame({f"vir_{c}": 1.0 - a[f"auc_{c}"] for c in CONDS})
    return m, a, vir


def layer1(m, a, vir, rows):
    print("=" * 78)
    print("1. PHENOTYPE LAYER")
    print("=" * 78)
    print(f"{'condition':10s} {'rho(mort,vir)':>14s} {'IQR mort':>9s} {'IQR vir':>8s} "
          f"{'uniq mort':>9s} {'uniq vir':>8s} {'pairs lost':>10s}")
    for c in CONDS:
        mo, vi = m[f"mort_{c}"], vir[f"vir_{c}"]
        rho = stats.spearmanr(mo, vi).statistic
        # strain pairs the endpoint cannot separate but the trajectory can
        lost = 0
        names = list(mo.index)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                if abs(mo[names[i]] - mo[names[j]]) <= TIE_MORT and \
                        abs(vi[names[i]] - vi[names[j]]) >= GAIN_AUC:
                    lost += 1
        print(f"{c:10s} {rho:>14.3f} {mo.quantile(.75)-mo.quantile(.25):>9.3f} "
              f"{vi.quantile(.75)-vi.quantile(.25):>8.3f} {mo.nunique():>9d} "
              f"{vi.nunique():>8d} {lost:>10d}")
        rows.append(dict(layer="phenotype", item=c, metric="rho_mort_vs_vir_auc",
                         value=round(rho, 4)))
        rows.append(dict(layer="phenotype", item=c, metric="pairs_indistinguishable_at_day10",
                         value=lost))


def load_hog(strains):
    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = (hog[hog.SampleID.isin(strains)].pivot(index="HOG_ID", columns="SampleID",
                                                 values="GeneCount")
           .reindex(columns=strains).fillna(0))
    n_pres = (mat > 0).sum(axis=1)
    var = mat[strains].nunique(axis=1) > 1
    return mat[(n_pres >= 3) & (n_pres <= len(strains) - 2) & var]


def decollinearize(df, mat, rcut=0.8):
    """Greedy: keep the best-q HOG, drop anything |r|>=rcut with something already kept."""
    kept = []
    for h in df.sort_values("q").HOG_ID:
        v = mat.loc[h].astype(float).values
        if all(abs(np.corrcoef(v, mat.loc[k].astype(float).values)[0, 1]) < rcut
               for k in kept):
            kept.append(h)
    return kept


def layer2(rows):
    mort = pd.read_csv(f"{GW}/hog_phenotype_assoc.tsv", sep="\t")
    auc = pd.read_csv(f"{GW}/hog_phenotype_assoc_auc.tsv", sep="\t")
    strains = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0).index.tolist()
    mat = load_hog(strains)

    print("\n" + "=" * 78)
    print("2. ASSOCIATION LAYER (HOG screen)")
    print("=" * 78)
    mt, at = mort[mort.tested], auc
    print(f"{'':26s} {'mortality':>12s} {'AUC':>12s}")
    for lbl, fn in (("rows in table", len), ("tests entering FDR", len),
                    ("min q", lambda d: d.q.min()), ("q<0.05", lambda d: int((d.q < 0.05).sum())),
                    ("q<0.10", lambda d: int((d.q < 0.10).sum())),
                    ("p<=1e-3", lambda d: int((d.p <= 1e-3).sum()))):
        a_ = len( mort) if lbl == "rows in table" else fn(mt)
        b_ = len(auc) if lbl == "rows in table" else fn(at)
        print(f"{lbl:26s} {a_:>12} {b_:>12}")
        rows.append(dict(layer="association", item=lbl, metric="mortality", value=a_))
        rows.append(dict(layer="association", item=lbl, metric="auc", value=b_))

    print("\nper matched trait pair: agreement of the 4,484 rho vectors")
    print(f"{'mortality trait':16s} {'AUC trait':16s} {'rho-of-rhos':>12s} {'same top-20':>12s}")
    pairs = [(f"mort_{c}", f"vir_{c}") for c in CONDS] + \
            [("mort_28_mean", "vir_28_mean"), ("mort_37_mean", "vir_37_mean"),
             ("temp_effect", "vir_temp_diff")]
    for mt_, at_ in pairs:
        x = mort[mort.trait == mt_].set_index("HOG_ID").rho
        y = auc[auc.trait == at_].set_index("HOG_ID").rho
        j = pd.concat([x, y], axis=1, keys=["m", "a"]).dropna()
        r = stats.spearmanr(j.m, j.a).statistic
        t20m = set(j.nsmallest(20, "m").index)
        t20a = set(j.nsmallest(20, "a").index)
        print(f"{mt_:16s} {at_:16s} {r:>12.3f} {len(t20m & t20a):>9d}/20")
        rows.append(dict(layer="association", item=f"{mt_}~{at_}", metric="rho_of_rhos",
                         value=round(r, 4)))
        rows.append(dict(layer="association", item=f"{mt_}~{at_}",
                         metric="top20_overlap", value=len(t20m & t20a)))

    print("\nhit sets at q<0.10, de-collinearised at |r|>=0.8")
    for name, d in (("mortality", mt[mt.q < 0.1]), ("AUC", at[at.q < 0.1])):
        if len(d):
            kept = decollinearize(d, mat)
            print(f"  {name:10s} {len(d):3d} rows -> {len(kept):2d} independent signals")
            rows.append(dict(layer="association", item=name, metric="hits_q10", value=len(d)))
            rows.append(dict(layer="association", item=name, metric="independent_signals",
                             value=len(kept)))
        else:
            print(f"  {name:10s}   0 rows ->  0 independent signals")
            rows.append(dict(layer="association", item=name, metric="hits_q10", value=0))

    print("\nleave-one-strain-out stability of each screen's top hit")
    ph = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
    ap = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    d_mort = ph["temp_effect"].reindex(strains)
    d_auc = ((1 - ap[[f"auc_{c}" for c in CONDS[:3]]].mean(1))
             - (1 - ap[[f"auc_{c}" for c in CONDS[3:]]].mean(1))).reindex(strains)

    def trait_vector(kind, trait):
        if kind == "mortality":
            if trait == "temp_effect":
                return d_mort.values
            if trait.startswith("mort_"):
                return ph[trait].reindex(strains).values
        else:
            if trait == "vir_temp_diff":
                return d_auc.values
            if trait.startswith("vir_"):
                return (1 - ap[f"auc_{trait[len('vir_'):]}"].reindex(strains).values)
        raise KeyError(trait)

    for name, d in (("mortality", mt), ("AUC", at)):
        top = d.nsmallest(1, "p")
        h, tr = top.HOG_ID.iloc[0], top.trait.iloc[0]
        y = trait_vector(name, tr)
        x = mat.loc[h].astype(float).values
        base = stats.spearmanr(x, y).statistic
        loo = [stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic
               for i in range(len(y))]
        print(f"  {name:10s} {h} vs {tr}: rho {base:+.3f}  LOO {min(loo):+.3f}..{max(loo):+.3f}")
        rows.append(dict(layer="association", item=f"{name}:{h}~{tr}", metric="rho",
                         value=round(base, 4)))
        rows.append(dict(layer="association", item=f"{name}:{h}~{tr}", metric="loo_abs_min",
                         value=round(min(abs(v) for v in loo), 4)))


def layer3(rows):
    ph = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
    ap = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    bias = pd.read_csv(f"{PA}/temperature_bias_cmh.tsv", sep="\t").set_index("strain")
    d_mort = ph["temp_effect"]
    d_auc = (1 - ap[[f"auc_{c}" for c in CONDS[:3]]].mean(1)) - \
            (1 - ap[[f"auc_{c}" for c in CONDS[3:]]].mean(1))
    # same sign convention as temp_effect: positive = 28 °C more lethal
    tab = pd.DataFrame({"d_mort": d_mort, "d_auc": d_auc}).dropna()
    r = stats.spearmanr(tab.d_mort, tab.d_auc).statistic
    print("\n" + "=" * 78)
    print("3. PREFERENCE LAYER")
    print("=" * 78)
    print(f"Spearman(d_mort, d_auc) = {r:.3f}")
    rows.append(dict(layer="preference", item="all", metric="rho_dmort_vs_dauc",
                     value=round(r, 4)))
    tab = tab.join(bias[["class", "P"]])
    tab["rank_mort"] = tab.d_mort.rank(ascending=False).astype(int)
    tab["rank_auc"] = tab.d_auc.rank(ascending=False).astype(int)
    tab["rank_shift"] = (tab.rank_mort - tab.rank_auc).abs()
    print("\nlargest rank shifts (mortality ranking vs AUC ranking, 1 = most 28 °C-biased):")
    print(tab.sort_values("rank_shift", ascending=False)
          .head(8)[["d_mort", "d_auc", "rank_mort", "rank_auc", "rank_shift", "class", "P"]]
          .round(4).to_string())
    for s, r_ in tab.sort_values("rank_shift", ascending=False).head(5).iterrows():
        rows.append(dict(layer="preference", item=s, metric="rank_shift",
                         value=int(r_.rank_shift)))
    n_same_sign = int((np.sign(tab.d_mort) == np.sign(tab.d_auc)).sum())
    print(f"\nsign agreement of the per-strain preference: {n_same_sign}/{len(tab)} strains")
    rows.append(dict(layer="preference", item="all", metric="sign_agreement",
                     value=f"{n_same_sign}/{len(tab)}"))



def _partial_spearman(x, y, Z):
    """Rank-based partial correlation of x and y controlling for the columns of Z."""
    rx, ry = rankdata(x), rankdata(y)
    A = np.column_stack([np.ones(len(y))] + [rankdata(Z[:, i]) for i in range(Z.shape[1])])
    def resid(v):
        b, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ b
    ex, ey = resid(rx), resid(ry)
    r = np.corrcoef(ex, ey)[0, 1]
    k = A.shape[1] - 1
    t = r * np.sqrt((len(y) - 2 - k) / (1 - r ** 2)) if abs(r) < 1 else np.inf
    return r, 2 * stats.t.sf(abs(t), df=len(y) - 2 - k)


def snp_pcs(strains, cache, n=3):
    """Genotype principal components, used as the population-structure covariates."""
    g = np.load(os.path.join(cache, "snp21.npy")).astype(np.float32)
    lab = [str(x) for x in np.load(os.path.join(cache, "snpdist_labels.npy"), allow_pickle=True)]
    g = g[:, [lab.index(s) for s in strains]]
    comp = (g >= 0).all(axis=1)
    p = g.mean(axis=1) / 2.0
    maf = np.minimum(p, 1 - p)
    g = g[comp & (maf >= 0.10)]
    X = (g - g.mean(0)) / np.where(g.std(0) == 0, 1, g.std(0))
    _, sv, Vt = np.linalg.svd(X, full_matrices=False)
    return (Vt[:n].T * sv[:n]), (sv[:5] ** 2 / np.sum(sv ** 2) * 100)


def layer4(rows, mat, strains, ph, ap):
    """Does the association survive control for population structure?"""
    print("\n" + "=" * 78)
    print("4. STRUCTURE CONTROL (rank partial correlation on SNP PC1-PC2)")
    print("=" * 78)
    pcs, var = snp_pcs(strains, os.path.join(GW, "snp_cache"))
    print(f"SNP PC variance explained (%): {np.round(var, 1)}")
    auc = pd.read_csv(f"{GW}/hog_phenotype_assoc_auc.tsv", sep="\t")
    mort = pd.read_csv(f"{GW}/hog_phenotype_assoc.tsv", sep="\t")
    mort = mort[mort.tested]
    d_auc = ((1 - ap[[f"auc_{c}" for c in CONDS[:3]]].mean(1))
             - (1 - ap[[f"auc_{c}" for c in CONDS[3:]]].mean(1)))

    def yvec(trait):
        if trait == "temp_effect":
            return ph["temp_effect"].reindex(strains).values
        if trait == "vir_temp_diff":
            return d_auc.reindex(strains).values
        if trait.startswith("mort_"):
            return ph[trait].reindex(strains).values
        return (1 - ap[f"auc_{trait[len('vir_'):]}"].reindex(strains).values)

    kept = []
    for h in auc[auc.q < 0.1].sort_values("q").HOG_ID:
        v = mat.loc[h].astype(float).values
        if all(abs(np.corrcoef(v, mat.loc[k].astype(float).values)[0, 1]) < 0.8 for k in kept):
            kept.append(h)
    print(f"{'screen':10s} {'HOG':16s} {'trait':14s} {'rho':>7s} {'P':>9s} "
          f"{'rho|PC':>8s} {'P|PC':>9s} {'survives':>9s}")
    for name, d, hits in (("AUC", auc, kept),
                          ("mortality", mort, list(mort.nsmallest(5, "p").HOG_ID))):
        for h in hits:
            row = d[d.HOG_ID == h].nsmallest(1, "p").iloc[0]
            y = yvec(row.trait)
            x = mat.loc[h].astype(float).values
            raw = stats.spearmanr(x, y)
            r, pv = _partial_spearman(x, y, pcs[:, :2])
            ok = "yes" if pv < 0.05 else "no"
            print(f"{name:10s} {h:16s} {row.trait:14s} {raw.statistic:+7.3f} "
                  f"{raw.pvalue:9.2e} {r:+8.3f} {pv:9.4f} {ok:>9s}")
            rows.append(dict(layer="structure", item=f"{name}:{h}~{row.trait}",
                             metric="partial_rho_pc12", value=round(r, 4)))
            rows.append(dict(layer="structure", item=f"{name}:{h}~{row.trait}",
                             metric="partial_p_pc12", value=round(pv, 6)))


def main():
    m, a, vir = pheno_frames()
    rows = []
    layer1(m, a, vir, rows)
    layer2(rows)
    layer3(rows)
    strains = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0).index.tolist()
    ph = pd.read_csv(f"{PA}/phenotype_metrics.tsv", sep="\t", index_col=0)
    ap = pd.read_csv(f"{PA}/phenotype_auc.tsv", sep="\t", index_col=0)
    layer4(rows, load_hog(strains), strains, ph, ap)
    out = f"{GW}/mortality_vs_auc_comparison.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
