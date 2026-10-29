#!/usr/bin/env python3
"""Supplementary Table S16 - outbreak-window lineages and their persistence.

Builds the workbook behind the final paragraph of the population section:

  summary            per species: cohort size, outbreak-window counts, Mantel test of
                     genetic distance against the interval between collection dates,
                     same-year versus different-year distances, and the year enrichment
                     among the closest 1% of clinical pairs
  nearest_neighbour   per clinical isolate: its closest relative, that relative's
                     collection year, the gap and the distance
  subpopulation_spans per ADMIXTURE component: size, how many were collected inside
                     and after the outbreak window, and the span of collection years
  definitions        the metric, the window, the permutation scheme and the software

Distance = simple-matching dissimilarity over non-missing core-SNP genotypes
(1 - matching proportion), computed from tmp/amova/<species>_dist_matrix.tsv. On this
scale an unrelated pair sits near 0.5 and near-identical isolates near 0; it is a
relative measure and is not a count of differing sites.

Only BMU isolates enter the year tests, because collection year exists only for them
(Supplementary Table S1); public reference genomes carry no collection year. The
outbreak window is 2000-2005 inclusive. Subpopulation labels follow the ADMIXTURE
column order for K = 3 (F. solani, F. falciformis) and K = 2 (F. keratoplastica) and
are therefore reported as "component N", not as the Pop1/Pop2/Pop3 labels of Fig. 2E.

  python3 04.population/outbreak_coexistence_tables.py
"""

import glob
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
META = os.path.join(ROOT, "00.manuscript", "supplementary", "suppl_tables.genome_part.0224-GS.xlsx")
DIST = os.path.join(ROOT, "tmp", "amova")
OUT = os.path.join(ROOT, "00.manuscript", "supplementary",
                   "suppl_table.S16_outbreak_lineages.xlsx")
SPECIES = [("falciformis", 3), ("solani", 3), ("keratoplastica", 2)]
WIN = (2000, 2005)
NPERM = 10000
SEED = 42


def load_meta():
    m = pd.read_excel(META, sheet_name="Table 1", header=1)
    m = m.rename(columns={"Isolate": "isolate", "Isolation Source": "source",
                          "Identified Species": "species", "Collection Year": "year"})
    m = m[["isolate", "source", "species", "year"]].dropna()
    m["year"] = m["year"].astype(int)
    return m


def load_components(sp, k):
    q = pd.read_csv(glob.glob(os.path.join(ROOT, "04.population", f"solani-short-{sp}",
                                           f"*maf5pct.recode.{k}.Q"))[0], sep=r"\s+", header=None)
    fam = glob.glob(os.path.join(ROOT, "04.population", f"solani-short-{sp}", "*maf5pct.recode.fam"))[0]
    fid = [l.split()[0] for l in open(fam)]
    iid = [l.split()[1] for l in open(fam)]
    names = iid if len(set(iid)) == len(iid) else fid
    names = [n.split(".")[0] for n in names]
    assert len(names) == len(q)
    return pd.Series(q.idxmax(axis=1).map(lambda c: f"component_{c+1}").values, index=names)


def main():
    meta = load_meta()
    rng = np.random.default_rng(SEED)
    summary, nbr, spans = [], [], []

    for sp, k in SPECIES:
        dm = pd.read_csv(os.path.join(DIST, f"{sp}_dist_matrix.tsv"), sep="\t", index_col=0)
        clin = meta[(meta.species == sp) & (meta.source == "med") & (meta.isolate.isin(dm.index))]
        ids, years = clin.isolate.tolist(), clin.year.to_numpy()
        sub = dm.loc[ids, ids].to_numpy()
        n = len(ids)
        iu = np.triu_indices(n, 1)
        d, dy = sub[iu], np.abs(years[:, None] - years[None, :])[iu]
        same = dy == 0

        # Mantel: distance vs |year difference|
        y = years.astype(float)
        obs = np.corrcoef(d, dy)[0, 1]
        null = np.empty(NPERM)
        for i in range(NPERM):
            p = rng.permutation(n)
            null[i] = np.corrcoef(d, np.abs(y[p][:, None] - y[p][None, :])[iu])[0, 1]
        mantel_p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (NPERM + 1)

        # same-year versus different-year
        obs_gap = d[~same].mean() - d[same].mean()
        null_gap = np.empty(NPERM)
        for i in range(NPERM):
            p = rng.permutation(years)
            m = (np.abs(p[:, None] - p[None, :])[iu]) == 0
            null_gap[i] = d[~m].mean() - d[m].mean()
        gap_p = (np.sum(np.abs(null_gap) >= abs(obs_gap)) + 1) / (NPERM + 1)

        # year enrichment among the closest 1%
        cut = np.quantile(d, 0.01)
        near = np.where(d <= cut)[0]
        ii, jj = iu[0][near], iu[1][near]
        same_near = int(np.sum(years[ii] == years[jj]))
        exp_near = len(near) * same.mean()
        null_same = np.empty(NPERM, dtype=int)
        for i in range(NPERM):
            p = rng.permutation(years)
            null_same[i] = int(np.sum(p[ii] == p[jj]))
        near_p = (np.sum(null_same >= same_near) + 1) / (NPERM + 1)

        inwin = (years >= WIN[0]) & (years <= WIN[1])
        comp = load_components(sp, k)
        comp = comp.reindex(ids)

        summary.append(dict(
            species=sp, n_clinical=n, year_min=int(years.min()), year_max=int(years.max()),
            n_in_window=int(inwin.sum()), n_after_window=int((years > WIN[1]).sum()),
            n_pairs=len(d), median_pairwise_distance=round(float(np.median(d)), 4),
            mantel_r=round(float(obs), 4), mantel_p=round(float(mantel_p), 4),
            mean_distance_same_year=round(float(d[same].mean()), 4),
            mean_distance_diff_year=round(float(d[~same].mean()), 4),
            same_vs_diff_perm_p=round(float(gap_p), 4),
            closest_pairs_cutoff=round(float(cut), 4), n_closest_pairs=len(near),
            n_closest_same_year=same_near, expected_same_year=round(float(exp_near), 2),
            closest_same_year_perm_p=round(float(near_p), 4),
            n_permutations=NPERM, seed=SEED,
        ))

        sub_nn = sub.copy()
        np.fill_diagonal(sub_nn, np.inf)
        nn = sub_nn.argmin(axis=1)
        for i, iso in enumerate(ids):
            nbr.append(dict(
                species=sp, isolate=iso, year=int(years[i]),
                component=comp.iloc[i],
                nearest_isolate=ids[nn[i]], nearest_year=int(years[nn[i]]),
                year_gap=int(abs(years[i] - years[nn[i]])),
                distance=round(float(sub_nn[i, nn[i]]), 4),
                in_outbreak_window=bool(inwin[i]),
                nearest_in_window=bool(WIN[0] <= years[nn[i]] <= WIN[1]),
            ))

        for c, g in pd.DataFrame({"iso": ids, "year": years, "component": comp.values}).groupby("component"):
            spans.append(dict(
                species=sp, component=c, n=len(g),
                n_in_window=int(g.year.between(*WIN).sum()),
                n_after_window=int((g.year > WIN[1]).sum()),
                year_min=int(g.year.min()), year_max=int(g.year.max()),
                year_span=int(g.year.max() - g.year.min()),
                contains_post_window_isolate=bool((g.year > WIN[1]).any()),
            ))

    definitions = pd.DataFrame({
        "item": ["distance metric", "metric scale", "outbreak window", "cohort for year tests",
                 "permutation scheme", "closest-pair set", "subpopulation labels",
                 "coexistence claim", "software"],
        "value": [
            "simple-matching dissimilarity over non-missing core-SNP genotypes (1 - matching proportion), "
            "from tmp/amova/<species>_dist_matrix.tsv",
            "relative measure: unrelated pairs sit near 0.5, near-identical isolates near 0; not a count of differing sites",
            f"{WIN[0]}-{WIN[1]} inclusive",
            "BMU isolates with a known collection year (Supplementary Table S1); public reference genomes carry none",
            f"{NPERM} label permutations, seed {SEED}; P = (1 + number of null values as extreme) / (1 + permutations)",
            "the 1% of clinical pairs with the smallest distance, per species",
            "ADMIXTURE column order for the K reported in the main text (K = 3 for F. solani and F. falciformis, "
            "K = 2 for F. keratoplastica); labelled component_N so as not to conflict with the Pop1/Pop2/Pop3 labels of Fig. 2E",
            "lineages are called persistent when they contain isolates collected on both sides of the window",
            "python 3.13, numpy, pandas, scipy; script 04.population/outbreak_coexistence_tables.py",
        ],
    })

    with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
        pd.DataFrame(summary).to_excel(xw, sheet_name="summary", index=False)
        pd.DataFrame(nbr).to_excel(xw, sheet_name="nearest_neighbour", index=False)
        pd.DataFrame(spans).to_excel(xw, sheet_name="subpopulation_spans", index=False)
        definitions.to_excel(xw, sheet_name="definitions", index=False)

    print(pd.DataFrame(summary).to_string(index=False))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
