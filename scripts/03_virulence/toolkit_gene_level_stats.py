#!/usr/bin/env python3
"""Every gene- and HOG-level statistic behind the toolkit paragraph (p154).

The paragraph contrasts the two strain-specific toolkits at two resolutions: the HOG-level
annotation that colours Supplementary Fig. S27, and the annotation of the genes the biased
strains actually carry. Each number is recomputed here from the current toolkit HOG sets,
so the text can be rewritten from a script instead of from ad-hoc arithmetic.

Sets come from `build_toolkit_hog_sets.py` (CMH grouping, sp37=15, sp28=18).

Reports
  1. set sizes and HOG-level CAZy composition, with the two-sided Fisher test
  2. chitin-related and infection-relevant HOGs in each set
  3. gene-level annotation of the biased strains' own copies: CAZy hits, Pfam coverage
  4. HOG-level Pfam enrichment per set against the 4,484-variable-HOG background (BH q)
  5. Pfam domains shared by the two sets
"""

import os
from collections import Counter

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GW = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc")
PER = os.path.join(ROOT, "02.annotation/per_isolates")
HOGCOUNT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv")
ANNOT = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
GVH = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv")
BIAS = os.path.join(ROOT, "07.lethal_groups/phenotype_analysis/temperature_bias_cmh.tsv")

CHIT = ["Glyco_hydro_18", "Glyco_hydro_19", "Glyco_hydro_20", "Glyco_hydro_12",
        "Chitin_synth", "chitin", "Glyco_hydro_28"]
CARRY_MIN, CARRY_MAX = 3, 19


def fisher_two_sided(a, b, c, d):
    from math import comb
    n = a + b + c + d
    row1, col1 = a + b, a + c
    pmf = lambda x: comb(row1, x) * comb(n - row1, col1 - x) / comb(n, col1)
    p_obs = pmf(a)
    lo, hi = max(0, col1 - (n - row1)), min(row1, col1)
    return min(1.0, sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= p_obs * (1 + 1e-9)))


def bh(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    q = [0.0] * len(pvals)
    m, prev = len(pvals), 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * m / (rank + 1))
        q[i] = prev
    return q


def hypergeom_sf(k, m, n, K):
    """P(X >= k) with X ~ Hypergeom(population m+n, successes m, draws K)."""
    from math import comb
    K = min(K, m + n)
    if k <= 0:
        return 1.0
    return sum(comb(m, x) * comb(n, K - x) for x in range(k, min(m, K) + 1)) / comb(m + n, K)


def load_gene_annot(iso):
    pf = pd.read_csv(f"{PER}/{iso}/{iso}.gene_pfam.tsv", sep="\t", header=None,
                     usecols=[0, 2], names=["g", "dom"])
    pfam = pf.groupby("g")["dom"].apply(lambda s: sorted(set(s))).to_dict()
    cz = pd.read_csv(f"{PER}/{iso}/{iso}.gene_dbcan.tsv", sep="\t", header=None,
                     names=["g", "v"])
    cazy = cz.groupby("g")["v"].apply(lambda s: sorted(set(s))).to_dict()
    return pfam, cazy


def main():
    pheno = pd.read_csv(os.path.join(ROOT, "07.lethal_groups/phenotype_analysis/"
                                           "phenotype_metrics.tsv"), sep="\t", index_col=0)
    strains = list(pheno.index)
    bias = pd.read_csv(BIAS, sep="\t").set_index("strain")
    B37 = sorted(bias.index[bias["class"] == "37C_biased"])
    B28 = sorted(bias.index[bias["class"] == "28C_biased"])

    sets = {}
    for name in ("sp37", "sp28"):
        df = pd.read_csv(os.path.join(GW, f"hog_{name}.tsv"), sep="\t", usecols=[0])
        sets[name] = list(df.iloc[:, 0])
    sp37, sp28 = sets["sp37"], sets["sp28"]

    annot = pd.read_csv(ANNOT, sep="\t", dtype=str).drop_duplicates("HOG_ID").set_index("HOG_ID")
    def doms(h):
        v = annot.at[h, "Pfam_Domain"] if h in annot.index else ""
        return [d for d in str(v).split(";") if d and d.lower() not in ("nan", "-")]
    def cazy_class(h):
        v = annot.at[h, "CAZy_Class"] if h in annot.index else ""
        return "" if str(v).lower() in ("nan", "-", "") else str(v)
    def gene_name(h):
        v = annot.at[h, "Gene_Name"] if h in annot.index else ""
        return "" if str(v).lower() in ("nan", "-", "") else str(v)

    print("=" * 78)
    print("1. set sizes and HOG-level CAZy composition")
    print("=" * 78)
    counts = {}
    for name, ids, group in (("sp37", sp37, B37), ("sp28", sp28, B28)):
        k = sum(1 for h in ids if cazy_class(h))
        counts[name] = (k, len(ids) - k)
        print(f"  {name}: {len(ids)} HOGs · {k} carry a CAZy class ({k / len(ids):.1%}) "
              f"· grouped strains {', '.join(group)}")
    p = fisher_two_sided(counts['sp37'][0], counts['sp37'][1],
                         counts['sp28'][0], counts['sp28'][1])
    print(f"  CAZy fraction sp37 vs sp28: two-sided Fisher P = {p:.4f}")

    print("\n" + "=" * 78)
    print("2. chitin-related / infection-relevant HOGs per set")
    print("=" * 78)
    INF = ["Hce2", "CFEM", "LysM", "HET", "Peptidase_S", "Effector"]
    for name, ids in (("sp37", sp37), ("sp28", sp28)):
        chit = [h for h in ids if any(c in ";".join(doms(h)) for c in CHIT)]
        inf = [h for h in ids if any(c in ";".join(doms(h)) for c in INF)]
        print(f"  {name}: chitin-related {len(chit)} · infection-relevant {len(inf)}")
        for h in chit:
            hit = [c for c in CHIT if c in ";".join(doms(h))]
            print(f"     {h}  {gene_name(h).split(';')[0]:12s} {','.join(hit)}")

    print("\n" + "=" * 78)
    print("3. gene level: the biased strains' own copies")
    print("=" * 78)
    gvh = pd.read_csv(GVH, sep="\t", usecols=["GeneID", "SampleID", "HOG_ID"])
    maps = {s: load_gene_annot(s) for s in sorted(set(B37 + B28))}
    gene_facts = {}
    for name, ids, group in (("sp37", sp37, B37), ("sp28", sp28, B28)):
        sub = gvh[gvh.HOG_ID.isin(ids) & gvh.SampleID.isin(group)]
        n = len(sub)
        n_cazy = sum(1 for r in sub.itertuples() if maps[r.SampleID][1].get(r.GeneID))
        n_pfam = sum(1 for r in sub.itertuples() if maps[r.SampleID][0].get(r.GeneID))
        dom_count = Counter(d for r in sub.itertuples()
                            for d in maps[r.SampleID][0].get(r.GeneID, []))
        gene_facts[name] = dict(n=n, n_cazy=n_cazy, n_pfam=n_pfam, doms=dom_count)
        print(f"  {name}: {n} genes in {len(group)} strains · dbCAN hit in {n_cazy} · "
              f"Pfam record in {n_pfam} ({n_pfam / n:.0%}) · no Pfam in {n - n_pfam}")
        print(f"     recurring Pfam: {dom_count.most_common(6)}")

    print("\n" + "=" * 78)
    print("4. HOG-level Pfam enrichment vs the 4,484-variable-HOG background")
    print("=" * 78)
    hog = pd.read_csv(HOGCOUNT, sep="\t")
    mat = hog[hog.iloc[:, 0].isin(strains)].pivot(index="HOG_ID", columns="SampleID",
                                                  values="GeneCount").reindex(columns=strains).fillna(0)
    n_pres = (mat > 0).sum(axis=1)
    var = mat[strains].nunique(axis=1) > 1
    bg = list(mat[(n_pres >= CARRY_MIN) & (n_pres <= CARRY_MAX) & var].index)
    print(f"  background: {len(bg)} HOGs")
    bg_counts = Counter(d for h in bg for d in set(doms(h)))
    N = len(bg)
    for name, ids in (("sp37", sp37), ("sp28", sp28)):
        set_counts = Counter(d for h in ids for d in set(doms(h)))
        rows = []
        for d, k in set_counts.items():
            m = bg_counts.get(d, 0)
            K = sum(1 for h in ids if d in set(doms(h)))
            pv = hypergeom_sf(k, m, N - m, len(ids))
            rows.append((d, K, m, pv))
        qs = bh([r[3] for r in rows])
        rows = sorted([(d, K, m, pv, q) for (d, K, m, pv), q in zip(rows, qs)],
                      key=lambda r: r[3])
        print(f"  {name}: {len(rows)} Pfam domains tested with K>=1 · "
              f"best q = {min(r[4] for r in rows):.3f}")
        for d, K, m, pv, q in rows[:5]:
            print(f"     {d:20s} {K:2d}/{len(ids)} in set · {m:4d}/{N} background · "
                  f"P = {pv:.4f} · q = {q:.3f}")

    print("\n" + "=" * 78)
    print("5. Pfam domains shared by the two sets")
    print("=" * 78)
    d37 = set(d for h in sp37 for d in doms(h))
    d28 = set(d for h in sp28 for d in doms(h))
    print(f"  sp37 domains {len(d37)} · sp28 domains {len(d28)} · shared {len(d37 & d28)}")
    print(f"  shared: {', '.join(sorted(d37 & d28))}")

    out = os.path.join(GW, "toolkit_gene_level_stats.tsv")
    with open(out, "w") as fh:
        fh.write("set\tn_hogs\tn_cazy_hog\tgene_copies\tgene_cazy\tgene_pfam\tgene_no_pfam\n")
        for name in ("sp37", "sp28"):
            g = gene_facts[name]
            fh.write(f"{name}\t{len(sets[name])}\t{counts[name][0]}\t{g['n']}\t"
                     f"{g['n_cazy']}\t{g['n_pfam']}\t{g['n'] - g['n_pfam']}\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
