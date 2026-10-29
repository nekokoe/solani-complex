#!/usr/bin/env python3
"""Per-strain temperature preference by Cochran-Mantel-Haenszel, stratified by inoculum dose.

WHY THIS EXISTS
    A hard cutoff on the continuous difference in mean mortality (an absolute threshold,
    or mean ± k·s.d.) was tried and rejected: in a unimodal distribution `mean ± 1 s.d.` flags ~1/3 of the
    cohort by construction, and both cutoffs put strains on either side of the line on
    margins of 0.002-0.02. The preference is a comparison of binomial counts (dead vs
    alive, 37 °C vs 28 °C) at three inoculum doses, so the exact test for it is a
    stratified CMH with dose as the stratum. It needs no user-chosen constant and handles
    the unequal number of larvae per cell (n0 = 9-14).

INPUT
    07.lethal_groups/caterpillar_21_narrow.txt
        Sample | Measure (= <temp>C<dose>_dayNN) | Count (alive) | Percent
        day00 gives n0 per cell, day10 the survivors.

OUTPUT
    07.lethal_groups/phenotype_analysis/temperature_bias_cmh.tsv
    07.lethal_groups/phenotype_analysis/temperature_bias.tsv   (class column replaced)

    class: 37C_biased if P < 0.05 and OR > 1; 28C_biased if P < 0.05 and OR < 1; else flat
    q: Benjamini-Hochberg over the 21 strains.
"""

import csv
from collections import defaultdict
from math import erfc, sqrt
from pathlib import Path

ROOT = Path("/share/pasteur/gaosh/Fusarium/phase2/solaniPathogenicityPaper")
RAW = ROOT / "07.lethal_groups/caterpillar_21_narrow.txt"
OUTDIR = ROOT / "07.lethal_groups/phenotype_analysis"
DOSES = ["e3", "e4", "e5"]
TEMPS = ["28C", "37C"]
ALPHA = 0.05


def norm_sf(z):
    """Upper tail of the standard normal."""
    return 0.5 * erfc(z / sqrt(2))


def read_counts():
    """-> {(strain, temp, dose): (n0, survivors_at_day10)}"""
    n0, surv = {}, {}
    for r in csv.DictReader(open(RAW), delimiter="\t"):
        cond, day = r["Measure"].split("_day")
        temp, dose = cond[:3], cond[3:]
        key = (r["Sample"], temp, dose)
        if day == "00":
            n0[key] = int(r["Count"])
        elif day == "10":
            surv[key] = int(r["Count"])
    missing = set(n0) ^ set(surv)
    assert not missing, f"day00/day10 mismatch: {sorted(missing)[:5]}"
    return {k: (n0[k], surv[k]) for k in n0}


def cmh_for_strain(counts, strain):
    """Stratified CMH over the three doses. Table per stratum: [[dead37, alive37], [dead28, alive28]]."""
    num, den, dev, var = 0.0, 0.0, 0.0, 0.0
    strata = []
    for dose in DOSES:
        n37, s37 = counts[(strain, "37C", dose)]
        n28, s28 = counts[(strain, "28C", dose)]
        a, b = n37 - s37, s37          # dead 37, alive 37
        c, d = n28 - s28, s28          # dead 28, alive 28
        n = a + b + c + d
        strata.append((dose, a, b, c, d, n))
        num += a * d / n
        den += b * c / n
        dev += a - (a + b) * (a + c) / n
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
    chi2 = (abs(dev) - 0.5) ** 2 / var      # continuity-corrected: the reported statistic
    p = 2 * norm_sf(sqrt(chi2))             # chi2(1) tail = 2 * normal upper tail
    return num / den, chi2, p, strata


def bh(pvals):
    """Benjamini-Hochberg q values, same order as input."""
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    q = [0.0] * len(pvals)
    m, prev = len(pvals), 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * m / (rank + 1))
        q[i] = prev
    return q


def main():
    counts = read_counts()
    strains = sorted({k[0] for k in counts})
    assert len(strains) == 21, f"expected 21 strains, got {len(strains)}"

    rows = []
    for s in strains:
        or_mh, chi2, p, strata = cmh_for_strain(counts, s)
        rows.append(dict(strain=s, OR=or_mh, chi2=chi2, P=p, strata=strata))

    qs = bh([r["P"] for r in rows])
    for r, q in zip(rows, qs):
        r["q"] = q
        if r["P"] < ALPHA and r["OR"] > 1:
            r["class"] = "37C_biased"
        elif r["P"] < ALPHA and r["OR"] < 1:
            r["class"] = "28C_biased"
        else:
            r["class"] = "flat"

    # continuous effect size, reported for continuity with phenotype_metrics.tsv but NOT used
    # for classification: mean of the three per-dose mortality rates, not a pooled ratio
    eff = {}
    for s in strains:
        m = []
        for temp in TEMPS:
            props = [(counts[(s, temp, d)][0] - counts[(s, temp, d)][1]) / counts[(s, temp, d)][0]
                     for d in DOSES]
            m.append(sum(props) / len(props))
        eff[s] = m[0] - m[1]           # mortality(28) - mortality(37), positive = 28-biased

    out = OUTDIR / "temperature_bias_cmh.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["strain", "temp_effect", "OR_37_vs_28", "chi2", "P", "q_BH",
                    "class", "dead37", "alive37", "dead28", "alive28", "strata"])
        for r in sorted(rows, key=lambda r: r["P"]):
            s = r["strain"]
            d37 = sum(x[1] for x in r["strata"])
            a37 = sum(x[2] for x in r["strata"])
            d28 = sum(x[3] for x in r["strata"])
            a28 = sum(x[4] for x in r["strata"])
            w.writerow([s, f"{eff[s]:.6f}", f"{r['OR']:.4f}", f"{r['chi2']:.3f}",
                        f"{r['P']:.5f}", f"{r['q']:.4f}", r["class"],
                        d37, a37, d28, a28,
                        ";".join(f"{x[0]}:{x[1]}/{x[2]}/{x[3]}/{x[4]}(n={x[5]})" for x in r["strata"])])

    print(f"{'strain':10s} {'delta':>8s} {'OR':>7s} {'chi2':>7s} {'P':>9s} {'q':>7s}  class")
    for r in sorted(rows, key=lambda r: r["P"]):
        print(f"{r['strain']:10s} {eff[r['strain']]:+8.4f} {r['OR']:7.3f} {r['chi2']:7.2f} "
              f"{r['P']:9.5f} {r['q']:7.4f}  {r['class']}")

    cls = {c: [r["strain"] for r in rows if r["class"] == c]
           for c in ("37C_biased", "28C_biased", "flat")}
    print(f"\n37C-biased ({len(cls['37C_biased'])}): {', '.join(sorted(cls['37C_biased']))}")
    print(f"28C-biased ({len(cls['28C_biased'])}): {', '.join(sorted(cls['28C_biased']))}")
    print(f"flat ({len(cls['flat'])}): {', '.join(sorted(cls['flat']))}")
    print(f"P<0.05 = {sum(r['P'] < 0.05 for r in rows)} · P<0.10 = {sum(r['P'] < 0.10 for r in rows)}"
          f" · q<0.05 = {sum(r['q'] < 0.05 for r in rows)} · q<0.10 = {sum(r['q'] < 0.10 for r in rows)}")
    print(f"\nwrote {out}")

    # keep the legacy filename in step so nothing downstream reads a stale class column
    legacy = OUTDIR / "temperature_bias.tsv"
    meta = {}
    if legacy.exists():
        for r in csv.DictReader(open(legacy), delimiter="\t"):
            meta[r["strain"]] = (r["host"], r["med_env"])
    with open(legacy, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["strain", "temp_effect", "host", "med_env", "class"])
        for s in sorted(strains, key=lambda s: eff[s]):
            host, me = meta.get(s, ("", ""))
            klass = next(r["class"] for r in rows if r["strain"] == s)
            w.writerow([s, eff[s], host, me, klass])
    print(f"updated {legacy}")


if __name__ == "__main__":
    main()
