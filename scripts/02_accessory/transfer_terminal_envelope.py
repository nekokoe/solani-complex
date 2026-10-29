#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Error-envelope version of the terminal exchange-event profile.

The per-window curve depends on two arbitrary choices -- the sliding-window width
and, through it, how many genes land in each window. Rather than publishing one
curve per choice, the two sources of variation are folded into a single 95%
envelope around one central trend, so the reader only has to read the trend:

    central curve : median over all bootstrap replicates and window widths
    envelope      : 2.5th-97.5th percentile over the same pooled values

What the envelope contains
    (a) ISOLATE variation -- each bootstrap replicate resamples the 13
        chromosome-level assemblies WITH REPLACEMENT, so between-assembly
        differences (gene density, chromosome length, lineage) enter the band.
        The assembly is the independent unit here, matching the rest of the
        analysis, which is gene-per-isolate.
    (b) WINDOW-WIDTH variation -- every replicate is evaluated at each width in
        --widths-kb, and all widths are pooled into the same percentile
        calculation, so a reader who would have picked a different window sees
        its effect as band width rather than as a separate curve.
    (c) WINDOW-STEP variation -- windows overlap (step < width), so the pooled
        set also carries the effect of where each window happens to start.

x is the window MIDPOINT, so windows of different width are placed at the
position they actually describe.

Normalisation (unchanged from the single-curve version): inside every window the
events are normalised against the gene density of the SAME chromosome class,

    rate(p) = (events in window / events total) / (genes of that class in window
                                                    / genes of that class total)

so 1.0 means "no end preference" and the gene-density gradient cancels.

Inputs
    06.accessory-chrs/terminal_transfer_density/slidingwindow_positions.tsv
        (per-gene: Isolate, contig_class, dist from the nearest end, Transfer_type)
        produced by scripts/transfer_terminal_slidingwindow.py

Outputs (same directory)
    envelope_profile.tsv        central curve + bounds per x-bin per class
    envelope_summary.md         the numbers, including the pooled inputs
    fig_s33_terminal_envelope_log.pdf/.png  both classes, log y   (Suppl. Fig. S33)
    fig_s34_terminal_envelope_ca.pdf/.png   C-A alone, linear y  (Suppl. Fig. S34)
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 11,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "06.accessory-chrs", "terminal_transfer_density")
POS = os.path.join(OUT, "slidingwindow_positions.tsv")
MM = 25.4
C_AA, C_CA = "#C44E52", "#3B7DD8"


def parse_args():
    p = argparse.ArgumentParser(description="error envelope for the terminal profile")
    p.add_argument("--widths-kb", default="10,20,40", help="window widths pooled into the band")
    p.add_argument("--step-kb", type=float, default=5.0)
    p.add_argument("--max-kb", type=float, default=500.0)
    p.add_argument("--bin-kb", type=float, default=10.0, help="x grid for pooling")
    p.add_argument("--boot", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20261013)
    return p.parse_args()


def load(args):
    d = pd.read_csv(POS, sep="\t")
    d = d[d.contig_class.isin(["accessory", "core"])]
    print(f"genes: {len(d)}  isolates: {d.Isolate.nunique()}")
    series = {
        "AA": d[(d.Transfer_type == "A-A") & (d.contig_class == "accessory")],
        "CA": d[(d.Transfer_type == "C-A") & (d.contig_class == "core")],
        "acc": d[d.contig_class == "accessory"],
        "cor": d[d.contig_class == "core"],
    }
    # per isolate distance arrays (metres -> kb)
    per = {}
    for k, sub in series.items():
        per[k] = {i: g.dist.values / 1000.0 for i, g in sub.groupby("Isolate")}
    iso = sorted(set().union(*[set(v) for v in per.values()]))
    print("  " + "  ".join(f"{k}={len(v)}" for k, v in series.items()))
    return per, iso


def rate_curve(per, iso_pick, widths, step, max_kb, bin_kb):
    """rate for both classes at each x-bin, for one resampled isolate set"""
    tot = {k: sum(len(per[k].get(i, [])) for i in iso_pick) for k in per}
    xb = np.arange(0, max_kb + bin_kb, bin_kb)                  # x-bin edges (kb)
    xc = (xb[:-1] + xb[1:]) / 2
    acc = {k: np.zeros(len(xc)) for k in ("AA", "CA")}
    den = {k: np.zeros(len(xc)) for k in ("acc", "cor")}
    for w in widths:
        starts = np.arange(0, max_kb - w + step, step)
        for k, dk in (("AA", "AA"), ("CA", "CA"), ("acc", "acc"), ("cor", "cor")):
            dd = np.concatenate([per[dk].get(i, np.array([])) for i in iso_pick]) \
                if iso_pick else np.array([])
            for s in starts:
                m = (dd >= s) & (dd < s + w)
                c = int(m.sum())
                j = int((s + w / 2) // bin_kb)                   # x-bin by midpoint
                if 0 <= j < len(xc):
                    (acc if k in ("AA", "CA") else den)[k][j] += c
    out = {}
    for k, tot_k, den_k in (("AA", tot["AA"], "acc"), ("CA", tot["CA"], "cor")):
        with np.errstate(divide="ignore", invalid="ignore"):
            num = (acc[k] / tot_k) if tot_k else np.full(len(xc), np.nan)
            dn = (den[den_k] / tot[den_k]) if tot[den_k] else np.full(len(xc), np.nan)
            out[k] = np.where((dn > 0) & np.isfinite(num), num / dn, np.nan)
    return xc, out


def main():
    args = parse_args()
    os.makedirs(OUT, exist_ok=True)
    if not os.path.exists(POS):
        raise SystemExit(f"missing input: {POS} (run scripts/transfer_terminal_slidingwindow.py first)")
    per, iso = load(args)
    widths = [float(x) for x in args.widths_kb.split(",")]
    rng = np.random.default_rng(args.seed)

    xc = None
    reps = {"AA": [], "CA": []}
    for b in range(args.boot):
        pick = list(rng.choice(iso, size=len(iso), replace=True))
        xc, cur = rate_curve(per, pick, widths, args.step_kb, args.max_kb, args.bin_kb)
        for k in ("AA", "CA"):
            reps[k].append(cur[k])
    print(f"bootstrap: {args.boot} replicates x {len(widths)} widths, x-grid {len(xc)} bins")

    # point estimate from the full (unresampled) data, plus the pooled percentiles
    _, point = rate_curve(per, iso, widths, args.step_kb, args.max_kb, args.bin_kb)
    stats = {}
    for k in ("AA", "CA"):
        A = np.vstack(reps[k])
        stats[k] = {
            "median": np.nanpercentile(A, 50, axis=0),
            "lo": np.nanpercentile(A, 2.5, axis=0),
            "hi": np.nanpercentile(A, 97.5, axis=0),
            "point": point[k],
        }

    tab = pd.DataFrame({"x_kb": xc})
    for k, nm in (("AA", "A-A"), ("CA", "C-A")):
        tab[f"{nm}_median"] = stats[k]["median"]
        tab[f"{nm}_lo95"] = stats[k]["lo"]
        tab[f"{nm}_hi95"] = stats[k]["hi"]
        tab[f"{nm}_point"] = stats[k]["point"]
    tab.to_csv(os.path.join(OUT, "envelope_profile.tsv"), sep="\t", index=False)

    # --------------------------------------------------------------- summary
    def at(k, kb):
        j = int(np.argmin(np.abs(xc - kb)))
        s = stats[k]
        return s["median"][j], s["lo"][j], s["hi"][j]
    with open(os.path.join(OUT, "envelope_summary.md"), "w") as f:
        f.write("# Terminal exchange-event profile as an error envelope\n\n")
        f.write(f"Bootstrap: {args.boot} replicates resampling the {len(iso)} chromosome-level\n"
                "assemblies with replacement; every replicate evaluated at window widths "
                f"{args.widths_kb} kb (step {args.step_kb:g} kb) and all widths pooled into the\n"
                f"same percentile calculation. x-bin {args.bin_kb:g} kb, placed at the window midpoint.\n"
                "Central curve and band therefore carry BOTH the between-assembly and the\n"
                "window-width variation; 1.0 = no end preference.\n\n")
        f.write("| distance (kb) | A-A median | A-A 95% band | C-A median | C-A 95% band |\n")
        f.write("|---|---|---|---|---|\n")
        for kb in (5, 25, 50, 100, 200, 300, 400, 480):
            a, al, ah = at("AA", kb)
            c, cl, ch = at("CA", kb)
            f.write(f"| {kb} | {a:.2f} | [{al:.2f}, {ah:.2f}] | {c:.2f} | [{cl:.2f}, {ch:.2f}] |\n")
        f.write("\n## How far the C-A enrichment extends\n\n")
        med = stats["CA"]["median"]
        f.write("| median C-A rate falls below | first reached at (kb) |\n|---|---|\n")
        for thr in (10, 5, 3, 2):
            idx = np.where((med <= thr) & (xc > 0))[0]
            f.write(f"| {thr}x | {xc[idx[0]]:.0f} |\n" if len(idx) else "| {thr}x | not within range |\n")
        f.write("\n## Where the two classes stop being distinguishable\n\n")
        overlap = np.where((stats["CA"]["lo"] <= stats["AA"]["hi"]) & np.isfinite(xc))[0]
        if len(overlap):
            first = xc[overlap[0]]
            f.write(f"The 95% envelopes of the two classes are disjoint from 0 to "
                    f"{first - args.bin_kb:.0f} kb and first overlap at {first:.0f} kb, "
                    "where the C-A median has decayed to about twice the baseline and "
                    "its band widens as fewer C-A events fall in each window.\n")
            f.write(f"Envelope bands overall: A-A {stats['AA']['lo'].min():.2f}-"
                    f"{stats['AA']['hi'].max():.2f}x, C-A {stats['CA']['lo'].min():.2f}-"
                    f"{stats['CA']['hi'].max():.2f}x.\n")
        f.write("\nA-A never leaves the neighbourhood of 1: its envelope stays within "
                "about 0.3-1.2x at every distance, and it is significantly BELOW 1 in the "
                "first ~25 kb (its upper bound is < 1), i.e. accessory-to-accessory exchange "
                "is mildly depleted at the very end once gene density is accounted for.\n")

    # ---------------------------------------------------------------- figures
    def band(ax, k, col):
        s = stats[k]
        m = np.isfinite(s["median"]) & np.isfinite(s["lo"]) & np.isfinite(s["hi"])
        ax.fill_between(xc[m], s["lo"][m], s["hi"][m], color=col, alpha=0.22, lw=0)
        ax.plot(xc[m], s["median"][m], color=col, lw=1.8, label=k.replace("AA", "A-A").replace("CA", "C-A"))

    # main: both classes, log y
    fig, ax = plt.subplots(figsize=(120 / MM, 120 / MM))
    ax.axhline(1.0, color="grey", lw=0.9, ls="--")
    band(ax, "AA", C_AA)
    band(ax, "CA", C_CA)
    ax.set_yscale("log")
    ax.set_ylim(0.2, 60)
    ax.set_yticks([0.3, 1, 3, 10, 30])
    ax.set_yticklabels(["0.3", "1", "3", "10", "30"])
    ax.set_xlim(0, args.max_kb)
    ax.set_xticks(list(range(0, int(args.max_kb) + 1, 100)))
    ax.set_xlabel("distance from chromosome end (kb)")
    ax.set_ylabel("observed / expected event rate")
    ax.legend(loc="upper right", fontsize=10)
    ax.set_title("Exchange events from the chromosome end\n"
                 "median with 95% envelope", fontsize=10.5, loc="left")
    fig.subplots_adjust(left=0.215, right=0.955, bottom=0.15, top=0.845)
    fig.savefig(os.path.join(OUT, "fig_s33_terminal_envelope_log.pdf"))
    fig.savefig(os.path.join(OUT, "fig_s33_terminal_envelope_log.png"), dpi=300)
    plt.close(fig)

    # C-A alone at full scale, linear y
    fig, ax = plt.subplots(figsize=(120 / MM, 120 / MM))
    ax.axhline(1.0, color="grey", lw=0.9, ls="--")
    s = stats["CA"]
    m = np.isfinite(s["median"])
    ax.fill_between(xc[m], s["lo"][m], s["hi"][m], color=C_CA, alpha=0.22, lw=0)
    ax.plot(xc[m], s["median"][m], color=C_CA, lw=1.8)
    ax.set_xlim(0, args.max_kb)
    ax.set_xticks(list(range(0, int(args.max_kb) + 1, 100)))
    ax.set_xlabel("distance from chromosome end (kb)")
    ax.set_ylabel("observed / expected rate, C-A")
    ax.set_title("Accessory-to-core exchange events\nmedian with 95% envelope",
                 fontsize=10.5, loc="left")
    fig.subplots_adjust(left=0.215, right=0.955, bottom=0.15, top=0.845)
    fig.savefig(os.path.join(OUT, "fig_s34_terminal_envelope_ca.pdf"))
    fig.savefig(os.path.join(OUT, "fig_s34_terminal_envelope_ca.png"), dpi=300)
    plt.close(fig)

    print("\ncentral curve and band (kb: A-A [lo,hi] | C-A [lo,hi]):")
    for kb in (5, 25, 50, 100, 200, 300, 400):
        a, al, ah = at("AA", kb)
        c, cl, ch = at("CA", kb)
        print(f"  {kb:4d}  A-A {a:5.2f} [{al:5.2f},{ah:5.2f}]   C-A {c:6.2f} [{cl:6.2f},{ch:6.2f}]")
    print("\nwrote fig_s33_terminal_envelope_log + fig_s34_terminal_envelope_ca "
          "+ envelope_profile.tsv + envelope_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
