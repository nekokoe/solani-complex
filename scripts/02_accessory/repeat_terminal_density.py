#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Repeat-TYPE composition as a function of distance from the chromosome end,
separately for core (CC) and accessory (AC) chromosomes, with bootstrap
envelopes -- the repeat analogue of the exchange-event terminal profile.

Question it answers: the earlier repeat work compared the TOTAL repeat content of
the terminal region with the accessory chromosomes. The hypothesis here is that
individual repeat TYPES, not repeats as a whole, are what sit at the ends, and
that CC and AC differ in which types those are.

Measurement
    RepeatMasker hits are pooled by type, nested hits dropped so nothing is
    double counted. A window of width W is slid from the chromosome end inward,
    always staying inside the half of the chromosome nearest that end (so the two
    ends never overlap). Inside the window,

        content(d)    = bp of that type in the window / window length
        enrichment(d) = content(d) / that type's content over the whole genome
                        of the SAME chromosome class

    so 1.0 means "this type behaves here as it does everywhere else". Dividing by
    the class-specific genome-wide content cancels the large abundance differences
    between types (DNA is 16% of AC but 2.4% of CC), which makes CC and AC
    readable on one axis and directly comparable.

Uncertainty
    95% envelope from bootstrapping the 13 chromosome-level assemblies WITH
    replacement -- the assembly is the independent unit, as elsewhere in this
    analysis. The window width is a free parameter, so every width in --widths-kb
    is pushed into the same percentile calculation: a reader who would have picked
    another width sees its effect as band width, not as a separate curve.

    Windows are binned by their START distance, not their midpoint. With
    midpoint binning only the narrowest window can populate the first bin (a
    40 kb window's midpoint is 20 kb away), which would hide the width
    uncertainty exactly where it matters most.

Implementation note
    Per-window coverage is precomputed ONCE per (isolate, chromosome, type, width)
    from a cumulative-coverage array at --res-bp resolution, then summed per
    isolate. A bootstrap replicate is therefore just ~780 array additions rather
    than a rescan of 350 chromosomes x windows, which is what makes 400
    replicates affordable.

Inputs (in-repo)
    06.accessory-chrs/tgs_assemblies/<ISO>.assembly.fasta.out         RepeatMasker
    02.annotation/per_isolates/<ISO>/repeatmasker/<ISO>.assembly.fasta.masked.length
    06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id

Outputs (06.accessory-chrs/repeat_terminal_density/)
    repeat_density_profile.tsv    per x-bin, per type, per class: median + 95% band
    repeat_density_abundance.tsv  genome-wide % of each type in each class
    repeat_density_summary.md     the numbers
    fig_s32_repeat_density_byfamily.pdf/.png  per-type distance profiles (Suppl. Fig. S32)
    repeat_density_overall.tsv    overall repeat content: absolute % and enrichment
    fig_s31_repeat_density_overall.pdf/.png   overall repeat content, sliding window (Suppl. Fig. S31)
"""

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "06.accessory-chrs", "repeat_terminal_density")
TGS = os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies")

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]
C_CC, C_AC = "#4A6FA5", "#C44E52"

GROUP = {"Unknown": "Unknown", "Simple_repeat": "Simple_repeat",
         "Low_complexity": "Low_complexity"}
for _f, _g in {
    "LTR/Gypsy": "LTR", "LTR/Copia": "LTR", "LTR/Pao": "LTR", "LTR/Ngaro": "LTR",
    "DNA/TcMar-Fot1": "DNA", "DNA/TcMar-Tc1": "DNA", "DNA/TcMar-Tc2": "DNA",
    "DNA/TcMar-Ant1": "DNA", "DNA/MULE-MuDR": "DNA", "DNA/hAT-Restless": "DNA",
    "DNA/hAT-Ac": "DNA", "DNA/hAT-Tag1": "DNA", "DNA/PiggyBac": "DNA",
    "DNA/Kolobok-H": "DNA", "DNA/PIF-Harbinger": "DNA", "DNA/IS3EU": "DNA",
    "DNA/Academ-H": "DNA", "RC/Helitron": "RC", "RC/Helitron-2": "RC",
    "LINE/Tad1": "LINE", "LINE/CRE": "LINE"}.items():
    GROUP[_f] = _g


def parse_args():
    p = argparse.ArgumentParser(description="repeat-type terminal profile")
    p.add_argument("--widths-kb", default="10,20,40")
    p.add_argument("--pool-with", default=None,
                   help="type used for the per-type panels (default: the 6 most abundant)")
    p.add_argument("--n-panels", type=int, default=6)
    p.add_argument("--step-kb", type=float, default=10.0)
    p.add_argument("--max-kb", type=float, default=500.0)
    p.add_argument("--bin-kb", type=float, default=10.0)
    p.add_argument("--res-bp", type=int, default=200,
                   help="resolution of the cumulative-coverage array")
    p.add_argument("--boot", type=int, default=400)
    p.add_argument("--seed", type=int, default=20261013)
    p.add_argument("--core-min-length", type=int, default=500000)
    p.add_argument("--min-pct", type=float, default=0.30,
                   help="keep a type only if it reaches this %% of either class")
    p.add_argument("--summary-kb", type=float, default=100.0,
                   help="terminal window used for the one-panel summary")
    return p.parse_args()


def cum_cov(intervals, L, res):
    """cumulative bp of merged intervals, sampled every `res` bp (index i = [0, i*res))"""
    n = L // res + 2
    d = np.zeros(L + 2, dtype=np.int32)
    for s, e in intervals:
        d[max(1, s)] += 1
        d[min(L, e) + 1] -= 1
    cov = np.cumsum(d)[:L + 1] > 0
    return np.concatenate([[0], np.cumsum(cov)[::res][:n]]).astype(np.int32)


def load_and_precompute(args):
    """Precompute, per isolate, the per-window window LENGTH (exposure) and the
    per-window bp of every repeat type.

    Exposure is accumulated for EVERY chromosome of a class, whether or not it
    carries the type being measured -- otherwise chromosomes that lack a type
    would drop out of the denominator and every enrichment would be inflated.
    """
    widths = [float(x) for x in args.widths_kb.split(",")]
    starts = {w: np.arange(0, args.max_kb - w + args.step_kb, args.step_kb) for w in widths}
    nw_of = {w: len(starts[w]) for w in widths}      # arange length depends on w
    EXP = {}      # (iso, class, w)          -> window length in bp per window
    BP = {}       # (iso, class, family, w)  -> bp of that family per window
    glen = defaultdict(lambda: defaultdict(int))
    gfam = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    for iso in ISOLATES:
        acc = set(open(os.path.join(TGS, f"{iso}.accessory_chrs.id")).read().split())
        L = {}
        with open(os.path.join(REPO, "02.annotation", "per_isolates", iso,
                               "repeatmasker", f"{iso}.assembly.fasta.masked.length")) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    L[f[0]] = int(f[1])
        keep = {c: ("accessory" if c in acc else
                    "core" if n >= args.core_min_length else None) for c, n in L.items()}
        fam_iv = defaultdict(lambda: defaultdict(list))
        with open(os.path.join(TGS, f"{iso}.assembly.fasta.out")) as fh:
            for i, line in enumerate(fh):
                if i < 3:
                    continue
                p = line.split()
                if len(p) < 15:
                    continue
                c = p[4]
                if keep.get(c) is None:
                    continue
                if len(p) > 15 and p[15] == "*":          # nested -> double counting
                    continue
                f = p[10]
                g = GROUP.get(f) or ("DNA" if f.startswith("DNA/") else
                                     "LTR" if f.startswith("LTR/") else
                                     "LINE" if f.startswith("LINE/") else
                                     "RC" if f.startswith("RC/") else "Unknown")
                fam_iv[c][g].append((int(p[5]), int(p[6])))

        # total repeat content: merge ACROSS types so overlapping hits of
        # different families are not counted twice (nested hits are already out)
        all_iv = defaultdict(list)
        for c, byfam in fam_iv.items():
            for g, iv in byfam.items():
                all_iv[c].extend(iv)

        for c, k in keep.items():
            if k is None:
                continue
            n = L[c]
            glen[iso][k] += n
            half_kb = n / 2000.0                          # nearest-half limit, in kb
            # exposure: every chromosome of the class contributes
            for w in widths:
                st = starts[w]
                ok = (st + w) <= half_kb
                key = (iso, k, w)
                if key not in EXP:
                    EXP[key] = np.zeros(nw_of[w])
                EXP[key][ok] += w * 1000.0
            series = dict(fam_iv[c])
            if all_iv.get(c):
                series["ALL"] = all_iv[c]
            for g, iv in series.items():
                gfam[iso][k][g] += sum(e - s + 1 for s, e in iv)
                P = cum_cov(iv, n, args.res_bp)
                last = len(P) - 1
                for w in widths:
                    st = starts[w]
                    ok = (st + w) <= half_kb
                    lo = np.minimum(((st[ok] * 1000) // args.res_bp).astype(int), last)
                    hi = np.minimum((((st[ok] + w) * 1000) // args.res_bp).astype(int), last)
                    bp = (P[hi] - P[lo]).astype(np.float64)
                    key = (iso, k, g, w)
                    if key not in BP:
                        BP[key] = np.zeros(nw_of[w])
                    BP[key][ok] += bp

    print(f"precomputed exposure for {len(EXP)} (isolate, class, width) series and "
          f"{len(BP)} (isolate, class, type, width) series")
    return EXP, BP, glen, gfam, starts, nw_of


def main():
    args = parse_args()
    os.makedirs(OUT, exist_ok=True)
    EXP, BP, glen, gfam, starts, nw_of = load_and_precompute(args)

    # ---- which types to profile
    tot = defaultdict(lambda: defaultdict(int))
    for iso in ISOLATES:
        for k in ("core", "accessory"):
            for g, v in gfam[iso][k].items():
                tot[k][g] += v
    ab = pd.DataFrame({k: {g: 100 * v / sum(glen[i][k] for i in ISOLATES)
                           for g, v in tot[k].items()} for k in ("core", "accessory")}).fillna(0)
    ab["max"] = ab.max(axis=1)
    ab = ab.sort_values("max", ascending=False)
    fams = list(ab[(ab.core >= args.min_pct) | (ab.accessory >= args.min_pct)].index)
    fams = [f for f in fams if f != "ALL"]        # profiled separately, as the overall curve
    print(f"\ntypes profiled ({len(fams)}): " + ", ".join(fams))
    print(ab.loc[fams].round(2).to_string())
    ab.to_csv(os.path.join(OUT, "repeat_density_abundance.tsv"), sep="\t")

    widths = [float(x) for x in args.widths_kb.split(",")]
    xb = np.arange(0, args.max_kb + args.bin_kb, args.bin_kb)
    xc = (xb[:-1] + xb[1:]) / 2

    def enrich(iso_pick, k, g, w, fam_total, class_len):
        """observed/expected content of type g at each x-bin, for one resample

        numerator   bp of g inside the window / window length   (pooled over windows
                    landing in the same x-bin, so overlapping windows are handled)
        denominator the same type's content over the whole genome of this class
        Exposure is summed over ALL sampled chromosomes of the class, so
        chromosomes that do not carry the type stay in the denominator.
        """
        if fam_total <= 0 or class_len <= 0:
            return np.full(len(xc), np.nan)
        b = np.zeros(nw_of[w])
        e = np.zeros(nw_of[w])
        for iso in iso_pick:
            kb = (iso, k, w)
            if kb in EXP:
                e += EXP[kb]
            kf = (iso, k, g, w)
            if kf in BP:
                b += BP[kf]
        st = starts[w]
        j = (st // args.bin_kb).astype(int)
        num = np.zeros(len(xc))
        den = np.zeros(len(xc))
        for jj, bb, ee in zip(j, b, e):
            if 0 <= jj < len(xc) and ee > 0:
                num[jj] += bb
                den[jj] += ee
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(den > 0,
                            (num / np.where(den > 0, den, 1)) / (fam_total / class_len),
                            np.nan)

    def content(iso_pick, k, g, w):
        """absolute content of type g inside the window, as a %% of window length"""
        b = np.zeros(nw_of[w])
        e = np.zeros(nw_of[w])
        for iso in iso_pick:
            kb = (iso, k, w)
            if kb in EXP:
                e += EXP[kb]
            kf = (iso, k, g, w)
            if kf in BP:
                b += BP[kf]
        st = starts[w]
        j = (st // args.bin_kb).astype(int)
        num = np.zeros(len(xc))
        den = np.zeros(len(xc))
        for jj, bb, ee in zip(j, b, e):
            if 0 <= jj < len(xc) and ee > 0:
                num[jj] += bb
                den[jj] += ee
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(den > 0, 100.0 * num / np.where(den > 0, den, 1), np.nan)

    rng = np.random.default_rng(args.seed)
    reps = {k: {g: [] for g in fams} for k in ("core", "accessory")}
    ov = {k: {"abs": [], "enr": []} for k in ("core", "accessory")}
    print(f"\nbootstrap: {args.boot} replicates")
    for b_i in range(args.boot):
        pick = list(rng.choice(ISOLATES, size=len(ISOLATES), replace=True))
        for k in ("core", "accessory"):
            Lsum = sum(glen[i][k] for i in pick) or 1
            for g in fams:
                ftot = sum(gfam[i][k].get(g, 0) for i in pick)
                # EVERY width is pushed into the replicate list, so the percentile
                # calculation pools the window-width choice as the docstring
                # promises. Using only widths[0] here (the original bug) hid the
                # window-width uncertainty from the envelope entirely.
                for w in widths:
                    reps[k][g].append(enrich(pick, k, g, w, ftot, Lsum))
            atot = sum(gfam[i][k].get("ALL", 0) for i in pick)
            for w in widths:
                ov[k]["abs"].append(content(pick, k, "ALL", w))
                ov[k]["enr"].append(enrich(pick, k, "ALL", w, atot, Lsum))

    stats = {k: {} for k in ("core", "accessory")}
    for k in ("core", "accessory"):
        for g in fams:
            A = np.vstack(reps[k][g])
            stats[k][g] = {"median": np.nanpercentile(A, 50, axis=0),
                           "lo": np.nanpercentile(A, 2.5, axis=0),
                           "hi": np.nanpercentile(A, 97.5, axis=0)}

    rd = {"x_kb": xc}
    for g in fams:
        for k, tag in (("core", "CC"), ("accessory", "AC")):
            s = stats[k][g]
            rd[f"{g}|{tag}_median"] = s["median"]
            rd[f"{g}|{tag}_lo"] = s["lo"]
            rd[f"{g}|{tag}_hi"] = s["hi"]
    pd.DataFrame(rd).to_csv(os.path.join(OUT, "repeat_density_profile.tsv"),
                            sep="\t", index=False)

    with open(os.path.join(OUT, "repeat_density_summary.md"), "w") as fh:
        fh.write("# Repeat-TYPE composition vs distance from the chromosome end\n\n")
        fh.write("Enrichment = (type content in the window) / (that type's content over the "
                 "whole genome of the SAME class). 1.0 = behaves as it does everywhere; the "
                 "class-specific denominator cancels the abundance difference between types, "
                 "so CC and AC are comparable.\n\n")
        fh.write(f"Window widths pooled: {args.widths_kb} kb, step {args.step_kb:g} kb, "
                 f"x-bin {args.bin_kb:g} kb, coverage resolution {args.res_bp} bp. "
                 f"95% envelope from {args.boot} bootstrap replicates resampling the "
                 f"{len(ISOLATES)} chromosome-level assemblies.\n\n")
        fh.write("## Genome-wide abundance (% of class length)\n\n")
        fh.write("| type | CC % | AC % |\n|---|---|---|\n")
        for g in fams:
            fh.write(f"| {g} | {ab.loc[g,'core']:.2f} | {ab.loc[g,'accessory']:.2f} |\n")
        fh.write("\nAn enrichment whose 95% envelope includes 1 is not distinguishable "
                 "from the genome-wide level for that class.\n")

    # ---- overall repeat content (the successor to the fixed-100-kb comparison)
    ovstats = {k: {} for k in ("core", "accessory")}
    for k in ("core", "accessory"):
        for key in ("abs", "enr"):
            A = np.vstack(ov[k][key])
            ovstats[k][key] = {"median": np.nanpercentile(A, 50, axis=0),
                               "lo": np.nanpercentile(A, 2.5, axis=0),
                               "hi": np.nanpercentile(A, 97.5, axis=0)}
    od = {"x_kb": xc}
    for k, tag in (("core", "CC"), ("accessory", "AC")):
        for key, nm in (("abs", "content_pct"), ("enr", "enrichment")):
            v = ovstats[k][key]
            od[f"ALL|{tag}_{nm}_median"] = v["median"]
            od[f"ALL|{tag}_{nm}_lo"] = v["lo"]
            od[f"ALL|{tag}_{nm}_hi"] = v["hi"]
    pd.DataFrame(od).to_csv(os.path.join(OUT, "repeat_density_overall.tsv"),
                            sep="\t", index=False)

    print("\noverall repeat content vs distance (median [95% envelope]):")
    print(f"{'kb':>5}  {'CC content %':>21}  {'AC content %':>21}  "
          f"{'CC x genome':>19}  {'AC x genome':>19}")
    for kb in (5, 25, 50, 100, 250, 450):
        jj = int(np.argmin(np.abs(xc - kb)))
        row = f"{xc[jj]:5.0f} "
        for k in ("core", "accessory"):
            a = ovstats[k]["abs"]
            row += f"  {a['median'][jj]:7.3f} [{a['lo'][jj]:6.3f},{a['hi'][jj]:6.3f}]"
        for k in ("core", "accessory"):
            e2 = ovstats[k]["enr"]
            if np.isfinite(e2["median"][jj]):
                row += f"  {e2['median'][jj]:6.2f} [{e2['lo'][jj]:5.2f},{e2['hi'][jj]:5.2f}]"
            else:
                row += f"  {'n/a':>17}"
        print(row)
    with open(os.path.join(OUT, "repeat_density_summary.md"), "a") as fh:
        fh.write("\n## Overall repeat content vs distance from the end\n\n")
        fh.write("Absolute content is the %% of the window occupied by any repeat; "
                 "the second pair is the same value over each class's own "
                 "genome-wide repeat content (1.0 = flat).\n\n")
        fh.write("| distance (kb) | CC content % [95%] | AC content % [95%] | "
                 "CC x genome [95%] | AC x genome [95%] |\n|---|---|---|---|---|\n")
        for kb in (5, 25, 50, 100, 250, 450):
            jj = int(np.argmin(np.abs(xc - kb)))
            ca, aa = ovstats["core"]["abs"], ovstats["accessory"]["abs"]
            ce, ae = ovstats["core"]["enr"], ovstats["accessory"]["enr"]
            fh.write(f"| {xc[jj]:.0f} | {ca['median'][jj]:.3f} "
                     f"[{ca['lo'][jj]:.3f}, {ca['hi'][jj]:.3f}] | "
                     f"{aa['median'][jj]:.3f} [{aa['lo'][jj]:.3f}, {aa['hi'][jj]:.3f}] | "
                     f"{ce['median'][jj]:.2f} [{ce['lo'][jj]:.2f}, {ce['hi'][jj]:.2f}] | "
                     f"{ae['median'][jj]:.2f} [{ae['lo'][jj]:.2f}, {ae['hi'][jj]:.2f}] |\n")

    # Each class's genome-wide repeat content, i.e. the level the profile returns to.
    # Computed from the same gfam/glen totals that panel b uses as its denominator, so
    # the dashed line and the 1.0 line of panel b mean the same thing.
    # clearance between a baseline label and whatever is above it, in y-axis units
    LABEL_LIFT = 2.5
    base_abs = {k: 100.0 * sum(gfam[i][k].get("ALL", 0) for i in ISOLATES)
                / sum(glen[i][k] for i in ISOLATES)
                for k in ("core", "accessory")}
    print("genome-wide repeat content (baseline drawn in panel a): "
          + "  ".join(f"{k} {v:.2f}%" for k, v in base_abs.items()))

    fig, axes = plt.subplots(1, 2, figsize=(184 / 25.4, 96 / 25.4))
    ax = axes[0]
    for k, tag, col in (("core", "CC", C_CC), ("accessory", "AC", C_AC)):
        a = ovstats[k]["abs"]
        m = np.isfinite(a["median"])
        ax.fill_between(xc[m], a["lo"][m], a["hi"][m], color=col, alpha=0.20, lw=0)
        ax.plot(xc[m], a["median"][m], color=col, lw=1.6, label=tag)
        # The class's own genome-wide level, labelled at the right edge. The label is
        # lifted clear of the line: both profiles sit BELOW their own baseline at this
        # end of the axis (the end that they decay back to), so raising the text moves
        # it away from the curve rather than into it.
        ax.axhline(base_abs[k], color=col, ls="--", lw=0.8, zorder=1)
        # A FIXED lift is not enough: near the right edge the profile bumps back up over
        # its own baseline (the tail is noisy), so a constant offset grazes the curve.
        # The label is placed above that profile's maximum over the band the TEXT
        # occupies, not over a guessed band: right-aligned at the edge, ~35 pt of text
        # spans roughly the rightmost quarter of the axis. An earlier 10% band was
        # narrower than the label and left a 1 pt overlap.
        band = np.isfinite(a["median"]) & (xc >= args.max_kb * 0.65)
        top = float(np.nanmax(a["median"][band])) if band.any() else base_abs[k]
        lab_y = max(base_abs[k], top) + LABEL_LIFT
        ax.text(args.max_kb * 0.995, lab_y, f"{tag} mean {base_abs[k]:.1f}%",
                color=col, fontsize=7.4, va="bottom", ha="right", zorder=3,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.set_xlim(0, args.max_kb)
    ax.set_xticks(range(0, int(args.max_kb) + 1, 100))
    ax.set_xlabel("distance from chromosome end (kb)", fontsize=8.5)
    ax.set_ylabel("repeat content (% of window)", fontsize=8.5)
    # upper right, not center right: the centre-right position sat on the AC curve
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("a  Overall repeat content", fontsize=9, loc="left")
    ax.tick_params(labelsize=8)
    ax = axes[1]
    ax.axhline(1.0, color="grey", lw=0.9, ls="--")
    for k, tag, col in (("core", "CC", C_CC), ("accessory", "AC", C_AC)):
        e2 = ovstats[k]["enr"]
        m = np.isfinite(e2["median"])
        ax.fill_between(xc[m], e2["lo"][m], e2["hi"][m], color=col, alpha=0.20, lw=0)
        ax.plot(xc[m], e2["median"][m], color=col, lw=1.6, label=tag)
    ax.set_xlim(0, args.max_kb)
    ax.set_xticks(range(0, int(args.max_kb) + 1, 100))
    ax.set_xlabel("distance from chromosome end (kb)", fontsize=8.5)
    ax.set_ylabel("over its genome-wide level", fontsize=8.5)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("b  Relative to each class's own level", fontsize=9, loc="left")
    ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.085, right=0.965, bottom=0.155, top=0.875, wspace=0.32)
    fig.savefig(os.path.join(OUT, "fig_s31_repeat_density_overall.pdf"))
    fig.savefig(os.path.join(OUT, "fig_s31_repeat_density_overall.png"), dpi=300)
    plt.close(fig)

    # ---- figure 1: per-type profiles
    panels = args.pool_with.split(",") if args.pool_with else fams[:args.n_panels]
    ncol = 2
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(180 / 25.4, 74 * nrow / 25.4),
                             sharex=True, squeeze=False)
    for ax, g in zip(axes.ravel(), panels):
        ax.axhline(1.0, color="grey", lw=0.8, ls="--")
        for k, tag, col in (("core", "CC", C_CC), ("accessory", "AC", C_AC)):
            s = stats[k][g]
            m = np.isfinite(s["median"])
            ax.fill_between(xc[m], s["lo"][m], s["hi"][m], color=col, alpha=0.20, lw=0)
            ax.plot(xc[m], s["median"][m], color=col, lw=1.4, label=tag)
        ax.set_yscale("log")
        ax.set_ylim(0.05, 200)
        ax.set_yticks([0.1, 1, 10, 100])
        ax.set_yticklabels(["0.1", "1", "10", "100"], fontsize=8)
        ax.set_xlim(0, args.max_kb)
        ax.set_xticks(range(0, int(args.max_kb) + 1, 100))
        ax.tick_params(axis="x", labelsize=8)
        ax.set_title(f"{g}   (CC {ab.loc[g,'core']:.2f}% · AC {ab.loc[g,'accessory']:.2f}%)",
                     fontsize=8.6, loc="left")
        if g == panels[0]:
            ax.legend(loc="upper right", fontsize=7.5)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    for ax in axes[-1]:
        # kept short: the label is centred on the RIGHT-hand panel, so a long
        # one runs past the canvas edge (the margin checker catches exactly this)
        ax.set_xlabel("distance from end (kb)", fontsize=8)
    for ax in axes[:, 0]:
        ax.set_ylabel("enrichment over\ngenome-wide (log)", fontsize=8)
    # right=0.95, not 0.985: the outermost x tick label ("500") is centred on the
    # axes' right edge, so a tight right margin pushes it past the canvas
    fig.subplots_adjust(left=0.105, right=0.950, bottom=0.085, top=0.955,
                        hspace=0.40, wspace=0.26)
    fig.savefig(os.path.join(OUT, "fig_s32_repeat_density_byfamily.pdf"))
    fig.savefig(os.path.join(OUT, "fig_s32_repeat_density_byfamily.png"), dpi=300)
    plt.close(fig)

    print(f"\nwrote 3 figures + 3 tables to {os.path.relpath(OUT, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
