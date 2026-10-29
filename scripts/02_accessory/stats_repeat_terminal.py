#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Accessory-chromosome repeat content analysis.

For each isolate with a RepeatMasker run and an accessory-chromosome call:
  1. per-chromosome repeat content (unique, non-nested coverage) for ALL
     chromosomes (core + accessory), overall and per repeat family
  2. repeat content of the terminal 10% of each core chromosome (first and
     last 10% of the sequence), overall and per repeat family
  3. statistical tests:
       - core-chromosome terminal 10% vs accessory chromosomes (repeat content
         per chromosome; pooled Mann-Whitney + per-isolate paired Wilcoxon)
       - core chromosome (whole) vs accessory (context)
       - per repeat family: terminal 10% vs accessory (paired per isolate)

Core chromosomes = non-accessory contigs >= --core-min-length bp; tiny
unplaced scaffolds are excluded from the core-chromosome definition.

Repeat content = fraction of non-nested (non-'*') RepeatMasker hits, merged,
relative to the chromosome length (nested hits would double-count).

Inputs (per isolate):
  RepeatMasker .out  : 06.accessory-chrs/tgs_assemblies/<iso>.assembly.fasta.out
  chromosome lengths : 02.annotation/per_isolates/<iso>/repeatmasker/<iso>.assembly.fasta.masked.length
  accessory list     : 06.accessory-chrs/tgs_assemblies/<iso>.accessory_chrs.id

Outputs (--output-dir):
  suppl_table.chromosome_repeat_content.tsv    per-chromosome overall values
  suppl_table.chromosome_repeat_by_family.tsv  long-format per-chromosome per-family
  suppl_table.repeat_terminal_tests.tsv        statistical test results
"""

import argparse
import os
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

FAMILIES = ["Unknown", "LTR/Gypsy", "LTR/Copia", "LTR/Pao", "LTR/Ngaro",
            "DNA/TcMar-Fot1", "DNA/TcMar-Tc1", "DNA/TcMar-Tc2", "DNA/TcMar-Ant1",
            "DNA/MULE-MuDR", "DNA/hAT-Restless", "DNA/hAT-Ac", "DNA/hAT-Tag1",
            "DNA/PiggyBac", "DNA/Kolobok-H", "DNA/PIF-Harbinger", "DNA/IS3EU",
            "DNA/Academ-H", "RC/Helitron", "RC/Helitron-2", "LINE/Tad1",
            "LINE/CRE", "Simple_repeat", "Low_complexity", "Other"]

FAM_GROUP = {
    "Unknown": "Unknown", "LTR/Gypsy": "LTR", "LTR/Copia": "LTR",
    "LTR/Pao": "LTR", "LTR/Ngaro": "LTR",
    "DNA/TcMar-Fot1": "DNA", "DNA/TcMar-Tc1": "DNA", "DNA/TcMar-Tc2": "DNA",
    "DNA/TcMar-Ant1": "DNA", "DNA/MULE-MuDR": "DNA", "DNA/hAT-Restless": "DNA",
    "DNA/hAT-Ac": "DNA", "DNA/hAT-Tag1": "DNA", "DNA/PiggyBac": "DNA",
    "DNA/Kolobok-H": "DNA", "DNA/PIF-Harbinger": "DNA", "DNA/IS3EU": "DNA",
    "DNA/Academ-H": "DNA", "RC/Helitron": "RC", "RC/Helitron-2": "RC",
    "LINE/Tad1": "LINE", "LINE/CRE": "LINE",
    "Simple_repeat": "Simple_repeat", "Low_complexity": "Low_complexity",
}


def parse_args():
    p = argparse.ArgumentParser(description="Accessory-chromosome repeat content")
    p.add_argument("--out-dir", required=True, help="dir with *.assembly.fasta.out + *.accessory_chrs.id")
    p.add_argument("--length-dir", required=True, help="dir with <iso>/repeatmasker/<iso>...masked.length")
    p.add_argument("--isolates", default=None,
                   help="comma-separated isolate list (default: all with .out and .id)")
    p.add_argument("--output-dir", required=True)
    p.add_argument("--core-min-length", type=int, default=500000)
    p.add_argument("--terminal-fraction", type=float, default=0.10)
    p.add_argument("--terminal-length", type=int, default=100000,
                   help="fixed terminal window per chromosome end (bp)")
    p.add_argument("--species-file", default=None,
                   help="isolate<tab>species table (prefix match on isolate)")
    return p.parse_args()


def load_lengths(path):
    L = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                L[parts[0]] = int(parts[1])
    return L


def parse_out(path, chroms):
    """Return dict chrom -> list of (begin, end, is_nested, family)."""
    out = defaultdict(list)
    with open(path) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if not parts or len(parts) < 15:
                continue
            chrom = parts[4]
            if chrom not in chroms:
                continue
            begin, end = int(parts[5]), int(parts[6])
            is_nested = len(parts) > 15 and parts[15] == "*"
            fam = parts[10]
            out[chrom].append((begin, end, is_nested, fam))
    return out


def merged_bp(intervals):
    """Total bp of merged (non-overlapping) intervals."""
    if not intervals:
        return 0
    iv = sorted(intervals)
    total = 0
    cur_s, cur_e = iv[0]
    for s, e in iv[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


def family_of(fam):
    if fam in FAM_GROUP:
        return FAM_GROUP[fam]
    for prefix in ("DNA/", "LTR/", "LINE/", "RC/", "Unknown"):
        if fam.startswith(prefix):
            return prefix.rstrip("/") if prefix != "Unknown" else "Unknown"
    return "Other"


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.isolates:
        isolates = [x.strip() for x in args.isolates.split(",") if x.strip()]
    else:
        isolates = sorted({f[:-len(".assembly.fasta.out")]
                           for f in os.listdir(args.out_dir)
                           if f.endswith(".assembly.fasta.out")})
    isolates = [i for i in isolates if os.path.exists(
        os.path.join(args.out_dir, f"{i}.accessory_chrs.id"))]
    print(f"isolates: {len(isolates)}")

    species_map = {}
    if args.species_file:
        with open(args.species_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    species_map[parts[0]] = parts[1]
    def species_of(iso):
        for k, v in species_map.items():
            if iso == k or k.startswith(iso + "_") or iso.startswith(k.split("_")[0] + "_") or k.startswith(iso):
                return v
        return "NA"

    rows, fam_rows = [], []
    for iso in isolates:
        out_path = os.path.join(args.out_dir, f"{iso}.assembly.fasta.out")
        len_path = os.path.join(args.length_dir, iso, "repeatmasker",
                                f"{iso}.assembly.fasta.masked.length")
        acc_path = os.path.join(args.out_dir, f"{iso}.accessory_chrs.id")
        lengths = load_lengths(len_path)
        acc = set(open(acc_path).read().split())
        acc = {c for c in acc if c in lengths}
        core = {c for c, l in lengths.items()
                if c not in acc and l >= args.core_min_length}
        hits = parse_out(out_path, set(lengths) | acc)
        tfrac = args.terminal_fraction

        for chrom, L in lengths.items():
            if chrom in acc:
                ctype = "accessory"
            elif chrom in core:
                ctype = "core"
            else:
                continue  # unplaced scaffold
            non_nest = [(b, e) for b, e, n, fam in hits.get(chrom, []) if not n]
            tot_bp = merged_bp(non_nest)
            # terminal regions of core chromosomes (fixed-length window,
            # the final definition; percentage windows are size-dependent
            # and were abandoned)
            term_fix = None
            term_fix_bp = None
            if ctype == "core":
                wfix = min(args.terminal_length, L / 2)
                term_fix_iv = [(b, e) for b, e in non_nest
                               if e < wfix or b > L - wfix]
                term_fix_bp = merged_bp(term_fix_iv)
                term_fix = term_fix_bp / (2 * wfix)
            rows.append({
                "Isolate": iso, "Species": species_of(iso),
                "Chromosome": chrom, "Type": ctype,
                "Length": L, "Repeat_bp": tot_bp,
                "Repeat_pct": tot_bp / L * 100,
                "Terminal100kb_pct": (term_fix * 100 if term_fix is not None
                                      else np.nan),
                "Terminal100kb_repeat_bp": (term_fix_bp if term_fix_bp is not None
                                            else np.nan),
            })
            # per-family (unique non-nested) whole + terminal
            fam_bp = defaultdict(int)
            fam_term_fix = defaultdict(int)
            for b, e, n, fam in hits.get(chrom, []):
                if n:
                    continue
                f = family_of(fam)
                if ctype == "core" and (e < wfix or b > L - wfix):
                    fam_term_fix[f] += e - b + 1
            for f in sorted(set(list(fam_bp) + list(fam_term_fix))):
                fam_rows.append({
                    "Isolate": iso, "Species": species_of(iso),
                    "Chromosome": chrom, "Type": ctype,
                    "Family": f, "Repeat_pct": fam_bp[f] / L * 100,
                    "Terminal100kb_pct": (fam_term_fix[f] / (2 * wfix) * 100
                                          if ctype == "core" else np.nan),
                })

    df = pd.DataFrame(rows)
    df_fam = pd.DataFrame(fam_rows)
    df.to_csv(os.path.join(args.output_dir,
                           "suppl_table.chromosome_repeat_content.tsv"),
              sep="\t", index=False)
    df_fam.to_csv(os.path.join(args.output_dir,
                               "suppl_table.chromosome_repeat_by_family.tsv"),
                  sep="\t", index=False)

    # ---------------- statistical tests (fixed 100 kb terminal definition) ----
    test_rows = []
    term_vals = df.loc[df["Type"] == "core", "Terminal100kb_pct"].dropna()
    acc_vals = df.loc[df["Type"] == "accessory", "Repeat_pct"]
    core_vals = df.loc[df["Type"] == "core", "Repeat_pct"]

    def add_test(name, a, b):
        a = a.dropna().astype(float)
        b = b.dropna().astype(float)
        row = {"Comparison": name, "n_a": len(a), "n_b": len(b),
               "median_a": a.median(), "median_b": b.median(),
               "mean_a": a.mean(), "mean_b": b.mean()}
        if len(a) >= 3 and len(b) >= 3:
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            row["MannWhitney_U"] = u.statistic
            row["MannWhitney_p"] = u.pvalue
        else:
            row["MannWhitney_U"] = np.nan
            row["MannWhitney_p"] = np.nan
        test_rows.append(row)

    add_test("Core_terminal100kb_vs_Accessory (pooled)", term_vals, acc_vals)
    add_test("Core_whole_vs_Accessory (pooled)", core_vals, acc_vals)
    add_test("Core_terminal100kb_vs_Core_whole (pooled)", term_vals, core_vals)

    # per-isolate paired Wilcoxon: isolate mean terminal vs isolate mean accessory
    iso_term = df[df["Type"] == "core"].groupby("Isolate")["Terminal100kb_pct"].mean()
    iso_acc = df[df["Type"] == "accessory"].groupby("Isolate")["Repeat_pct"].mean()
    common = sorted(set(iso_term.index) & set(iso_acc.index))
    if len(common) >= 3:
        w = stats.wilcoxon(iso_term[common], iso_acc[common])
        test_rows.append({"Comparison": "Core_terminal100kb_vs_Accessory (paired by isolate)",
                          "n_a": len(common), "n_b": len(common),
                          "median_a": iso_term[common].median(),
                          "median_b": iso_acc[common].median(),
                          "mean_a": iso_term[common].mean(),
                          "mean_b": iso_acc[common].mean(),
                          "MannWhitney_U": w.statistic,
                          "MannWhitney_p": w.pvalue})
    else:
        test_rows.append({"Comparison": "Core_terminal100kb_vs_Accessory (paired by isolate)",
                          "n_a": len(common), "n_b": len(common)})

    # per-family: terminal10 vs accessory, pooled + paired
    for fam in sorted(df_fam["Family"].unique()):
        t = df_fam[(df_fam["Family"] == fam) & df_fam["Type"].eq("core") &
                   df_fam["Terminal100kb_pct"].notna()]["Terminal100kb_pct"]
        a = df_fam[(df_fam["Family"] == fam) & (df_fam["Type"] == "accessory")]["Repeat_pct"]
        if len(t) < 3 or len(a) < 3:
            continue
        u = stats.mannwhitneyu(t, a, alternative="two-sided")
        test_rows.append({"Comparison": f"Family_{fam}: terminal100kb vs accessory",
                          "n_a": len(t), "n_b": len(a),
                          "median_a": t.median(), "median_b": a.median(),
                          "mean_a": t.mean(), "mean_b": a.mean(),
                          "MannWhitney_U": u.statistic,
                          "MannWhitney_p": u.pvalue})

    tests = pd.DataFrame(test_rows)
    tests.to_csv(os.path.join(args.output_dir,
                              "suppl_table.repeat_terminal_tests.tsv"),
                 sep="\t", index=False)
    print(f"chromosomes: {len(df)} (core {len(df[df.Type=='core'])}, "
          f"accessory {len(df[df.Type=='accessory'])})")
    print(f"tables -> {args.output_dir}")


if __name__ == "__main__":
    main()
