#!/usr/bin/env python3
"""Where do the clinical-associated candidate orthogroups sit in the pan-genome, and
how does their functional composition compare with the two ends of that gradient?

Answering two questions the Results subsection asks:

  1 the distribution of the candidate HOGs over the five pan-genome presence intervals.
    Authoritative source = `05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv`
    (presence rate over the 211 FSSC genomes, 29,080 orthogroups). Earlier candidate tables
    binned against different universes (12,842 variable orthogroups, or a species cohort)
    and therefore disagreed; this script uses the single 211-genome definition.

  2 the functional composition of the candidates versus the 0-20% interval and the
    80-100% interval of the same pan-genome, using the identical Pfam classifier that
    produced the 0-20% analysis (`classify_pfam` from scripts/stats_020_genes.py).

Outputs  07.lethal_groups/genomewide_assoc/candidate_presence_intervals.tsv
         07.lethal_groups/genomewide_assoc/candidate_function_vs_bins.tsv
"""

import os
import sys

import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from stats_020_genes import classify_pfam                      # noqa: E402

IV = os.path.join(ROOT, "05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv")
ANN = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv")
PAV = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
OUT = os.path.join(ROOT, "07.lethal_groups/genomewide_assoc")
BINS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]


def main():
    iv = pd.read_csv(IV, sep="\t")
    iv.columns = ["HOG", "rate", "bin"]
    bg = iv.bin.value_counts().reindex(BINS)
    T = int(bg.sum())
    print(f"pan-genome background: {T} orthogroups over 211 genomes")
    for b in BINS:
        print(f"   {b:9s} {bg[b]:6d}  {bg[b] / T * 100:5.1f}%")

    rows = []
    cands = {}
    for sp, label in (("solani", "F. solani"), ("falciformis", "F. falciformis")):
        ids = pd.read_csv(f"{PAV}/rf_pav_{sp}_features.tsv", sep="\t").HOG_ID
        cands[label] = set(ids)
        b = iv.set_index("HOG").bin.reindex(ids)
        assert b.notna().all(), f"{int(b.isna().sum())} candidates have no interval"
        for binname in BINS:
            n = int((b == binname).sum())
            rows.append(dict(dataset=label, interval=binname, n=n,
                             pct=round(n / len(ids) * 100, 1),
                             background_pct=round(bg[binname] / T * 100, 1),
                             fold=round((n / len(ids)) / (bg[binname] / T), 2)))
    dist = pd.DataFrame(rows)
    print("\n" + dist.pivot(index="interval", columns="dataset", values="pct")
          .reindex(BINS).to_string())
    dist.to_csv(f"{OUT}/candidate_presence_intervals.tsv", sep="\t", index=False)

    # ---------- functional composition vs the two ends of the gradient ----------
    ann = pd.read_csv(ANN, sep="\t", usecols=["HOG_ID", "Pfam_Domain"])
    ann = ann.drop_duplicates("HOG_ID").set_index("HOG_ID")["Pfam_Domain"].dropna()
    cls = ann.apply(lambda v: classify_pfam(str(v).split(";")))
    binh = iv.set_index("HOG")["bin"]
    groups = {
        "pan 0%-20%": set(binh[binh == "0%-20%"].index),
        "pan 80%-100%": set(binh[binh == "80%-100%"].index),
        "F. solani candidates": cands["F. solani"],
        "F. falciformis candidates": cands["F. falciformis"],
    }
    comp = pd.DataFrame({k: cls.reindex(sorted(v)).value_counts()
                         for k, v in groups.items()}).fillna(0)
    comp_pct = (comp / comp.sum() * 100).round(1)
    print("\nPfam functional class (%):")
    print(comp_pct.sort_values("F. solani candidates", ascending=False).head(12).to_string())

    # 候选合并集 vs 0-20% / 80-100% 的 2xk 检验（只取三类都 ≥3 的类别）
    both = cands["F. solani"] | cands["F. falciformis"]
    comp_b = pd.DataFrame({"candidates": cls.reindex(sorted(both)).value_counts(),
                           "pan 0%-20%": comp["pan 0%-20%"],
                           "pan 80%-100%": comp["pan 80%-100%"]}).fillna(0)
    tests = {}
    for other in ("pan 0%-20%", "pan 80%-100%"):
        keep = comp_b[(comp_b.candidates >= 3) & (comp_b[other] >= 3)]
        if len(keep) >= 2:
            p = stats.chi2_contingency(keep[["candidates", other]].values)[1]
            tests[other] = (len(keep), p)
            print(f"\ncandidates (n={len(both)}) vs {other}: {len(keep)} shared classes, "
                  f"chi2 P = {p:.3g}")
    comp_b["pct_candidates"] = (comp_b.candidates / comp_b.candidates.sum() * 100).round(1)
    comp_b["pct_0_20"] = (comp_b["pan 0%-20%"] / comp_b["pan 0%-20%"].sum() * 100).round(1)
    comp_b["pct_80_100"] = (comp_b["pan 80%-100%"] / comp_b["pan 80%-100%"].sum() * 100).round(1)
    comp_b.reset_index().rename(columns={"index": "PfamClass"}).to_csv(
        f"{OUT}/candidate_function_vs_bins.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT}/candidate_presence_intervals.tsv and candidate_function_vs_bins.tsv")


if __name__ == "__main__":
    main()
