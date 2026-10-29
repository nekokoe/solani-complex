#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gene transfer vs repeat mediation - 211-genome clustering version.

Transfer detection uses an mmseqs protein clustering over ALL 211 FSSC
isolates (--clusters; either the archived 100%-identity file or the locally
re-run 99% one), so a gene's
"non-transferred" status is reliable (not limited to 13 genomes). The
flanking analysis is restricted to the 13 TGS assemblies (the only genomes
with RepeatMasker annotations and accessory-chromosome calls).

Definitions:
  AC-related transfer candidate = gene on the 13 TGS assemblies whose cluster
    (over all 211 genomes) spans >=2 distinct chromosomes AND involves at
    least one accessory chromosome of the 13 TGS isolates.
    Subtypes: A-A (all cluster chromosomes accessory), C-A (core + accessory).
  Background = all other genes on the 13 TGS assemblies (whole-genome
    comparison, as requested).

Outputs (06.accessory-chrs/gene_transfer_211/):
  suppl_table.gene_transfer_211.tsv           per-gene table (13 TGS genes)
  suppl_table.repeat_enrichment_211.tsv       per-family flanking enrichment
                                              vs whole genome, ranked
  suppl_table.repeat_enrichment_211.fine.tsv  per fine-family enrichment
  suppl_table.transfer_summary_211.tsv        counts
"""

import argparse
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]

FAM_GROUPS = {
    "Unknown": "Unknown", "LTR/Gypsy": "LTR", "LTR/Copia": "LTR",
    "LTR/Pao": "LTR", "LTR/Ngaro": "LTR", "DNA/TcMar-Fot1": "DNA",
    "DNA/TcMar-Tc1": "DNA", "DNA/TcMar-Tc2": "DNA", "DNA/TcMar-Ant1": "DNA",
    "DNA/MULE-MuDR": "DNA", "DNA/hAT-Restless": "DNA", "DNA/hAT-Ac": "DNA",
    "DNA/hAT-Tag1": "DNA", "DNA/PiggyBac": "DNA", "DNA/Kolobok-H": "DNA",
    "DNA/PIF-Harbinger": "DNA", "DNA/IS3EU": "DNA", "DNA/Academ-H": "DNA",
    "RC/Helitron": "RC", "RC/Helitron-2": "RC", "LINE/Tad1": "LINE",
    "LINE/CRE": "LINE", "Simple_repeat": "Simple_repeat",
    "Low_complexity": "Low_complexity",
}


def fam_group(fam):
    if fam in FAM_GROUPS:
        return FAM_GROUPS[fam]
    for pre in ("DNA/", "LTR/", "LINE/", "RC/"):
        if fam.startswith(pre):
            return pre.rstrip("/")
    return "Unknown" if fam.startswith("Unknown") else "Other"


def parse_args():
    p = argparse.ArgumentParser(description="211-genome transfer vs repeat")
    # Any mmseqs identity level: 06.accessory-chrs/tgs_assemblies/all_solaniSC.cluster100_cluster.tsv
    # (archived, 100%) or tmp/mmseqs99/out99_cluster.tsv (re-run locally, 99%).
    p.add_argument("--clusters", required=True)
    p.add_argument("--isolate-list", required=True)
    p.add_argument("--annotation-root", required=True)
    p.add_argument("--accessory-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def build_affiliation(isolate_list, annotation_root):
    """transcript -> contig for ALL 211 isolates."""
    iso_set = set(isolate_list)
    m = {}
    for iso in iso_set:
        gff = os.path.join(annotation_root, iso, "funannotate",
                           "predict_results", "fusarium.gff3")
        with open(gff) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "mRNA":
                    continue
                tid = None
                for a in parts[8].split(";"):
                    if a.startswith("ID="):
                        tid = a[3:]
                        break
                if tid:
                    m[tid] = parts[0]
    return m


def split_member(member, iso_names):
    """(iso, transcript) for either id form.

    The archived clustering writes "ISOLATE|transcript"; a fresh mmseqs run writes the
    bare transcript id "ISOLATE_transcript". Both name the same gff3 transcript, but a
    parser that assumes one form matches NOTHING for the other without any error, so
    both are handled here. 19 isolate names contain an underscore (FSSC_6, IMI_366586,
    NRRL_20438, ...), so the isolate is taken by LONGEST-prefix match, never by
    splitting on "_".
    """
    if "|" in member:
        iso, _, trans = member.partition("|")
        return iso, trans
    for n in iso_names:
        if member.startswith(n + "_"):
            return n, member
    return None, member


def load_clusters_211(path, affiliation, tgs_set, iso_names):
    """Return gene(iso|transcript) -> dict with cluster info:
    (n_contigs, n_members, acc_related, ttype)."""
    rep_members = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            iso, trans = split_member(member, iso_names)
            contig = affiliation.get(trans)
            if contig is None:
                continue
            rep_members[rep].append((iso, trans, contig))

    info = {}
    for rep, members in rep_members.items():
        contigs = sorted({c for _, _, c in members})
        n_contig = len(contigs)
        if n_contig == 1:
            for iso, trans, c in members:
                if iso in tgs_set:
                    info[f"{iso}|{trans}"] = (1, len(members), False, "single")
            continue
        # accessory involvement: accessory chromosomes exist only on the 13 TGS
        has_acc = any(ACC_TYPE.get(c) == "accessory" for c in contigs)
        if not has_acc:
            for iso, trans, c in members:
                if iso in tgs_set:
                    info[f"{iso}|{trans}"] = (n_contig, len(members), False, "no-acc")
            continue
        has_core = any(ACC_TYPE.get(c) == "core" for c in contigs)
        ttype = "C-A" if has_core else "A-A"
        for iso, trans, c in members:
            if iso in tgs_set:
                info[f"{iso}|{trans}"] = (n_contig, len(members), True, ttype)
    return info


def load_repeats(out_path):
    rep = defaultdict(list)
    with open(out_path) as f:
        for i, line in enumerate(f):
            if i < 3:
                continue
            parts = line.split()
            if not parts or len(parts) < 15:
                continue
            if len(parts) > 15 and parts[15] == "*":
                continue
            rep[parts[4]].append((int(parts[5]), int(parts[6]),
                                  fam_group(parts[10]), parts[10]))
    for c in rep:
        rep[c].sort()
    return rep


def nearest_repeat(intervals, start, end):
    best_d, best_f = None, None
    for s, e, f, fine in intervals:
        d = 0 if (s <= end and e >= start) else (start - e if e < start else s - end)
        if best_d is None or d < best_d:
            best_d, best_f = d, f
            if d == 0:
                break
    return best_d, best_f


def flanked_by(intervals, start, end, flank):
    fs = set()
    for s, e, f, fine in intervals:
        if e < start - flank:
            continue
        if s > end + flank:
            break
        fs.add(f)
    return fs


def flanked_by_fine(intervals, start, end, flank):
    fs = set()
    for s, e, f, fine in intervals:
        if e < start - flank:
            continue
        if s > end + flank:
            break
        fs.add(fine)
    return fs


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    iso_list = [x.strip() for x in open(args.isolate_list) if x.strip()]
    tgs_set = set(ISOLATES)

    affiliation = build_affiliation(iso_list, args.annotation_root)
    print(f"affiliation (all {len(iso_list)} isolates): {len(affiliation)} transcripts")

    # contig types for the 13 TGS (only they have accessory calls)
    global ACC_TYPE
    ACC_TYPE = {}
    for iso in ISOLATES:
        acc = set(open(os.path.join(args.accessory_dir,
                                    f"{iso}.accessory_chrs.id")).read().split())
        with open(os.path.join(args.length_dir, iso, "repeatmasker",
                               f"{iso}.assembly.fasta.masked.length")) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    c, L = parts[0], int(parts[1])
                    ACC_TYPE[c] = ("accessory" if c in acc
                                   else "core" if L >= 500000 else "scaffold")

    repeats = {}
    for iso in ISOLATES:
        out = os.path.join(args.accessory_dir, f"{iso}.assembly.fasta.out")
        repeats[iso] = load_repeats(out)

    # 13-TGS genes with coordinates
    genes = []
    for iso in ISOLATES:
        gff = os.path.join(args.annotation_root, iso, "funannotate",
                           "predict_results", "fusarium.gff3")
        with open(gff) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9 or parts[2] != "mRNA":
                    continue
                tid = None
                for a in parts[8].split(";"):
                    if a.startswith("ID="):
                        tid = a[3:]
                        break
                if tid:
                    genes.append((iso, tid, parts[0], int(parts[3]),
                                  int(parts[4])))
    print(f"13-TGS genes: {len(genes)}")

    iso_names = sorted(iso_list, key=len, reverse=True)
    cluster_info = load_clusters_211(args.clusters, affiliation, tgs_set, iso_names)
    print(f"13-TGS genes with cluster info: {len(cluster_info)}")

    rows = []
    counter = Counter()
    for iso, tid, contig, gs, ge in genes:
        key = f"{iso}|{tid}"
        n_contig, n_member, acc_related, ttype = cluster_info.get(
            key, (0, 0, False, "unclustered"))
        counter[ttype] += 1
        rp = repeats[iso].get(contig, [])
        nd, nf = nearest_repeat(rp, gs, ge)
        f1 = flanked_by(rp, gs, ge, 1000)
        f5 = flanked_by(rp, gs, ge, 5000)
        f1f = flanked_by_fine(rp, gs, ge, 1000)
        f5f = flanked_by_fine(rp, gs, ge, 5000)
        rows.append({
            "Isolate": iso, "Transcript": tid, "Contig": contig,
            "Contig_type": ACC_TYPE.get(contig, "scaffold"),
            "Cluster_contigs_211": n_contig, "Cluster_members_211": n_member,
            "AC_related_transfer": acc_related, "Transfer_type": ttype,
            "Nearest_repeat_dist": nd, "Nearest_repeat_family": nf,
            "Flank1kb_families": ";".join(sorted(f1)),
            "Flank5kb_families": ";".join(sorted(f5)),
            "Flank1kb_fine": ";".join(sorted(f1f)),
            "Flank5kb_fine": ";".join(sorted(f5f)),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out_dir, "suppl_table.gene_transfer_211.tsv"),
              sep="\t", index=False)
    print("transfer types:", dict(counter))

    # ---- enrichment vs whole genome (user-requested comparison) ----
    trans = df[df["AC_related_transfer"]]
    back = df[~df["AC_related_transfer"]]
    print(f"\ntransfer candidates: {len(trans)}; whole-genome background: {len(back)}")

    grouped = ["DNA", "LTR", "LINE", "RC", "Unknown", "Simple_repeat",
               "Low_complexity", "Other"]
    summ = []
    for f in grouped:
        for col, lab in [("Flank1kb_families", "flank_1kb"),
                         ("Flank5kb_families", "flank_5kb")]:
            a = trans[col].fillna("").str.contains(f).sum()
            b = back[col].fillna("").str.contains(f).sum()
            or_, p = stats.fisher_exact([[a, len(trans) - a],
                                         [b, len(back) - b]])
            summ.append({"Family": f, "Metric": lab,
                         "Transfer_frac": a / len(trans),
                         "Background_frac": b / len(back),
                         "Odds_ratio": or_, "Fisher_p": p})
    sdf = pd.DataFrame(summ)
    sdf.to_csv(os.path.join(args.out_dir,
                            "suppl_table.repeat_enrichment_211.tsv"),
               sep="\t", index=False)
    print("\n=== grouped families, transfer vs whole genome ===")
    print(sdf.round(4).to_string())

    # fine families ranked (flank 1 kb)
    fine = defaultdict(lambda: [0, 0])
    for col in ["Flank1kb_families"]:
        pass
    fam_counts = Counter()
    for fams in df["Flank1kb_families"].fillna("").str.split(";"):
        for f in fams:
            if f:
                fam_counts[f] += 1
    fine_rows = []
    for col, met in [("Flank1kb_fine", "flank_1kb"),
                     ("Flank5kb_fine", "flank_5kb")]:
        fam_counts2 = Counter()
        for fams in df[col].fillna("").str.split(";"):
            for f in fams:
                if f:
                    fam_counts2[f] += 1
        for f in fam_counts2:
            a = trans[col].fillna("").str.contains(f).sum()
            b = back[col].fillna("").str.contains(f).sum()
            if a == 0:
                continue
            or_, p = stats.fisher_exact([[a, len(trans) - a],
                                         [b, len(back) - b]])
            fine_rows.append({"Fine_family": f, "Metric": met,
                              "Transfer_frac": a / len(trans),
                              "Background_frac": b / len(back),
                              "Odds_ratio": or_, "Fisher_p": p,
                              "n_transfer": a})
    fdf = pd.DataFrame(fine_rows)
    fdf = fdf.sort_values("Odds_ratio", ascending=False)
    fdf.to_csv(os.path.join(args.out_dir,
                            "suppl_table.repeat_enrichment_211.fine.tsv"),
               sep="\t", index=False)
    print("\n=== top 15 fine families by OR (flank 1 kb) ===")
    print(fdf.head(15).round(3).to_string())


if __name__ == "__main__":
    main()
