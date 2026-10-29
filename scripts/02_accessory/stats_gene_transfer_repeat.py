#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gene transfer vs repeat-mediation analysis (13 TGS assemblies).

Question: are genes shared/transferred between chromosomes (mmseqs clusters
whose copies sit on >=2 distinct chromosomes) associated with specific repeat
families (repeat-mediated transfer)?

Transfer definition (cluster-level, members from the 13 TGS isolates):
  single           - cluster confined to one contig (background)
  A-A              - cluster spans >=2 contigs, all accessory chromosomes
  C-A              - cluster spans >=1 core + >=1 accessory chromosome
  C-C_intra        - no accessory; >=2 core contigs in the SAME isolate
                     (intra-genome core spread / translocation)
  C-C_ortho        - no accessory; >=2 core contigs but only ONE contig per
                     isolate (cross-isolate orthologs on homologous core
                     chromosomes -> NOT a transfer, background)
Transferred = A-A, C-A (only genes related to accessory chromosomes).
Background = single, C-C_ortho, C-C_intra.

For every gene: repeat context from RepeatMasker (non-nested hits):
  nearest repeat distance + family, flanking families within 1 kb / 5 kb.

Outputs (per cluster identity threshold):
  suppl_table.gene_transfer_repeat.<id>.tsv
  suppl_table.repeat_mediation_summary.<id>.tsv
  suppl_table.transfer_summary.<id>.tsv
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
    p = argparse.ArgumentParser(description="Gene transfer vs repeat mediation")
    p.add_argument("--cluster90", required=True)
    p.add_argument("--cluster100", required=True)
    p.add_argument("--affiliation", required=True)
    p.add_argument("--accessory-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--gff-dir", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def load_affiliation(path):
    m = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                m[parts[0]] = parts[2]
    return m


def load_clusters(path, affiliation, isolates, contig_type):
    """Return gene -> (cluster_contigs, cluster_members, transfer_type)."""
    iso_set = set(isolates)
    rep_members = defaultdict(list)      # rep -> [(iso, transcript, contig)]
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            iso, _, trans = member.partition("|")
            if iso not in iso_set:
                continue
            contig = affiliation.get(trans)
            if contig is None:
                continue
            rep_members[rep].append((iso, trans, contig))

    transfer = {}
    for rep, members in rep_members.items():
        contigs = sorted({c for _, _, c in members})
        n_contig = len(contigs)
        if n_contig == 1:
            ttype = "single"
        else:
            ctypes = {contig_type.get(c, "scaffold") for c in contigs}
            has_acc = "accessory" in ctypes
            has_core = "core" in ctypes
            # intra-genome spread: any isolate with members on >=2 contigs
            by_iso = defaultdict(set)
            for iso, _, c in members:
                by_iso[iso].add(c)
            intra = any(len(v) >= 2 for v in by_iso.values())
            if has_acc and not has_core:
                ttype = "A-A"
            elif has_acc and has_core:
                ttype = "C-A"
            elif has_core and intra:
                ttype = "C-C_intra"
            elif has_core:
                ttype = "C-C_ortho"
            else:
                ttype = "scaffold"
        for iso, trans, contig in members:
            transfer[f"{iso}|{trans}"] = (n_contig, len(members), ttype)
    return transfer


def load_gff(path):
    rows = []
    with open(path) as f:
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
                rows.append((tid, parts[0], int(parts[3]), int(parts[4])))
    return rows


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
                                  fam_group(parts[10])))
    for c in rep:
        rep[c].sort()
    return rep


def nearest_repeat(intervals, start, end):
    best_d, best_f = None, None
    for s, e, f in intervals:
        d = 0 if (s <= end and e >= start) else (start - e if e < start else s - end)
        if best_d is None or d < best_d:
            best_d, best_f = d, f
            if d == 0:
                break
    return best_d, best_f


def flanked_by(intervals, start, end, flank):
    fs = set()
    for s, e, f in intervals:
        if e < start - flank:
            continue
        if s > end + flank:
            break
        fs.add(f)
    return fs


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    affiliation = load_affiliation(args.affiliation)
    print(f"affiliation transcripts: {len(affiliation)}")

    contig_type, contig_len = {}, {}
    for iso in ISOLATES:
        acc = set(open(os.path.join(args.accessory_dir,
                                    f"{iso}.accessory_chrs.id")).read().split())
        with open(os.path.join(args.length_dir, iso, "repeatmasker",
                               f"{iso}.assembly.fasta.masked.length")) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    c, L = parts[0], int(parts[1])
                    contig_len[c] = L
                    contig_type[c] = ("accessory" if c in acc
                                      else "core" if L >= 500000 else "scaffold")

    repeats = {}
    for iso in ISOLATES:
        out = os.path.join(args.accessory_dir, f"{iso}.assembly.fasta.out")
        repeats[iso] = load_repeats(out)

    genes = []
    for iso in ISOLATES:
        gff = os.path.join(args.gff_dir, iso, "funannotate",
                           "predict_results", "fusarium.gff3")
        for tid, contig, start, end in load_gff(gff):
            genes.append({"isolate": iso, "transcript": tid,
                          "contig": contig, "gstart": start, "gend": end})
    print(f"genes loaded: {len(genes)}")

    for label, clpath in [("cluster90", args.cluster90),
                          ("cluster100", args.cluster100)]:
        print(f"\n=== {label} ===")
        transfer = load_clusters(clpath, affiliation, ISOLATES, contig_type)
        print(f"genes with cluster info: {len(transfer)}")
        ttype_counter = Counter()
        rows = []
        for g in genes:
            key = f"{g['isolate']}|{g['transcript']}"
            n_contig, n_member, ttype = transfer.get(key, (0, 0, "unclustered"))
            ttype_counter[ttype] += 1
            rp = repeats[g["isolate"]].get(g["contig"], [])
            nd, nf = nearest_repeat(rp, g["gstart"], g["gend"])
            flank1 = flanked_by(rp, g["gstart"], g["gend"], 1000)
            flank5 = flanked_by(rp, g["gstart"], g["gend"], 5000)
            rows.append({
                "Isolate": g["isolate"], "Transcript": g["transcript"],
                "Contig": g["contig"],
                "Contig_type": contig_type.get(g["contig"], "scaffold"),
                "Contig_length": contig_len.get(g["contig"], np.nan),
                "Cluster_contigs": n_contig, "Cluster_members": n_member,
                "Transfer_type": ttype,
                "Nearest_repeat_dist": nd,
                "Nearest_repeat_family": nf,
                "Flank1kb_families": ";".join(sorted(flank1)),
                "Flank5kb_families": ";".join(sorted(flank5)),
            })
        df = pd.DataFrame(rows)
        df.to_csv(os.path.join(args.out_dir,
                               f"suppl_table.gene_transfer_repeat.{label}.tsv"),
                  sep="\t", index=False)
        print("transfer types:", dict(ttype_counter))

        # ---- mediation summary: transferred vs background ----
        trans = df[df["Transfer_type"].isin(["A-A", "C-A"])]
        back = df[df["Transfer_type"].isin(["single", "C-C_ortho", "C-C_intra"])]
        fams = ["DNA", "LTR", "LINE", "RC", "Unknown", "Simple_repeat",
                "Low_complexity", "Other"]
        summ = []
        for f in fams:
            for col, lab in [("Flank1kb_families", "flank_1kb"),
                             ("Flank5kb_families", "flank_5kb")]:
                a = trans[col].fillna("").str.contains(f).sum()
                b = back[col].fillna("").str.contains(f).sum()
                fa, fb = a / len(trans), b / len(back)
                or_, p = stats.fisher_exact([[a, len(trans) - a],
                                             [b, len(back) - b]])
                summ.append({"Family": f, "Metric": lab,
                             "Transferred_frac": fa, "Background_frac": fb,
                             "Odds_ratio": or_, "Fisher_p": p})
        sdf = pd.DataFrame(summ)
        sdf.to_csv(os.path.join(args.out_dir,
                                f"suppl_table.repeat_mediation_summary.{label}.tsv"),
                   sep="\t", index=False)
        pd.DataFrame([{"Transfer_type": k, "n_genes": v}
                      for k, v in ttype_counter.items()]).to_csv(
            os.path.join(args.out_dir, f"suppl_table.transfer_summary.{label}.tsv"),
            sep="\t", index=False)
        print(f"\ntransferred n={len(trans)}, background n={len(back)}")
        print(sdf.round(4).to_string())


if __name__ == "__main__":
    main()
