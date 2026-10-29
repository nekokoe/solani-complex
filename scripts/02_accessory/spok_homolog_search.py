#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
C-2: do the 13 chromosome-level FSSC assemblies carry spok (spore killer
meiotic driver) homologs, and if so, are they on accessory chromosomes?

Background
----------
Puebla Planas et al. (bioRxiv 2026.10.04.756499) show that in *F. oxysporum*
Fol4287 a spok gene sits inside a giant Starship transposon on a dispensable
chromosome and maintains it. They identify Fol4287 Spok by homology to **Spok3
of Podospora anserina** (tripartite spore-killer domain structure; the
killing-essential lysine at position 240 is conserved). A survey cited there
reports spok homologs accumulating in accessory regions of several Fusarium
lineages.

WHY NOT USE A PFAM DOMAIN
-------------------------
An earlier attempt used Pfam PF09207 ("Yeast-kill-tox", the K1/K2 killer toxin
of *S. cerevisiae*) as a spok proxy. That is the WRONG family: spok belongs to
the spore-killer meiotic drivers, not to the yeast killer toxins. For the record
PF09207 occurs in only 6 genes across the whole 211-isolate panel and in none of
the 13 chromosome-level assemblies -- but that fact does not answer this
question, so the answer is produced here by sequence homology instead.

Query set: 06.accessory-chrs/starship_ac/spok_reference.faa
    SPOK3_PODAS  A0A516F180   <- the paper's own reference
    SPOK1_PODCO  A0A447CCJ8
    SPOK2_PODAN  B2AFA8
    SPOK4_PODAS  A0A516F183
    SPOK_PODCOMM A0AAV9FY09   (Podospora aff. communis)
    SPOK_PODCO2  A0ABY6S949   (Podospora comata)

GATE: the query is first searched against ITSELF. Every query protein must have
itself as its best hit at ~100% identity. If that fails, diamond or the DB is
broken and no downstream number means anything.

Targets: 02.annotation/per_isolates/<ISO>/funannotate/predict_results/fusarium.proteins.fa
Contig classes: same convention as the rest of the pipeline
    accessory : contig in 06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id
    core      : other contig with masked length >= --core-min-length
    scaffold  : the rest

Outputs (06.accessory-chrs/starship_ac/)
    spok_hits.tsv          every hit passing --min-identity/--min-qcov
    spok_summary.md        human-readable summary
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]

OUTFMT = "6 qseqid sseqid pident length qlen slen evalue bitscore"


def parse_args():
    p = argparse.ArgumentParser(description="spok homolog search in the 13 FSSC assemblies")
    p.add_argument("--annotation-root",
                   default=os.path.join(REPO, "02.annotation", "per_isolates"))
    p.add_argument("--accessory-dir",
                   default=os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies"))
    p.add_argument("--out-dir",
                   default=os.path.join(REPO, "06.accessory-chrs", "starship_ac"))
    p.add_argument("--query", default=os.path.join(
        REPO, "06.accessory-chrs", "starship_ac", "spok_reference.faa"))
    p.add_argument("--diamond", default="diamond")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument("--min-identity", type=float, default=25.0)
    p.add_argument("--min-qcov", type=float, default=40.0)
    p.add_argument("--core-min-length", type=int, default=500000)
    p.add_argument("--skip-gate", action="store_true")
    return p.parse_args()


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"command failed ({r.returncode}): {' '.join(cmd)}\n{r.stderr[-3000:]}")
    return r.stdout


def read_lengths(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                d[p[0]] = int(p[1])
    return d


def read_gene_contig(gff):
    d = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "mRNA":
                continue
            for a in p[8].split(";"):
                if a.startswith("ID="):
                    d[a[3:]] = p[0]
                    break
    return d


def read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    args = parse_args()
    for exe in (args.diamond,):
        if shutil.which(exe) is None:
            raise SystemExit(f"required tool not found on PATH: {exe}")
    if not os.path.exists(args.query):
        raise SystemExit(f"query FASTA not found: {args.query}")
    os.makedirs(args.out_dir, exist_ok=True)

    q_ids = read_fasta_ids(args.query)
    print(f"query: {args.query}  ({len(q_ids)} proteins: {', '.join(q_ids)})")

    tmp = tempfile.mkdtemp(prefix="spok_")
    try:
        # ---------------------------------------------------------- gate
        if not args.skip_gate:
            print("\n--- SELF-VALIDATION GATE (query vs itself) ---")
            qdb = os.path.join(tmp, "q")
            run([args.diamond, "makedb", "--in", args.query, "-d", qdb])
            out = os.path.join(tmp, "self.tsv")
            run([args.diamond, "blastp", "-q", args.query, "-d", qdb,
                 "-o", out, "--outfmt", *OUTFMT.split(), "--evalue", "1e-3",
                 "--max-target-seqs", "1", "-p", str(args.threads), "--quiet"])
            bad = []
            for line in open(out):
                f = line.rstrip("\n").split("\t")
                q, s, pid = f[0], f[1], float(f[2])
                ok = (q == s and pid >= 99.0)
                print(f"  {q:<14} -> {s:<14} pident={pid:6.1f}  {'OK' if ok else 'FAIL'}")
                if not ok:
                    bad.append(q)
            missing = set(q_ids) - {l.split("\t")[0] for l in open(out)}
            if missing:
                bad.append(f"no self-hit: {sorted(missing)}")
            if bad:
                print("\nGATE FAILED: " + "; ".join(map(str, bad)))
                print("The aligner or the database is broken; downstream numbers "
                      "would be meaningless.")
                return 1
            print("GATE PASSED -- the aligner and DB round-trip correctly\n")
        else:
            print("!! gate SKIPPED (--skip-gate)\n")

        # ------------------------------------------- build the target DB
        combined = os.path.join(tmp, "fscc13.faa")
        n_prot = 0
        with open(combined, "w") as out:
            for iso in ISOLATES:
                fa = os.path.join(args.annotation_root, iso, "funannotate",
                                  "predict_results", "fusarium.proteins.fa")
                if not os.path.exists(fa):
                    raise SystemExit(f"missing proteome: {fa}")
                n_prot += len(read_fasta_ids(fa))
                with open(fa) as f:
                    out.write(f.read())
        print(f"target DB: {len(ISOLATES)} proteomes, {n_prot} proteins")

        db = os.path.join(tmp, "db")
        run([args.diamond, "makedb", "--in", combined, "-d", db])
        hits_tsv = os.path.join(tmp, "hits.tsv")
        run([args.diamond, "blastp", "-q", args.query, "-d", db, "-o", hits_tsv,
             "--outfmt", *OUTFMT.split(), "--evalue", str(args.evalue),
             "--max-target-seqs", "20", "-p", str(args.threads), "--quiet",
             "--very-sensitive"])

        # ------------------------------------------------- classify + emit
        gene_contig, contig_class = {}, {}
        class_len = defaultdict(int)
        for iso in ISOLATES:
            ann = os.path.join(args.annotation_root, iso)
            gene_contig.update(read_gene_contig(
                os.path.join(ann, "funannotate", "predict_results", "fusarium.gff3")))
            acc = set(open(os.path.join(args.accessory_dir,
                                        f"{iso}.accessory_chrs.id")).read().split())
            L = read_lengths(os.path.join(
                ann, "repeatmasker", f"{iso}.assembly.fasta.masked.length"))
            for c, ln in L.items():
                cl = ("accessory" if c in acc
                      else "core" if ln >= args.core_min_length
                      else "scaffold")
                contig_class[c] = cl
                class_len[cl] += ln
            for c in acc:
                contig_class.setdefault(c, "accessory")

        rows, seen = [], set()
        for line in open(hits_tsv):
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            q, s, pid, aln, qlen, slen, ev, bits = f
            pid, aln, qlen, ev, bits = float(pid), int(aln), int(qlen), float(ev), float(bits)
            qcov = 100.0 * aln / qlen if qlen else 0.0
            if pid < args.min_identity or qcov < args.min_qcov:
                continue
            iso = s.split("_")[0] if s.startswith("BMU") else s.split("_")[0]
            if iso not in ISOLATES:
                iso = next((i for i in ISOLATES if s.startswith(i)), "?")
            ct = gene_contig.get(s, "")
            rows.append({"query": q, "isolate": iso, "protein": s,
                         "pident": round(pid, 1), "aln_len": aln, "qlen": qlen,
                         "qcov_pct": round(qcov, 1), "evalue": ev, "bitscore": bits,
                         "contig": ct,
                         "contig_class": contig_class.get(ct, "unmapped")})
            seen.add(s)

        import pandas as pd
        df = pd.DataFrame(rows)
        out_tsv = os.path.join(args.out_dir, "spok_hits.tsv")
        (df if len(df) else pd.DataFrame(
            columns=["query", "isolate", "protein", "pident", "aln_len", "qlen",
                     "qcov_pct", "evalue", "bitscore", "contig", "contig_class"])
         ).to_csv(out_tsv, sep="\t", index=False)

        # Deduplicate to unique genes: one FSSC protein is often the best hit of
        # several queries, so raw hit counts overstate the gene number badly.
        uniq = df.drop_duplicates(subset="protein") if len(df) else df
        uniq_tsv = os.path.join(args.out_dir, "spok_unique_genes.tsv")
        uniq.to_csv(uniq_tsv, sep="\t", index=False)
        n_uniq = len(uniq)
        if n_uniq:
            print(f"\n{len(df)} hits -> {n_uniq} UNIQUE genes "
                  f"({df.duplicated(subset='protein').sum()} redundant query-protein pairs)")

        # A count that depends on a threshold must ship its sensitivity curve.
        p_ac = class_len["accessory"] / (class_len["accessory"] + class_len["core"]) \
            if (class_len["accessory"] + class_len["core"]) else 0.0
        sens = []
        for th in [25, 30, 35, 40, 45, 50, 55, 60]:
            u = uniq[uniq.pident >= th] if n_uniq else uniq
            placed = u[u.contig_class.isin(["accessory", "core"])] if len(u) else u
            n = len(placed)
            k = int((placed.contig_class == "accessory").sum()) if n else 0
            pv = stats.binomtest(k, n, p_ac, alternative="greater").pvalue if n else 1.0
            sens.append({"pident_threshold": th, "unique_genes": len(u),
                         "on_accessory": k, "on_core": n - k,
                         "AC_pct": round(100.0 * k / n, 1) if n else 0.0,
                         "expected_on_AC": round(n * p_ac, 1),
                         "binom_p_one_sided_greater": round(pv, 4),
                         "density_ratio_AC_over_core":
                             round((k / class_len["accessory"]) /
                                   ((n - k) / class_len["core"]), 2)
                             if k and (n - k) and class_len["core"] else None})
        sens_df = pd.DataFrame(sens)
        sens_df.to_csv(os.path.join(args.out_dir, "spok_threshold_sensitivity.tsv"),
                       sep="\t", index=False)

        # ------------------------------------------------------- reporting
        lines = ["# spok (spore killer) homolog search in the 13 FSSC chromosome-level assemblies", "",
                 f"Query: {os.path.basename(args.query)} "
                 f"({', '.join(q_ids)}) -- Spok3/A0A516F180 is the reference used by the paper.", "",
                 f"Filters: pident >= {args.min_identity}%, query coverage >= {args.min_qcov}%, "
                 f"e-value <= {args.evalue}.", "",
                 f"Targets: {len(ISOLATES)} proteomes, {n_prot} proteins.", "",
                 f"Length denominator: accessory {class_len['accessory']/1e6:.1f} Mb, "
                 f"core {class_len['core']/1e6:.1f} Mb "
                 f"=> AC share {p_ac*100:.1f}%.", ""]
        if not len(df):
            lines += ["**No hit passes the filters.** No spok homolog was detected in any of the "
                      "13 chromosome-level FSSC assemblies under these thresholds.", ""]
            print("\nNO HITS pass the filters.")
        else:
            print(f"\n{len(df)} hits pass the filters:")
            print(df.sort_values("bitscore", ascending=False)
                  .head(25).to_string(index=False))
            cnt = uniq.contig_class.value_counts().to_dict()
            lines += ["## Unique spok-like genes by contig class", "",
                      f"**{len(df)} query-protein pairs reduce to {n_uniq} unique genes.** "
                      "The raw pair count must NOT be quoted: one FSSC protein is frequently the "
                      "best hit of several Spok queries.", "",
                      "| class | unique genes |", "|---|---|"]
            for k, v in cnt.items():
                lines.append(f"| {k} | {v} |")
            lines += ["", "## Threshold sensitivity (this is a threshold-defined set)", "",
                      "| pident >= | unique genes | on AC | on core | AC % | expected on AC | binom p (one-sided) | density ratio |",
                      "|---|---|---|---|---|---|---|---|"]
            for _, r in sens_df.iterrows():
                lines.append(f"| {r.pident_threshold}% | {r.unique_genes} | {r.on_accessory} | "
                             f"{r.on_core} | {r.AC_pct} | {r.expected_on_AC} | "
                             f"{r.binom_p_one_sided_greater} | {r.density_ratio_AC_over_core} |")
            lines += ["", "## Best hit per query protein", "",
                      "| query | protein | isolate | %id | qcov% | e-value | contig | class |",
                      "|---|---|---|---|---|---|---|---|"]
            for q in q_ids:
                sub = df[df["query"] == q].sort_values("bitscore", ascending=False)
                if not len(sub):
                    lines.append(f"| {q} | (none) | | | | | | |")
                    continue
                r = sub.iloc[0]
                lines.append(f"| {q} | {r.protein} | {r.isolate} | {r.pident} | "
                             f"{r.qcov_pct} | {r.evalue:.2e} | {r.contig} | {r.contig_class} |")
            lines += ["", "## Unique spok-like genes per isolate", "",
                      "| isolate | unique genes | on AC |", "|---|---|---|"]
            for iso in ISOLATES:
                s = uniq[uniq.isolate == iso]
                lines.append(f"| {iso} | {len(s)} | {int((s.contig_class == 'accessory').sum())} |")
        lines += ["", "## Interpretation guidance", "",
                  "Spore killer genes are highly divergent at the sequence level, so a negative "
                  "result at these thresholds is weaker evidence than a positive one. Report the "
                  "best hit and its identity rather than a bare yes/no, and note that the paper's "
                  "own comparison relied on alignment to Spok3 (which is the query used here).", "",
                  "Co-localization with the Starship Captain marker (DUF3435) has to be tested "
                  "separately: sharing a contig is weak evidence, because these contigs are "
                  "multi-megabase chromosomes and a Captain gene elsewhere on the same chromosome "
                  "does not imply a Starship insertion next to the spok gene.", ""]
        with open(os.path.join(args.out_dir, "spok_summary.md"), "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"\nwrote {out_tsv}, {uniq_tsv}, spok_threshold_sensitivity.tsv and spok_summary.md")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
