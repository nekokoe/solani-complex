#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S20e 滑窗分析 v2：TGS 临床 AC 沿染色体位置 10kb 滑窗 × BMU02483 共有 k-mer 比例。

- 每窗口一个独立 fasta 文件（kmer-db 要求每样本独立文件）
- kmer-db build (k=19) 建窗口库 -> one2all (query = BMU02483) 得共有 k-mer 数
- 归一化: prop = common / min(win_kmers, b2483_kmers)
"""
import os
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO

WINDOW = 10000
K = 19
OUT = "06.accessory-chrs/BMU02483_gene_content/sliding"
WIN_DIR = f"{OUT}/windows"
os.makedirs(WIN_DIR, exist_ok=True)

CLIN_ACS = [
    "BMU03267_0024", "BMU03267_0018", "BMU03267_0022", "BMU03267_0028",
    "BMU03267_0019", "BMU03267_0020", "BMU03267_0025",
    "BMU01933_0015", "BMU01933_0017", "BMU01933_0019", "BMU01933_0021",
]
BMU2483 = ("/share/pasteur/gaosh/Fusarium/phase2/pan-genome-kmer/"
           "Fusarium_solani_Species_Complex/genome/BMU02483.fasta")

# ---- 1. 切窗口为独立文件 ----
meta = []  # (ac, win_idx, start, end, frac_pos)
seqs_cache = {}
for ac in CLIN_ACS:
    iso = ac.split("_")[0]
    if iso not in seqs_cache:
        seqs_cache[iso] = {r.id: str(r.seq)
                           for r in SeqIO.parse(
                               f"06.accessory-chrs/tgs_assemblies/{iso}.accessory_chrs.fasta",
                               "fasta")}
    seq = seqs_cache[iso].get(ac)
    if not seq:
        print(f"WARN skip {ac}"); continue
    L = len(seq)
    n_w = L // WINDOW
    for i in range(n_w):
        start = i * WINDOW
        win = seq[start:start + WINDOW]
        wid = f"{ac}_w{i:03d}"
        with open(f"{WIN_DIR}/{wid}.fa", "w") as f:
            f.write(f">{wid}\n{win}\n")
        meta.append((ac, i, start, start + len(win), (start + len(win) / 2) / L))
print(f"windows written: {len(meta)}")

# ---- 2. build 窗口库 ----
slist = f"{OUT}/window_list.txt"
with open(slist, "w") as f:
    for m in meta:
        f.write(f"{WIN_DIR}/{m[0]}_w{m[1]:03d}.fa\n")
db = f"{OUT}/ac_windows.kmer-db"
r = subprocess.run(["kmer-db", "build", "-k", str(K), "-t", "8", slist, db],
                   capture_output=True, text=True)
print("build rc:", r.returncode)
if r.returncode != 0:
    print(r.stderr[-500:]); raise SystemExit(1)

# ---- 3. one2all ----
common = f"{OUT}/common_2483.txt"
r = subprocess.run(["kmer-db", "one2all", db, BMU2483, common],
                   capture_output=True, text=True)
print("one2all rc:", r.returncode)
if r.returncode != 0:
    print(r.stderr[-500:]); raise SystemExit(1)

# ---- 4. 归一化 ----
# common table 3 行: 第1行 header, 第2行 db-samples, 第3行 query (第1列=query 文件名, 第2列起=每窗口共有 kmer)
raw = pd.read_csv(common, header=None)
query_row = raw.iloc[2]
vals = query_row.iloc[1:].astype(float).tolist()
b2483_len = sum(len(r) for r in SeqIO.parse(BMU2483, "fasta"))
b2483_kmers = max(1, b2483_len - K + 1)
props = []
for (ac, i, start, end, fp), v in zip(meta, vals):
    wk = max(1, (end - start) - K + 1)
    props.append(v / min(wk, b2483_kmers))
df = pd.DataFrame({
    "AC": [m[0] for m in meta],
    "window": [m[1] for m in meta],
    "start": [m[2] for m in meta],
    "end": [m[3] for m in meta],
    "frac_pos": [m[4] for m in meta],
    "common_kmer": vals,
    "prop": props,
})
df.to_csv(f"{OUT}/sliding_2483_prop.tsv", sep="\t", index=False)
agg = df.groupby("AC")["prop"].agg(["mean", "max", "count"])
print("\nper-AC sliding stats:")
print(agg.round(3).to_string())
print(f"\nsaved -> {OUT}/sliding_2483_prop.tsv")
