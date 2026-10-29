#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""T2 审稿回应 #61 — 准备 AMOVA 输入。
用 plink 将 strict_core_snps 基因型转为 012 矩阵 + 样本来源/物种分组，
供 R pegas::amova 使用。
"""
import subprocess, os, sys
import pandas as pd

BASE = "04.population"
GROUPS = {
    "solani": "solani-short-solani/solani-short-solani.raw.snps.combined.genotype.strict_core_snps.plink",
    "falciformis": "solani-short-falciformis/solani-short-falciformis.raw.snps.combined.genotype.strict_core_snps.plink",
    "keratoplastica": "solani-short-keratoplastica/solani-short-keratoplastica.raw.snps.combined.genotype.strict_core_snps.plink",
}
med_env = pd.read_csv("05.stats_pan_by_species_shorter/bmu_ref_groups.med_env.txt",
                      sep="\t", header=None, names=["sample", "source"])
species = pd.read_csv("05.stats_pan_by_species_shorter/bmu_ref_groups.species.txt",
                      sep="\t", header=None, names=["sample", "species"])

os.makedirs("tmp/amova", exist_ok=True)
for sp, prefix in GROUPS.items():
    bed = os.path.join(BASE, prefix + ".bed")
    if not os.path.exists(bed):
        print(f"skip {sp}: no bed ({bed})"); continue
    # plink -> 012 (missing as -9... use --recode A which gives 0/1/2/NA)
    out = f"tmp/amova/{sp}_geno"
    r = subprocess.run(
        ["plink", "--bfile", os.path.join(BASE, prefix), "--recode", "A",
         "--out", out, "--allow-extra-chr", "--double-id"],
        capture_output=True, text=True)
    if r.returncode != 0:
        print(f"plink failed {sp}:", r.stderr[-300:]); continue
    raw = pd.read_csv(out + ".raw", sep=r"\s+")
    geno = raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    geno = geno.set_index(raw["IID"])
    geno.to_csv(f"tmp/amova/{sp}_geno_matrix.tsv", sep="\t")
    # group file: sample -> species, source
    grp = raw[["IID"]].merge(med_env, left_on="IID", right_on="sample", how="left")
    grp = grp.merge(species, on="sample", how="left")
    grp["source2"] = grp["source"].replace("unknown", "env")  # unknown treated as env for completeness
    grp.to_csv(f"tmp/amova/{sp}_groups.tsv", sep="\t", index=False)
    print(f"{sp}: {geno.shape[0]} samples x {geno.shape[1]} SNPs; "
          f"med={grp['source'].eq('med').sum()} env={grp['source'].eq('env').sum()} "
          f"unknown={grp['source'].eq('unknown').sum()}")
