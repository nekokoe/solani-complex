#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Characterize the 0-20% presence-interval (strain-specific/rare) HOGs of the
FSSC pan-genome.

What this script computes:

1. Carrier-source statistics for every 0-20% HOG:
   which isolates (clinical / environmental / unknown) carry it, and an
   exclusivity class: clinical-only / environmental-only / shared / unknown-only.

2. Effector enrichment: fraction of effector genes and effector HOGs per
   presence interval (0-20% ... 80-100%), using the EffectorP master table
   (effector_master.gene_level.tsv). Genes without an EffectorP prediction
   (non-secreted or from the 4 isolates lacking effectorP output) are counted
   as non-effectors.

3. Functional annotation profile per interval and for 0-20% in detail:
   - annotation rates (eggNOG / Pfam / CAZy / PHI / DFVF / KEGG / COG)
   - COG category distribution
   - top Pfam domains
   - PHI phenotype categories
   - KEGG pathway enrichment (hypergeometric test, BH correction) of the
     0-20% HOG set against the whole-pan-genome HOG background

Outputs (all written to --output-dir):
   020_genes.hog_carrier_stats.tsv        per-HOG carrier/source stats
   suppl_table.020_overview_by_interval.tsv   per-interval overview
   suppl_table.020_COG.tsv                COG distribution, 0-20% vs all HOGs
   suppl_table.020_Pfam_top.tsv           top Pfam domains, 0-20% vs all HOGs
   suppl_table.020_PHI_DFVF.tsv           PHI/DFVF hit rates and PHI phenotypes
   suppl_table.020_KEGG_enrichment.tsv    KEGG pathway enrichment (BH-corrected)
   figure_data.020.tsv                    data for the main-figure sub-panel

Usage:
    python3 stats_020_genes.py \
        --interval-file 05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv \
        --gene-vs-hog 05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv \
        --effector-master 08.pan_020_genes/effector_master.gene_level.tsv \
        --group-file bmu_ref_groups.211.txt \
        --annotation-file 05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv \
        --pathway-name-file 07.lethal_groups/pathway.list \
        --output-dir 08.pan_020_genes
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

INTERVALS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]

# Pfam domain -> functional class keywords (case-insensitive substring match on
# the Pfam_Domain entry, HET matched as prefix). Priority order decides the
# single primary class when several match; generic/unknown domains fall to
# 'Other Pfam'. NB: bare 'transp' would also match transposase/transpeptidase
# domains (mu-transpos_c, transpeptidase), so transporter keywords are anchored
# to transporter-specific tokens.
PFAM_CLASSES = [
    ("Transporter", ["mfs", "permease", "transp_", "transport", "cotransp",
                     "sugar_tr", "exchanger"]),
    ("Transcription factor", ["fungal_trans", "zn_clus", "bzip", "homeobox"]),
    ("Oxidoreductase", ["adh_", "p450", "oxido", "reductase", "dehydrogenase",
                        "monooxygenase", "nad_binding", "fad_binding", "pyr_redox",
                        "flavin", "ferredoxin", "cytochrom"]),
    ("HET domain", ["^het"]),
    ("Ankyrin repeat", ["ank"]),
    ("Protein kinase", ["pkinase", "kinase", "hiska"]),
    ("Peptidase", ["peptidase", "protease", "trypsin"]),
    ("Methyltransferase", ["methyltransf", "methylase"]),
    ("WD40/TPR repeat", ["wd40", "tpr"]),
]


def classify_pfam(domains):
    """Assign a single primary functional class to a HOG from its Pfam domains."""
    for name, kws in PFAM_CLASSES:
        for dom in domains:
            dd = dom.strip().lower()
            for kw in kws:
                if kw == "^het":
                    if dd.startswith("het"):
                        return name
                elif kw in dd:
                    return name
    return "Other Pfam"


# Standard COG functional category names (single-letter code -> description)
COG_NAMES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "X": "Mobilome: prophages, transposons",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


def parse_arguments():
    parser = argparse.ArgumentParser(description="Characterize 0-20% interval HOGs of the FSSC pan-genome")
    parser.add_argument("--interval-file", required=True, help="HOG -> presence-rate interval TSV")
    parser.add_argument("--gene-vs-hog", required=True, help="N0.gene_vs_HOG.tsv")
    parser.add_argument("--effector-master", required=True, help="effector_master.gene_level.tsv")
    parser.add_argument("--group-file", required=True, help="sample -> MedEnv mapping (bmu_ref_groups.211.txt)")
    parser.add_argument("--annotation-file", required=True, help="N0.HOG_annotation.v2.tsv")
    parser.add_argument("--pathway-name-file", required=True, help="KEGG pathway ID -> name")
    parser.add_argument("--kegg-category-file", default=None,
                        help="KEGG BRITE pathway -> category map (kegg_pathway_categories.tsv)")
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def split_set(series, sep=";"):
    """Return dict from value -> set of HOGs carrying that annotation value."""
    out = {}
    for hog, val in series.items():
        if pd.isna(val) or val == "-":
            continue
        for item in str(val).split(sep):
            item = item.strip()
            if item:
                out.setdefault(item, set()).add(hog)
    return out


def bh_correct(pvalues):
    """Benjamini-Hochberg FDR correction (monotonicity enforced on the
    sorted values, then mapped back to input order)."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q_sorted = p[order] * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n)
    q[order] = q_sorted
    return np.minimum(q, 1.0)


def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    # ---------- 1. intervals + genes ----------
    interval_df = pd.read_csv(args.interval_file, sep="\t")
    interval_df.columns = [c.strip() for c in interval_df.columns]
    interval_df = interval_df[interval_df["HOG"].str.startswith("N0.HOG", na=False)]
    interval_map = dict(zip(interval_df["HOG"], interval_df["Presence_Rate_Interval"]))
    print(f"HOGs with interval: {len(interval_map)}")

    gvh = pd.read_csv(args.gene_vs_hog, sep="\t", usecols=["GeneID", "SampleID", "HOG_ID"])
    gvh = gvh.drop_duplicates(subset=["GeneID"])
    gvh["Interval"] = gvh["HOG_ID"].map(interval_map)
    print(f"genes mapped to intervals: {gvh['Interval'].notna().sum()}")

    # ---------- 2. groups ----------
    groups = pd.read_csv(args.group_file, sep="\t")
    groups.columns = [c.strip() for c in groups.columns]
    groups = groups.rename(columns={groups.columns[0]: "Sample", "MedEnv": "MedEnv"})
    med_env = dict(zip(groups["Sample"], groups["MedEnv"]))

    # ---------- 3. effector master ----------
    eff = pd.read_csv(args.effector_master, sep="\t", usecols=["GeneID", "IsEffector", "EffectorClass"])
    eff_map = dict(zip(eff["GeneID"], eff["IsEffector"]))
    eff_class_map = dict(zip(eff["GeneID"], eff["EffectorClass"]))
    gvh["IsEffector"] = gvh["GeneID"].map(eff_map).fillna(False)
    gvh["EffectorClass"] = gvh["GeneID"].map(eff_class_map).fillna("non-effector")
    # which samples lack effectorP output (their genes are non-predicted)
    gvh["EffectorPredicted"] = gvh["GeneID"].map(lambda g: g in eff_map)

    # ---------- 4. per-HOG carrier-source stats (all intervals) ----------
    gvh["Source"] = gvh["SampleID"].map(med_env).fillna("unknown")
    carrier = (
        gvh.groupby(["HOG_ID"])
        .agg(
            Interval=("Interval", "first"),
            n_carriers=("SampleID", "nunique"),
            n_clinical=("Source", lambda s: s.eq("med").sum()),
            n_environmental=("Source", lambda s: s.eq("env").sum()),
            n_unknown=("Source", lambda s: s.eq("unknown").sum()),
        )
        .reset_index()
    )
    carrier["Clinical"] = carrier["n_clinical"] > 0
    carrier["Environmental"] = carrier["n_environmental"] > 0
    def exclusivity(row):
        if row["Clinical"] and row["Environmental"]:
            return "shared"
        if row["Clinical"]:
            return "clinical-only"
        if row["Environmental"]:
            return "environmental-only"
        return "unknown-only"
    carrier["Exclusivity"] = carrier.apply(exclusivity, axis=1)
    carrier["Presence_Rate"] = carrier["HOG_ID"].map(
        dict(zip(interval_df["HOG"], interval_df["Presence_Rate"])))
    carrier_out = os.path.join(args.output_dir, "020_genes.hog_carrier_stats.tsv")
    carrier.to_csv(carrier_out, sep="\t", index=False)
    sub_hogs = set(carrier.loc[carrier["Interval"] == "0%-20%", "HOG_ID"])
    print(f"carrier stats for {carrier['HOG_ID'].nunique()} HOGs "
          f"({len(sub_hogs)} in 0-20%) -> {carrier_out}")

    # ---------- 5. overview per interval ----------
    rows = []
    for interval in INTERVALS:
        iv = gvh[gvh["Interval"] == interval]
        n_hog = int(interval_map and sum(1 for v in interval_map.values() if v == interval))
        genes = len(iv)
        eff_genes = int(iv["IsEffector"].sum())
        eff_hogs = iv.loc[iv["IsEffector"], "HOG_ID"].nunique()
        n_eff_pred_genes = int(iv["EffectorPredicted"].sum())
        eff_by_class = iv.loc[iv["IsEffector"], "EffectorClass"].value_counts()
        rows.append({
            "Interval": interval,
            "n_HOGs": n_hog,
            "n_genes": genes,
            "n_effector_genes": eff_genes,
            "pct_effector_genes": eff_genes / genes * 100 if genes else np.nan,
            "n_effector_HOGs": eff_hogs,
            "pct_effector_HOGs": eff_hogs / n_hog * 100 if n_hog else np.nan,
            "n_genes_with_effectorP": n_eff_pred_genes,
            "pct_genes_with_effectorP": n_eff_pred_genes / genes * 100 if genes else np.nan,
            "pct_cytoplasmic": eff_by_class.get("cytoplasmic", 0) / genes * 100 if genes else np.nan,
            "pct_apoplastic": eff_by_class.get("apoplastic", 0) / genes * 100 if genes else np.nan,
            "pct_both": eff_by_class.get("both", 0) / genes * 100 if genes else np.nan,
        })
    overview = pd.DataFrame(rows)

    # annotation rates per interval
    ann = pd.read_csv(args.annotation_file, sep="\t", usecols=[
        "HOG_ID", "CAZy_Class", "COG_Class", "KEGG_Pathway", "Pfam_Domain",
        "NOG_Description", "PHI_ID", "DFVF_UniProtID"])
    ann = ann.set_index("HOG_ID")
    ann = ann[~ann.index.duplicated(keep="first")]
    for col, label in [("CAZy_Class", "CAZy"), ("COG_Class", "COG"), ("KEGG_Pathway", "KEGG"),
                       ("Pfam_Domain", "Pfam"), ("NOG_Description", "eggNOG"),
                       ("PHI_ID", "PHI"), ("DFVF_UniProtID", "DFVF")]:
        hits = ann[ann[col].notna() & (ann[col] != "-")]
        overview[f"pct_HOGs_{label}"] = overview["Interval"].map(
            lambda iv: hits.index.intersection(
                [h for h, v in interval_map.items() if v == iv]).size / sum(1 for v in interval_map.values() if v == iv) * 100)
    overview_out = os.path.join(args.output_dir, "suppl_table.020_overview_by_interval.tsv")
    overview.to_csv(overview_out, sep="\t", index=False)
    print(f"overview -> {overview_out}")

    # ---------- 6. COG distribution (0-20% vs all) ----------
    all_hogs = set(interval_map.keys())
    cog_sets = split_set(ann["COG_Class"].dropna().to_dict())
    cog_rows = []
    for letter, hogset in cog_sets.items():
        if letter not in COG_NAMES:
            continue
        n_sub = len(hogset & sub_hogs)
        n_all = len(hogset & all_hogs)
        if n_sub or n_all:
            cog_rows.append({
                "COG_letter": letter,
                "COG_category": COG_NAMES[letter],
                "n_HOGs_0_20": n_sub,
                "pct_of_0_20_HOGs": n_sub / len(sub_hogs) * 100,
                "n_HOGs_all": n_all,
                "pct_of_all_HOGs": n_all / len(all_hogs) * 100,
            })
    cog_df = pd.DataFrame(cog_rows).sort_values("n_HOGs_0_20", ascending=False)
    cog_out = os.path.join(args.output_dir, "suppl_table.020_COG.tsv")
    cog_df.to_csv(cog_out, sep="\t", index=False)
    print(f"COG -> {cog_out}")

    # ---------- 7. top Pfam domains (0-20% vs all) ----------
    pfam_sets = split_set(ann["Pfam_Domain"].dropna().to_dict())
    pfam_rows = []
    for dom, hogset in pfam_sets.items():
        n_sub = len(hogset & sub_hogs)
        if n_sub == 0:
            continue
        pfam_rows.append({
            "Pfam_domain": dom,
            "n_HOGs_0_20": n_sub,
            "pct_of_0_20_HOGs": n_sub / len(sub_hogs) * 100,
            "n_HOGs_all": len(hogset & all_hogs),
            "pct_of_all_HOGs": len(hogset & all_hogs) / len(all_hogs) * 100,
        })
    pfam_df = pd.DataFrame(pfam_rows).sort_values("n_HOGs_0_20", ascending=False).head(30)
    pfam_out = os.path.join(args.output_dir, "suppl_table.020_Pfam_top.tsv")
    pfam_df.to_csv(pfam_out, sep="\t", index=False)
    print(f"Pfam -> {pfam_out}")

    # ---------- 8. PHI / DFVF ----------
    phi_hogs = set(ann.index[ann["PHI_ID"].notna() & (ann["PHI_ID"] != "-")])
    dfvf_hogs = set(ann.index[ann["DFVF_UniProtID"].notna() & (ann["DFVF_UniProtID"] != "-")])
    phi_cat = ann[ann["PHI_ID"].notna() & (ann["PHI_ID"] != "-")]["PHI_ID"] \
        if "PHI_ID" in ann else pd.Series(dtype=object)
    pheno = split_set(ann["PHI_Category"].dropna().to_dict(), sep=";") if "PHI_Category" in ann.columns else {}
    ann2 = pd.read_csv(args.annotation_file, sep="\t", usecols=["HOG_ID", "PHI_Category"]).set_index("HOG_ID")
    pheno = split_set(ann2["PHI_Category"].dropna().to_dict())
    phi_rows = []
    for phen, hogset in sorted(pheno.items(), key=lambda kv: -len(kv[1] & sub_hogs)):
        if not (hogset & sub_hogs):
            continue
        phi_rows.append({
            "PHI_phenotype": phen,
            "n_HOGs_0_20": len(hogset & sub_hogs),
            "pct_of_0_20_PHI_HOGs": len(hogset & sub_hogs) / max(1, len(phi_hogs & sub_hogs)) * 100,
            "n_HOGs_all": len(hogset & all_hogs),
        })
    phi_df = pd.DataFrame(phi_rows)
    summary_phi = pd.DataFrame([
        {"category": "PHI hits", "n_HOGs_0_20": len(phi_hogs & sub_hogs),
         "pct_of_0_20_HOGs": len(phi_hogs & sub_hogs) / len(sub_hogs) * 100,
         "n_HOGs_all": len(phi_hogs & all_hogs),
         "pct_of_all_HOGs": len(phi_hogs & all_hogs) / len(all_hogs) * 100},
        {"category": "DFVF hits", "n_HOGs_0_20": len(dfvf_hogs & sub_hogs),
         "pct_of_0_20_HOGs": len(dfvf_hogs & sub_hogs) / len(sub_hogs) * 100,
         "n_HOGs_all": len(dfvf_hogs & all_hogs),
         "pct_of_all_HOGs": len(dfvf_hogs & all_hogs) / len(all_hogs) * 100},
    ])
    phi_out = os.path.join(args.output_dir, "suppl_table.020_PHI_DFVF.tsv")
    with open(phi_out, "w") as fh:
        summary_phi.to_csv(fh, sep="\t", index=False)
        fh.write("\n# PHI phenotype breakdown (0-20% interval)\n")
        if len(phi_df):
            phi_df.to_csv(fh, sep="\t", index=False)
    print(f"PHI/DFVF -> {phi_out}")

    # ---------- 9. KEGG pathway enrichment (hypergeometric) ----------
    # KEGG_Pathway values are comma-separated numeric IDs ("00010,00030,...")
    # Because 0-20% HOGs have a much lower KEGG annotation rate than the rest of
    # the pan-genome, the enrichment test is run on the *annotated* background
    # (HOGs with >=1 KEGG pathway), and both over- and under-representation
    # tails are reported (two-sided).
    path_sets = split_set(ann["KEGG_Pathway"].dropna().to_dict(), sep=",")
    path_names = {}
    if os.path.isfile(args.pathway_name_file):
        with open(args.pathway_name_file) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    path_names[parts[0].replace("ko", "")] = parts[1]

    annotated_all = set().union(*path_sets.values()) & all_hogs
    annotated_sub = set().union(*path_sets.values()) & sub_hogs
    M = len(annotated_all)                  # annotated HOGs in background
    n = len(annotated_sub)                  # annotated 0-20% HOGs
    kegg_rows = []
    for path, hogset in path_sets.items():
        K = len(hogset & annotated_all)     # annotated background HOGs with pathway
        if K < 5:
            continue
        k = len(hogset & annotated_sub)     # annotated 0-20% HOGs with pathway
        if k == 0:
            continue
        p_enrich = hypergeom.sf(k - 1, M, K, n)   # over-representation tail
        p_deplete = hypergeom.cdf(k, M, K, n)     # under-representation tail
        if p_enrich <= p_deplete:
            direction, pval = "enriched", p_enrich
        else:
            direction, pval = "depleted", p_deplete
        kegg_rows.append({
            "KEGG_pathway": path,
            "Pathway_name": path_names.get(path, "-"),
            "n_HOGs_0_20": k,
            "n_HOGs_all_annotated": K,
            "expected_in_0_20": n * K / M,
            "pct_of_annotated_0_20_HOGs": k / n * 100,
            "Direction": direction,
            "pvalue": pval,
        })
    kegg_df = pd.DataFrame(kegg_rows)
    if len(kegg_df):
        kegg_df["padj_BH"] = bh_correct(kegg_df["pvalue"])
        kegg_df = kegg_df.sort_values(["Direction", "padj_BH"])
    kegg_out = os.path.join(args.output_dir, "suppl_table.020_KEGG_enrichment.tsv")
    kegg_df.to_csv(kegg_out, sep="\t", index=False)
    print(f"KEGG enrichment: {len(kegg_df)} pathways tested (annotated background) -> {kegg_out}")
    print("  top enriched (by padj):")
    top_enr = kegg_df[kegg_df["Direction"] == "enriched"].head(10)
    if len(top_enr):
        print(top_enr[["KEGG_pathway", "Pathway_name", "n_HOGs_0_20",
                       "n_HOGs_all_annotated", "pvalue", "padj_BH"]].to_string(index=False))
    print("  top depleted (by padj):")
    top_dep = kegg_df[kegg_df["Direction"] == "depleted"].head(10)
    if len(top_dep):
        print(top_dep[["KEGG_pathway", "Pathway_name", "n_HOGs_0_20",
                       "n_HOGs_all_annotated", "pvalue", "padj_BH"]].to_string(index=False))

    # ---------- 10. figure data ----------
    fig_rows = []

    # exclusivity composition per interval (for pie charts)
    excl_iv = (
        carrier.groupby("Interval")["Exclusivity"]
        .value_counts()
        .unstack(fill_value=0)
        .reindex(INTERVALS)
    )
    for interval in INTERVALS:
        for cls in ["clinical-only", "environmental-only", "shared", "unknown-only"]:
            n = int(excl_iv.loc[interval].get(cls, 0)) if interval in excl_iv.index else 0
            tot = int(excl_iv.loc[interval].sum()) if interval in excl_iv.index else 0
            fig_rows.append({"panel": "exclusivity_by_interval", "interval": interval,
                             "category": cls, "n_HOGs": n,
                             "pct": n / tot * 100 if tot else 0.0})
        fig_rows.append({"panel": "exclusivity_by_interval", "interval": interval,
                         "category": "total_HOGs",
                         "n_HOGs": int(excl_iv.loc[interval].sum()) if interval in excl_iv.index else 0,
                         "pct": 100.0})

    # effector class split per interval (for stacked bars)
    for _, row in overview.iterrows():
        for cls in ["cytoplasmic", "apoplastic", "both"]:
            fig_rows.append({"panel": "effector_class_by_interval",
                             "interval": row["Interval"], "category": cls,
                             "n_HOGs": row[f"n_effector_genes"], "pct": row[f"pct_{cls}"]})
        fig_rows.append({"panel": "effector_by_interval", "interval": row["Interval"],
                         "category": "total_effector",
                         "n_HOGs": row["n_effector_genes"], "pct": row["pct_effector_genes"]})

    # carrier pairs (clinical vs environmental (HOG x isolate) pairs) per interval
    n_clin_iso = int(sum(1 for v in med_env.values() if v == "med"))
    n_env_iso = int(sum(1 for v in med_env.values() if v == "env"))
    pair_iv = carrier.groupby("Interval")[["n_clinical", "n_environmental"]].sum()
    for interval in INTERVALS:
        if interval in pair_iv.index:
            clin_p = int(pair_iv.loc[interval, "n_clinical"])
            env_p = int(pair_iv.loc[interval, "n_environmental"])
        else:
            clin_p = env_p = 0
        fig_rows.append({"panel": "carrier_pairs_by_interval", "interval": interval,
                         "category": "clinical", "n_HOGs": clin_p,
                         "pct": clin_p / (clin_p + env_p) * 100 if (clin_p + env_p) else 0.0})
        fig_rows.append({"panel": "carrier_pairs_by_interval", "interval": interval,
                         "category": "environmental", "n_HOGs": env_p,
                         "pct": env_p / (clin_p + env_p) * 100 if (clin_p + env_p) else 0.0})
    # per-isolate normalized carrier rate per interval
    iso_iv = carrier.groupby("Interval")[["n_clinical", "n_environmental"]].mean()
    for interval in INTERVALS:
        if interval in iso_iv.index:
            clin_r = float(iso_iv.loc[interval, "n_clinical"]) / n_clin_iso * 100
            env_r = float(iso_iv.loc[interval, "n_environmental"]) / n_env_iso * 100
        else:
            clin_r = env_r = 0.0
        fig_rows.append({"panel": "per_isolate_rate_by_interval", "interval": interval,
                         "category": "clinical", "n_HOGs": clin_r, "pct": clin_r})
        fig_rows.append({"panel": "per_isolate_rate_by_interval", "interval": interval,
                         "category": "environmental", "n_HOGs": env_r, "pct": env_r})

    # 0-20% COG primary functional categories (priority rule, one class per HOG)
    cog_pri = ann[ann["COG_Class"].notna() & (ann["COG_Class"] != "-")]["COG_Class"].to_dict()
    cog_pri = {h: c for h, c in cog_pri.items() if h in sub_hogs}
    COG_PRIORITY = list("JKLDVTMNZWUOCGEFHIPQABYRS")
    def primary_cog(cls):
        for letter in COG_PRIORITY:
            if letter in cls:
                return letter
        return "-"
    cog_counts = pd.Series([primary_cog(c) for c in cog_pri.values()]).value_counts()
    for letter, n in cog_counts.items():
        fig_rows.append({"panel": "cog_primary_020", "interval": "0%-20%",
                         "category": letter, "n_HOGs": int(n),
                         "pct": int(n) / len(sub_hogs) * 100})
    fig_rows.append({"panel": "cog_primary_020", "interval": "0%-20%",
                     "category": "no_COG",
                     "n_HOGs": len(sub_hogs) - sum(cog_counts.values),
                     "pct": (len(sub_hogs) - sum(cog_counts.values)) / len(sub_hogs) * 100})

    # 0-20% KEGG primary functional categories (first pathway per HOG -> BRITE level-2)
    kegg_ann = ann.loc[ann.index.isin(sub_hogs) & ann["KEGG_Pathway"].notna()
                       & (ann["KEGG_Pathway"] != "-"), "KEGG_Pathway"]
    kegg_first = kegg_ann.apply(lambda v: str(v).split(",")[0].strip())
    cat_map = {}
    if args.kegg_category_file and os.path.isfile(args.kegg_category_file):
        cat_df = pd.read_csv(args.kegg_category_file, sep="\t", dtype={"pathway_id": str})
        cat_map = {pid.zfill(5): l2 for pid, l2 in zip(cat_df["pathway_id"], cat_df["level2"])}
    def level2_of(pid):
        v = cat_map.get(pid)
        if isinstance(v, str) and " " in v:
            return v.split(" ", 1)[1]          # drop the 091xx numeric prefix
        return v or "Unmapped"
    kegg_cat = kegg_first.map(level2_of)
    n_kegg_ann = len(kegg_cat)
    for name, n in kegg_cat.value_counts().items():
        fig_rows.append({"panel": "kegg_primary_020", "interval": "0%-20%",
                         "category": name, "n_HOGs": int(n),
                         "pct": int(n) / n_kegg_ann * 100})
    fig_rows.append({"panel": "kegg_primary_020", "interval": "0%-20%",
                     "category": "no_KEGG",
                     "n_HOGs": len(sub_hogs) - n_kegg_ann,
                     "pct": (len(sub_hogs) - n_kegg_ann) / len(sub_hogs) * 100})

    # 0-20% Pfam primary functional classes (priority rule over all domains);
    # the 'Other Pfam' catch-all is merged into the 'Others' slice at plot time
    pfam_ann = ann.loc[ann.index.isin(sub_hogs) & ann["Pfam_Domain"].notna()
                       & (ann["Pfam_Domain"] != "-"), "Pfam_Domain"]
    pfam_classes = pfam_ann.apply(lambda v: classify_pfam(str(v).split(";")))
    n_pfam_ann = len(pfam_classes)
    for name, n in pfam_classes.value_counts().items():
        fig_rows.append({"panel": "pfam_primary_020", "interval": "0%-20%",
                         "category": name, "n_HOGs": int(n),
                         "pct": int(n) / n_pfam_ann * 100})
    fig_rows.append({"panel": "pfam_primary_020", "interval": "0%-20%",
                     "category": "no_Pfam",
                     "n_HOGs": len(sub_hogs) - n_pfam_ann,
                     "pct": (len(sub_hogs) - n_pfam_ann) / len(sub_hogs) * 100})

    # 0-20% KEGG / Pfam composition of the clinical-only and non-clinical-only
    # (environmental-only) HOG subsets
    excl_map = dict(zip(carrier["HOG_ID"], carrier["Exclusivity"]))
    for tag, excl_val in [("clinonly", "clinical-only"),
                          ("nonclinonly", "environmental-only")]:
        sub_h = set(h for h, e in excl_map.items() if e == excl_val)
        n_sub = len(sub_h)
        k_sub = kegg_cat[kegg_cat.index.isin(sub_h)]
        n_k_sub = len(k_sub)
        for name, n in k_sub.value_counts().items():
            fig_rows.append({"panel": f"kegg_{tag}_020", "interval": "0%-20%",
                             "category": name, "n_HOGs": int(n),
                             "pct": int(n) / n_k_sub * 100 if n_k_sub else 0.0})
        fig_rows.append({"panel": f"kegg_{tag}_020", "interval": "0%-20%",
                         "category": "no_KEGG", "n_HOGs": n_sub - n_k_sub,
                         "pct": (n_sub - n_k_sub) / n_sub * 100 if n_sub else 0.0})
        p_sub = pfam_classes[pfam_classes.index.isin(sub_h)]
        n_p_sub = len(p_sub)
        for name, n in p_sub.value_counts().items():
            fig_rows.append({"panel": f"pfam_{tag}_020", "interval": "0%-20%",
                             "category": name, "n_HOGs": int(n),
                             "pct": int(n) / n_p_sub * 100 if n_p_sub else 0.0})
        fig_rows.append({"panel": f"pfam_{tag}_020", "interval": "0%-20%",
                         "category": "no_Pfam", "n_HOGs": n_sub - n_p_sub,
                         "pct": (n_sub - n_p_sub) / n_sub * 100 if n_sub else 0.0})

    # 0-20% effector-HOG proportion of the clinical-only / non-clinical-only subsets
    # (HOG-level: HOG is an effector HOG if any of its member genes is an effector;
    #  class = both if it contains both cytoplasmic and apoplastic effectors)
    eff_genes = gvh[gvh["IsEffector"]]
    hog_cyto = set(eff_genes.loc[eff_genes["EffectorClass"].isin(["cytoplasmic", "both"]), "HOG_ID"])
    hog_apo = set(eff_genes.loc[eff_genes["EffectorClass"].isin(["apoplastic", "both"]), "HOG_ID"])
    for tag, excl_val in [("clinonly", "clinical-only"),
                          ("nonclinonly", "environmental-only")]:
        sub_h = set(h for h, e in excl_map.items() if e == excl_val)
        n_sub = len(sub_h)
        counts = {"non-effector": 0, "cytoplasmic": 0, "apoplastic": 0, "both": 0}
        for h in sub_h:
            if h not in hog_cyto and h not in hog_apo:
                counts["non-effector"] += 1
            elif h in hog_cyto and h in hog_apo:
                counts["both"] += 1
            elif h in hog_cyto:
                counts["cytoplasmic"] += 1
            else:
                counts["apoplastic"] += 1
        for cat, n in counts.items():
            fig_rows.append({"panel": f"effector_{tag}_020", "interval": "0%-20%",
                             "category": cat, "n_HOGs": int(n),
                             "pct": n / n_sub * 100 if n_sub else 0.0})

    # carrier pairs + per-isolate rate for 0-20% (quick-reference numbers)
    car020 = carrier[carrier["Interval"] == "0%-20%"]
    med_carrier_pairs = int(car020["n_clinical"].sum())
    env_carrier_pairs = int(car020["n_environmental"].sum())
    per_iso = car020[["n_clinical", "n_environmental"]].mean()
    fig_rows.append({"panel": "carrier_pairs", "interval": "0%-20%", "category": "clinical",
                     "n_HOGs": med_carrier_pairs,
                     "pct": med_carrier_pairs / (med_carrier_pairs + env_carrier_pairs) * 100})
    fig_rows.append({"panel": "carrier_pairs", "interval": "0%-20%", "category": "environmental",
                     "n_HOGs": env_carrier_pairs,
                     "pct": env_carrier_pairs / (med_carrier_pairs + env_carrier_pairs) * 100})
    fig_rows.append({"panel": "per_isolate_rate", "interval": "0%-20%", "category": "clinical",
                     "n_HOGs": per_iso["n_clinical"], "pct": per_iso["n_clinical"] / n_clin_iso})
    fig_rows.append({"panel": "per_isolate_rate", "interval": "0%-20%", "category": "environmental",
                     "n_HOGs": per_iso["n_environmental"], "pct": per_iso["n_environmental"] / n_env_iso})

    fig_df = pd.DataFrame(fig_rows)
    fig_out = os.path.join(args.output_dir, "figure_data.020.tsv")
    fig_df.to_csv(fig_out, sep="\t", index=False)
    print(f"figure data -> {fig_out}")

    # quick report
    print("\n=== exclusivity per interval ===")
    print(excl_iv.to_string())
    print("\n=== 0-20% carrier pairs ===")
    print(f"clinical={med_carrier_pairs}, environmental={env_carrier_pairs} "
          f"({med_carrier_pairs / (med_carrier_pairs + env_carrier_pairs) * 100:.1f}% clinical)")
    print("\n=== effector fraction by interval (gene level) ===")
    print(overview[["Interval", "n_genes", "n_effector_genes", "pct_effector_genes"]].to_string(index=False))
    print("\n=== 0-20% COG primary categories ===")
    print(cog_counts.to_string())


if __name__ == "__main__":
    main()
