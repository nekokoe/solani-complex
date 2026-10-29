#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assess the practical significance of the RF-selected clinical-associated HOGs
(Random-Forest feature importance + Wilcoxon test, per-species: F. solani and
F. falciformis) used in the pan-genome section of the manuscript.

Answers, with numbers ready for small main-figure panels:

1. Fraction of the pan-genome: n selected / total HOGs
2. Per-genome representation: how many selected HOGs a single genome carries,
   and what fraction of its gene content that corresponds to
   (split by clinical / non-clinical source, per training species and overall)
3. Effector content of selected HOGs vs. the pan-genome background
4. Kinase (protein kinase domain) content of selected HOGs vs. background
5. Most frequent functions (Pfam functional classes, KEGG BRITE level-2
   categories) with hypergeometric enrichment vs. background
6. Clinical vs. non-clinical preference of the selected HOGs
   (presence-rate difference per HOG, direction counts, signed-rank test)
7. Significance: how many selected HOGs survive Bonferroni correction

Outputs (--output-dir):
   rf_genes.hog_stats.tsv        per-HOG detail table
   suppl_table.rf_overview.tsv   headline numbers (fractions, tests)
   suppl_table.rf_functions.tsv  Pfam/KEGG composition + enrichment
   figure_data.rf.tsv            data for the small figure panels

The two candidate tables may come from either screening round; the column names
are mapped here so the rest of the script never has to care:

    count screen  classify_*/top_significant_features.txt
                  Feature, Importance, Positive_Rate, Negative_Rate,
                  P_value, Bonferroni_p_value
    PAV screen    pav_rerun/rf_pav_<species>_features.tsv   <-- current
                  HOG_ID, importance, p_pav, q_pav, clinical_rate,
                  nonclinical_rate, rate_delta ...

Bonferroni is recomputed when the table only carries the raw p (p x the number of
HOGs actually tested, taken from pav_rerun/rf_pav_summary.tsv).

Usage:
    python3 scripts/stats_rf_genes.py \
        --rf-solani 09.rf_clinical_genes/pav_rerun/rf_pav_solani_features.tsv \
        --rf-falciformis 09.rf_clinical_genes/pav_rerun/rf_pav_falciformis_features.tsv \
        --hog-count 05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv \
        --group-file bmu_ref_groups.211.txt \
        --effector-hog 08.pan_020_genes/effector_master.hog_level.tsv \
        --annotation-file 05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv \
        --interval-file 05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv \
        --kegg-category-file 08.pan_020_genes/kegg_pathway_categories.tsv \
        --output-dir 09.rf_clinical_genes
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, mannwhitneyu, wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_020_genes import classify_pfam  # noqa: E402


def parse_arguments():
    p = argparse.ArgumentParser(description="Analyse RF-selected clinical-associated HOGs")
    p.add_argument("--rf-solani", required=True)
    p.add_argument("--rf-falciformis", required=True)
    p.add_argument("--combined-list", default=None,
                   help="union list; derived from the two candidate tables when omitted")
    p.add_argument("--hog-count", required=True)
    p.add_argument("--group-file", required=True)
    p.add_argument("--effector-hog", required=True)
    p.add_argument("--effector-coverage", default=None,
                   help="effector_master.coverage.tsv (isolates lacking EffectorP output)")
    p.add_argument("--effector-gene", default=None,
                   help="effector_master.gene_level.tsv (per-gene EffectorP status)")
    p.add_argument("--gene-vs-hog", default=None, help="N0.gene_vs_HOG.tsv")
    p.add_argument("--gene-pfam", default=None, help="all.fusarium.gene_pfam.tsv (gene-level Pfam)")
    p.add_argument("--solani-cohort", required=True, help="all.solani.sample.list")
    p.add_argument("--falciformis-cohort", required=True, help="all.falciformis.sample.list")
    p.add_argument("--annotation-file", required=True)
    p.add_argument("--interval-file", required=True)
    p.add_argument("--kegg-category-file", default=None)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


CANDIDATE_ALIASES = {
    "Feature": ("HOG_ID",),
    "Importance": ("importance",),
    "Positive_Rate": ("clinical_rate",),
    "Negative_Rate": ("nonclinical_rate",),
    "P_value": ("p_pav", "p_value"),
    "Bonferroni_p_value": ("bonferroni_p_value", "p_bonferroni"),
}


def normalise_candidates(df, species, n_tested=None):
    """Accept either the count-screen or the PAV-screen candidate table.

    The two rounds name their columns differently; mapping them here keeps the
    downstream analysis identical so the regenerated supplementary tables are
    comparable with the ones they replace.
    """
    df = df.copy()
    for target, sources in CANDIDATE_ALIASES.items():
        if target in df.columns:
            continue
        for src in sources:
            if src in df.columns:
                df[target] = df[src]
                break
        else:
            if target == "Bonferroni_p_value" and n_tested and "P_value" in df.columns:
                df[target] = (df["P_value"] * n_tested).clip(upper=1.0)
            else:
                raise SystemExit(
                    f"{species}: cannot find a column for {target!r}; "
                    f"have {list(df.columns)}")
    return df


def tested_hog_counts(summary_path):
    """How many HOGs each species screened, needed for the Bonferroni column."""
    if not os.path.isfile(summary_path):
        return {}
    sm = pd.read_csv(summary_path, sep="\t")
    return dict(zip(sm["species"], sm["n_variable_hogs"]))


def bh_correct(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q_sorted = p[order] * n / np.arange(1, n + 1)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty(n)
    q[order] = q_sorted
    return np.minimum(q, 1.0)


def hypergeom_p(k, K, n, M):
    """P(X >= k) drawing n from M with K marked."""
    return hypergeom.sf(k - 1, M, K, n)


def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    # ---------- selected HOG sets ----------
    n_tested = tested_hog_counts(os.path.join(
        os.path.dirname(os.path.abspath(args.rf_solani)), "rf_pav_summary.tsv"))
    sol = normalise_candidates(pd.read_csv(args.rf_solani, sep="\t"), "solani",
                               n_tested.get("solani"))
    fal = normalise_candidates(pd.read_csv(args.rf_falciformis, sep="\t"), "falciformis",
                               n_tested.get("falciformis"))
    sol_set, fal_set = set(sol["Feature"]), set(fal["Feature"])
    if args.combined_list and os.path.isfile(args.combined_list):
        with open(args.combined_list) as fh:
            comb = [l.strip() for l in fh if l.strip()]
    else:
        comb = sorted(sol_set | fal_set)
    comb_set = set(comb)
    n_total_hogs = 29080
    print(f"selected: solani {len(sol_set)}, falciformis {len(fal_set)}, "
          f"combined {len(comb_set)} (overlap {len(sol_set & fal_set)})")

    # ---------- annotation / effector / intervals ----------
    ann = pd.read_csv(args.annotation_file, sep="\t", usecols=[
        "HOG_ID", "Pfam_Domain", "KEGG_Pathway", "NOG_Description"])
    ann = ann.set_index("HOG_ID")
    ann = ann[~ann.index.duplicated(keep="first")]
    iv = pd.read_csv(args.interval_file, sep="\t")
    interval_map = dict(zip(iv["HOG"], iv["Presence_Rate_Interval"]))
    eff = pd.read_csv(args.effector_hog, sep="\t",
                      usecols=["HOG_ID", "EffectorHOG", "n_cytoplasmic", "n_apoplastic"])
    eff_set = set(eff.loc[eff["EffectorHOG"], "HOG_ID"])
    eff_cyto = dict(zip(eff["HOG_ID"], eff["n_cytoplasmic"]))
    eff_apo = dict(zip(eff["HOG_ID"], eff["n_apoplastic"]))
    all_hogs = set(interval_map.keys())

    def effector_class(h):
        cy, ap = eff_cyto.get(h, 0), eff_apo.get(h, 0)
        if cy > 0 and ap > 0:
            return "both"
        if cy > 0:
            return "cytoplasmic"
        if ap > 0:
            return "apoplastic"
        return "non-effector"

    # kinase: protein kinase domains
    kin = ann["Pfam_Domain"].dropna().str.contains("pkinase|hiska", case=False, regex=True)
    kin_set = set(kin[kin].index)

    # Pfam functional class per HOG
    pfam_class = ann["Pfam_Domain"].dropna().apply(
        lambda v: classify_pfam(str(v).split(";")))
    pfam_class = pfam_class[pfam_class.index.isin(all_hogs)]
    pfam_all = set(pfam_class.index)
    # actual (first/best-scoring) Pfam domain per HOG, for domain-level pies
    pfam_domain = ann["Pfam_Domain"].dropna().apply(
        lambda v: str(v).split(";")[0].strip())
    pfam_domain = pfam_domain[pfam_domain.index.isin(all_hogs)]
    pfam_domain = pfam_domain[pfam_domain != "-"]

    # KEGG level-2 category per HOG (first pathway)
    cat_map = {}
    if args.kegg_category_file and os.path.isfile(args.kegg_category_file):
        cat_df = pd.read_csv(args.kegg_category_file, sep="\t", dtype={"pathway_id": str})
        cat_map = {pid.zfill(5): l2 for pid, l2 in zip(cat_df["pathway_id"], cat_df["level2"])}
    def level2_of(pid):
        v = cat_map.get(str(pid).zfill(5))
        if isinstance(v, str) and " " in v:
            return v.split(" ", 1)[1]
        return v or "Unmapped"
    kegg_ann = ann["KEGG_Pathway"].dropna()
    kegg_ann = kegg_ann[kegg_ann != "-"]
    kegg_cat = kegg_ann.apply(lambda v: level2_of(str(v).split(",")[0].strip()))
    kegg_cat = kegg_cat[kegg_cat.index.isin(all_hogs)]
    kegg_all = set(kegg_cat.index)

    # ---------- per-HOG detail table ----------
    rows = []
    for hog in comb:
        in_sol = hog in sol_set
        in_fal = hog in fal_set
        s = sol[sol["Feature"] == hog]
        f = fal[fal["Feature"] == hog]
        imp = max(float(s["Importance"].iloc[0]) if len(s) else 0,
                  float(f["Importance"].iloc[0]) if len(f) else 0)
        pos_r = float(s["Positive_Rate"].iloc[0]) if len(s) else (
            float(f["Positive_Rate"].iloc[0]) if len(f) else np.nan)
        neg_r = float(s["Negative_Rate"].iloc[0]) if len(s) else (
            float(f["Negative_Rate"].iloc[0]) if len(f) else np.nan)
        pv = float(s["P_value"].iloc[0]) if len(s) else (
            float(f["P_value"].iloc[0]) if len(f) else np.nan)
        bpv = float(s["Bonferroni_p_value"].iloc[0]) if len(s) else (
            float(f["Bonferroni_p_value"].iloc[0]) if len(f) else np.nan)
        rows.append({
            "HOG_ID": hog,
            "Species": "both" if (in_sol and in_fal) else ("solani" if in_sol else "falciformis"),
            "Importance": imp,
            "P_value": pv,
            "Bonferroni_p": bpv,
            "Clinical_rate": pos_r,
            "NonClinical_rate": neg_r,
            "Rate_delta": pos_r - neg_r,
            "IsEffector": hog in eff_set,
            "EffectorClass": effector_class(hog),
            "IsKinase": hog in kin_set,
            "PfamClass": pfam_class.get(hog, "No Pfam"),
            "Pfam_Domain": pfam_domain.get(hog, "No Pfam"),
            "KEGG_Category": kegg_cat.get(hog, "No KEGG"),
            "Presence_Interval": interval_map.get(hog, "NA"),
        })
    detail = pd.DataFrame(rows).set_index("HOG_ID")
    detail_out = os.path.join(args.output_dir, "rf_genes.hog_stats.tsv")
    detail.to_csv(detail_out, sep="\t")
    print(f"per-HOG detail -> {detail_out}")

    # ---------- overview numbers ----------
    ov = {}
    ov["n_selected"] = len(comb_set)
    ov["n_total_HOGs"] = n_total_hogs
    ov["pct_of_HOGs"] = len(comb_set) / n_total_hogs * 100
    ov["n_solani"], ov["n_falciformis"] = len(sol_set), len(fal_set)
    ov["n_effector_selected"] = int((detail["IsEffector"]).sum())
    ov["pct_effector_selected"] = detail["IsEffector"].mean() * 100
    ov["n_effector_background"] = len(eff_set & all_hogs)
    ov["pct_effector_background"] = len(eff_set & all_hogs) / n_total_hogs * 100
    ov["effector_hypergeom_p"] = hypergeom_p(
        ov["n_effector_selected"], ov["n_effector_background"],
        len(comb_set), n_total_hogs)
    ov["n_kinase_selected"] = int(detail["IsKinase"].sum())
    ov["pct_kinase_selected"] = detail["IsKinase"].mean() * 100
    ov["n_kinase_background"] = len(kin_set & all_hogs)
    ov["pct_kinase_background"] = len(kin_set & all_hogs) / n_total_hogs * 100
    ov["kinase_hypergeom_p"] = hypergeom_p(
        ov["n_kinase_selected"], ov["n_kinase_background"], len(comb_set), n_total_hogs)
    ov["n_bonferroni_solani"] = int((sol["Bonferroni_p_value"] < 0.05).sum())
    ov["n_bonferroni_falciformis"] = int((fal["Bonferroni_p_value"] < 0.05).sum())
    ov["n_bonferroni_combined"] = int((detail["Bonferroni_p"] < 0.05).sum())

    # ---------- per-genome representation ----------
    hc = pd.read_csv(args.hog_count, sep="\t")
    hc = hc[pd.to_numeric(hc["GeneCount"], errors="coerce").fillna(0) > 0]
    groups = pd.read_csv(args.group_file, sep="\t")
    groups.columns = [c.strip() for c in groups.columns]
    groups = groups.rename(columns={groups.columns[0]: "Sample"})
    med_env = dict(zip(groups["Sample"], groups["MedEnv"]))
    species_map = dict(zip(groups["Sample"], groups["Species"])) if "Species" in groups.columns else {}

    hc["IsSel"] = hc["HOG_ID"].isin(comb_set)
    per_samp = hc.groupby("SampleID").size().rename("n_hogs_total").reset_index()
    sel_samp = hc[hc["IsSel"]].groupby("SampleID").size().rename("n_hogs_selected").reset_index()
    gene_sums = hc.groupby("SampleID")["GeneCount"].sum().rename("n_genes").reset_index()
    gdf = per_samp.merge(sel_samp, on="SampleID", how="left").fillna(0)
    gdf = gdf.merge(gene_sums, on="SampleID")
    gdf["sel_frac_hogs"] = gdf["n_hogs_selected"] / gdf["n_hogs_total"] * 100
    gdf["sel_frac_genes"] = gdf["n_hogs_selected"] / gdf["n_genes"] * 100
    gdf["Source"] = gdf["SampleID"].map(med_env).fillna("unknown")
    gdf["Species"] = gdf["SampleID"].map(species_map).fillna("")
    ov["genome_mean_selected"] = gdf["n_hogs_selected"].mean()
    ov["genome_mean_frac_genes"] = gdf["sel_frac_genes"].mean()
    # clinical vs non-clinical within each training species
    genome_sp = {}
    for sp_name in ["solani", "falciformis"]:
        sub = gdf[gdf["Species"] == sp_name]
        med = sub.loc[sub["Source"] == "med", "n_hogs_selected"].values
        env = sub.loc[sub["Source"] == "env", "n_hogs_selected"].values
        p = mannwhitneyu(med, env).pvalue if (len(med) >= 8 and len(env) >= 8) else np.nan
        genome_sp[sp_name] = {
            "median_clinical": float(np.median(med)) if len(med) else np.nan,
            "median_nonclinical": float(np.median(env)) if len(env) else np.nan,
            "n_clinical": len(med), "n_nonclinical": len(env),
            "clin_vs_nonclin_p": p,
        }
    ov["genome_sp"] = genome_sp
    # per-genome detail (for boxplots)
    pg = gdf[["SampleID", "Source", "Species", "n_hogs_total", "n_hogs_selected",
              "sel_frac_hogs", "sel_frac_genes"]]

    # ---------- per-genome effector / kinase fractions (GENE level, paired within genome) ----------
    # For each isolate of the two training species: fraction of effector (or
    # kinase) GENES among the genes belonging to the SELECTED HOGs present in
    # that genome vs among ALL genes of that genome. Effector status comes from
    # the isolate's own EffectorP output (genes without a prediction = non-
    # effector); kinase status from gene-level Pfam (Pkinase/hiska domains).
    # Isolates without EffectorP output are excluded from the effector
    # comparison.
    with open(args.solani_cohort) as fh:
        sol_cohort = set(l.strip() for l in fh if l.strip())
    with open(args.falciformis_cohort) as fh:
        fal_cohort = set(l.strip() for l in fh if l.strip())
    no_effector = set()
    if args.effector_coverage and os.path.isfile(args.effector_coverage):
        cov_df = pd.read_csv(args.effector_coverage, sep="\t")
        no_effector = set(cov_df.loc[~cov_df["HasEffectorP"], "Sample"])

    if not (args.effector_gene and args.gene_vs_hog):
        raise SystemExit("gene-level effector/kinase comparison requires "
                         "--effector-gene and --gene-vs-hog")
    gvh = pd.read_csv(args.gene_vs_hog, sep="\t", usecols=["GeneID", "SampleID", "HOG_ID"])
    gvh = gvh.drop_duplicates(subset=["GeneID"])
    all_genes = gvh.groupby("SampleID")["GeneID"].apply(set)
    sel_genes = {}
    for sp_name, cohort, sel in [("solani", sol_cohort, sol_set),
                                 ("falciformis", fal_cohort, fal_set)]:
        sub = gvh[gvh["SampleID"].isin(cohort) & gvh["HOG_ID"].isin(sel)]
        sel_genes.update(sub.groupby("SampleID")["GeneID"].apply(set))

    eff_gene = pd.read_csv(args.effector_gene, sep="\t", usecols=["GeneID", "IsEffector"])
    eff_gene_set = set(eff_gene.loc[eff_gene["IsEffector"], "GeneID"])
    kin_gene_set = set()
    if args.gene_pfam and os.path.isfile(args.gene_pfam):
        pfam_g = pd.read_csv(args.gene_pfam, sep="\t", header=None, usecols=[0, 2],
                             names=["GeneID", "Dom"])
        kin_gene_set = set(pfam_g.loc[pfam_g["Dom"].str.contains(
            "pkinase|hiska", case=False, regex=True), "GeneID"])

    eff_rows, kin_rows = [], []
    for sample in sorted(sol_cohort | fal_cohort):
        if sample not in all_genes:
            continue
        genes_all = all_genes[sample]
        genes_sel = sel_genes.get(sample, set())
        if not genes_sel:
            continue
        sp_name = species_map.get(sample, "")
        kin_rows.append((sample, sp_name,
                         len(genes_sel & kin_gene_set) / len(genes_sel),
                         len(genes_all & kin_gene_set) / len(genes_all)))
        if sample not in no_effector:
            eff_rows.append((sample, sp_name,
                             len(genes_sel & eff_gene_set) / len(genes_sel),
                             len(genes_all & eff_gene_set) / len(genes_all)))

    pg_eff = pd.DataFrame(eff_rows, columns=["SampleID", "Species", "SelectedEffectorFrac",
                                             "AllEffectorFrac"])
    pg_kin = pd.DataFrame(kin_rows, columns=["SampleID", "Species", "SelectedKinaseFrac",
                                             "AllKinaseFrac"])
    pg = pg.merge(pg_eff, on="SampleID", how="left").merge(pg_kin, on="SampleID", how="left")
    ov["n_genomes_effector"] = len(pg_eff)

    def paired_by_species(df, sel_col, all_col, label):
        out = {}
        for sp_name in ["solani", "falciformis"]:
            sub = df[df["Species"] == sp_name]
            if len(sub) < 8:
                out[sp_name] = {"mean_selected": np.nan, "mean_all": np.nan, "p": np.nan, "n": len(sub)}
                continue
            p = wilcoxon(sub[sel_col], sub[all_col]).pvalue
            out[sp_name] = {"mean_selected": float(sub[sel_col].mean()),
                            "mean_all": float(sub[all_col].mean()),
                            "p": float(p), "n": len(sub)}
        ov[label] = out

    paired_by_species(pg_eff, "SelectedEffectorFrac", "AllEffectorFrac", "effector_genome")
    paired_by_species(pg_kin, "SelectedKinaseFrac", "AllKinaseFrac", "kinase_genome")

    pg_out = os.path.join(args.output_dir, "rf_genes.per_genome.tsv")
    pg.to_csv(pg_out, sep="\t", index=False)
    print(f"per-genome -> {pg_out}")

    pg_out = os.path.join(args.output_dir, "rf_genes.per_genome.tsv")
    pg.to_csv(pg_out, sep="\t", index=False)
    print(f"per-genome -> {pg_out}")

    # ---------- preference at the HOG level (per species) ----------
    def pref_stats(sp_df, sp_name):
        d = sp_df["Positive_Rate"] - sp_df["Negative_Rate"]
        clin_b = int((d > 0.1).sum())
        nonclin_b = int((d < -0.1).sum())
        pair = sp_df[["Positive_Rate", "Negative_Rate"]].dropna()
        p = np.nan
        if len(pair) >= 8:
            try:
                p = float(wilcoxon(pair["Positive_Rate"], pair["Negative_Rate"]).pvalue)
            except ValueError:
                p = np.nan
        return {
            "species": sp_name, "n": len(sp_df),
            "clinical_biased": clin_b, "nonclinical_biased": nonclin_b,
            "balanced": len(d) - clin_b - nonclin_b,
            "mean_clinical_rate": float(pair["Positive_Rate"].mean()),
            "mean_nonclinical_rate": float(pair["Negative_Rate"].mean()),
            "signedrank_p": p,
            "n_bonferroni": int((sp_df["Bonferroni_p_value"] < 0.05).sum()),
        }

    pref = {sp_name: pref_stats(df, sp_name)
            for sp_name, df in [("solani", sol), ("falciformis", fal)]}

    ov_df = pd.DataFrame([{"Metric": k, "Value": v} for k, v in ov.items()])
    ov_out = os.path.join(args.output_dir, "suppl_table.rf_overview.tsv")
    ov_df.to_csv(ov_out, sep="\t", index=False)
    print(f"overview -> {ov_out}")

    # ---------- functional composition + enrichment (per species) ----------
    fun_rows = []
    for sp_name in ["solani", "falciformis"]:
        sel_detail = detail[detail["Species"] == sp_name]
        n_sel_sp = len(sel_detail)
        if n_sel_sp == 0:
            continue
        for cls in ["Transporter", "Transcription factor", "Oxidoreductase", "HET domain",
                    "Ankyrin repeat", "Protein kinase", "Peptidase", "Methyltransferase",
                    "WD40/TPR repeat", "Other Pfam"]:
            bg_set = set(pfam_class[pfam_class == cls].index)
            k = int((sel_detail["PfamClass"] == cls).sum())
            if k == 0 and len(bg_set) < 5:
                continue
            fun_rows.append({
                "Species": sp_name, "Source": "Pfam", "Category": cls,
                "n_selected": k,
                "pct_of_selected": k / n_sel_sp * 100,
                "n_background": len(bg_set),
                "pct_of_background": len(bg_set) / n_total_hogs * 100,
                "hypergeom_p": hypergeom_p(k, len(bg_set), n_sel_sp, n_total_hogs),
            })
        for cat, n_bg in kegg_cat.value_counts().items():
            k = int((sel_detail["KEGG_Category"] == cat).sum())
            if k == 0:
                continue
            fun_rows.append({
                "Species": sp_name, "Source": "KEGG", "Category": cat,
                "n_selected": k,
                "pct_of_selected": k / n_sel_sp * 100,
                "n_background": n_bg,
                "pct_of_background": n_bg / n_total_hogs * 100,
                "hypergeom_p": hypergeom_p(k, n_bg, n_sel_sp, n_total_hogs),
            })
    fun_df = pd.DataFrame(fun_rows)
    if len(fun_df):
        fun_df["padj_BH"] = bh_correct(fun_df["hypergeom_p"])
        fun_df["log10_padj"] = -np.log10(fun_df["padj_BH"].clip(lower=1e-300))
        fun_df = fun_df.sort_values(["Species", "padj_BH"])
    fun_out = os.path.join(args.output_dir, "suppl_table.rf_functions.tsv")
    fun_df.to_csv(fun_out, sep="\t", index=False)
    print(f"functions -> {fun_out}")

    # ---------- functional-composition preference by category (LEfSe-style, GENE level) ----------
    # Within the genes of each genome that belong to the set in scope (selected
    # HOGs for the RF set; 0-20% accessory HOGs for the pan-genome set), compute
    # the proportion of genes in each functional category, per genome. Compare
    # clinical vs non-clinical genomes (Mann-Whitney U). Effect = mean
    # (clinical genome fraction - non-clinical genome fraction), in %.
    pan020_hogs = set(h for h, v in interval_map.items() if v == "0%-20%")
    pan020_genes = gvh[gvh["HOG_ID"].isin(pan020_hogs)].groupby("SampleID")["GeneID"].apply(set)

    # gene-level functional categories
    pfam_by_gene = pfam_g.groupby("GeneID")["Dom"].apply(
        lambda ds: classify_pfam(list(ds)))
    gene_pfam_cat = pfam_by_gene.to_dict()
    gene2hog = dict(zip(gvh["GeneID"], gvh["HOG_ID"]))
    hog2kegg = kegg_cat.to_dict()
    gene_kegg_cat = {g: hog2kegg.get(h, "No KEGG") for g, h in gene2hog.items()}

    n_clin_iso = int(sum(1 for v in med_env.values() if v == "med"))
    n_env_iso = int(sum(1 for v in med_env.values() if v == "env"))

    def pref_cat_for_set(set_name, genes_by_genome, genomes):
        from collections import Counter
        rows = []
        per_iso_rows = []
        for src, gene_cat, default_cat in [("Pfam", gene_pfam_cat, "No Pfam"),
                                           ("KEGG", gene_kegg_cat, "No KEGG")]:
            # fractions are computed among ANNOTATED in-scope genes (the
            # default 'No X' category is excluded from the denominator) so that
            # annotation-coverage differences do not drive the comparison
            cat_counts = {}                       # cat -> {'med': [(sample, frac)], 'env': [...]}
            no_annot = {"med": [], "env": []}     # per-genome unannotated share
            for g in genomes:
                genes = genes_by_genome.get(g)
                if not genes:
                    continue
                src_g = med_env.get(g)
                if src_g not in ("med", "env"):
                    continue
                n_scope = len(genes)
                cnt = Counter(gene_cat.get(gene, default_cat) for gene in genes)
                n_annot = n_scope - cnt.get(default_cat, 0)
                no_annot[src_g].append(cnt.get(default_cat, 0) / n_scope)
                if n_annot == 0:
                    continue
                for cat, c in cnt.items():
                    if cat == default_cat:
                        continue
                    cat_counts.setdefault(cat, {"med": [], "env": []})[src_g].append((g, c / n_annot))
            for cat, grp in cat_counts.items():
                # 'Other Pfam' (Pfam catch-all) and 'Unmapped' (KEGG pathway
                # not in the BRITE map) are not meaningful functional
                # categories and are excluded from the statistical comparison;
                # their genes still count towards the annotated denominator.
                if cat in ("Other Pfam", "Unmapped"):
                    continue
                cf = [f for _, f in grp["med"]]
                nf = [f for _, f in grp["env"]]
                if len(cf) < 5 or len(nf) < 5:
                    continue
                effect = (float(np.mean(cf)) - float(np.mean(nf))) * 100
                p = np.nan
                try:
                    p = float(mannwhitneyu(cf, nf).pvalue)
                except ValueError:
                    p = np.nan
                for sample, frac in grp["med"]:
                    per_iso_rows.append({"Set": set_name, "Source": src,
                                         "Category": cat, "SampleID": sample,
                                         "Group": "med", "Fraction": frac})
                for sample, frac in grp["env"]:
                    per_iso_rows.append({"Set": set_name, "Source": src,
                                         "Category": cat, "SampleID": sample,
                                         "Group": "env", "Fraction": frac})
                rows.append({
                    "Set": set_name, "Source": src, "Category": cat,
                    "n_Genes": sum(len(genes_by_genome.get(g, [])) for g in genomes),
                    "Effect_pp": effect,
                    "Mean_Clinical_pct": float(np.mean(cf)) * 100,
                    "Mean_NonClinical_pct": float(np.mean(nf)) * 100,
                    "Direction": "clinical" if effect > 0 else "non-clinical",
                    "pvalue": p,
                    "n_clin_genomes": len(cf),
                    "n_nonclin_genomes": len(nf),
                })
            if no_annot["med"] and no_annot["env"] and len(no_annot["med"]) >= 5 and len(no_annot["env"]) >= 5:
                p = np.nan
                try:
                    p = float(mannwhitneyu(no_annot["med"], no_annot["env"]).pvalue)
                except ValueError:
                    p = np.nan
                rows.append({
                    "Set": set_name, "Source": src, "Category": default_cat,
                    "n_Genes": sum(len(genes_by_genome.get(g, [])) for g in genomes),
                    "Effect_pp": (float(np.mean(no_annot["med"])) - float(np.mean(no_annot["env"]))) * 100,
                    "Direction": "clinical" if (np.mean(no_annot["med"]) - np.mean(no_annot["env"])) > 0 else "non-clinical",
                    "pvalue": p,
                    "n_clin_genomes": len(no_annot["med"]),
                    "n_nonclin_genomes": len(no_annot["env"]),
                })
        return pd.DataFrame(rows), pd.DataFrame(per_iso_rows)

    pref_rf, per_iso_rf = pref_cat_for_set("rf_selected", sel_genes,
                                           sorted(sol_cohort | fal_cohort))
    pref_pan, per_iso_pan = pref_cat_for_set("pan_0_20", pan020_genes,
                                             sorted(all_genes.keys()))
    pref_cat = pd.concat([pref_rf, pref_pan], ignore_index=True)
    pref_cat_out = os.path.join(args.output_dir, "suppl_table.rf_clinical_preference.tsv")
    pref_cat.to_csv(pref_cat_out, sep="\t", index=False)
    print(f"functional-composition preference (gene level) -> {pref_cat_out}")
    print(pref_cat.sort_values("Effect_pp").to_string(index=False))

    # per-isolate fractions (for pairwise-difference distributions)
    per_iso = pd.concat([per_iso_rf, per_iso_pan], ignore_index=True)
    per_iso_out = os.path.join(args.output_dir, "suppl_table.rf_preference_per_isolate.tsv")
    per_iso.to_csv(per_iso_out, sep="\t", index=False)
    print(f"per-isolate fractions -> {per_iso_out}")

    # ---------- figure data (species-aware) ----------
    fig_rows = []

    def frow(panel, species, category, value, extra=""):
        fig_rows.append({"panel": panel, "species": species, "category": category,
                         "value": value, "n": extra})

    # fraction of pan-genome (per species + combined)
    for sp_name, n_sp in [("solani", len(sol_set)), ("falciformis", len(fal_set))]:
        frow("fraction", sp_name, "selected", n_sp / n_total_hogs * 100, n_sp)
    frow("fraction", "combined", "selected", len(comb_set) / n_total_hogs * 100, len(comb_set))
    frow("fraction", "combined", "total", 100.0, n_total_hogs)

    # per-genome representation means (plot reads raw values from per_genome.tsv)
    for sp_name in ["solani", "falciformis"]:
        sub = gdf[gdf["Species"] == sp_name]
        frow("per_genome", sp_name, "mean", float(sub["n_hogs_selected"].mean()), len(sub))
        frow("per_genome", sp_name, "median_clinical", genome_sp[sp_name]["median_clinical"], genome_sp[sp_name]["n_clinical"])
        frow("per_genome", sp_name, "median_nonclinical", genome_sp[sp_name]["median_nonclinical"], genome_sp[sp_name]["n_nonclinical"])
        frow("per_genome", sp_name, "clin_vs_nonclin_p", genome_sp[sp_name]["clin_vs_nonclin_p"], None)
    frow("per_genome", "combined", "frac_genes", ov["genome_mean_frac_genes"], None)

    # per-genome paired effector / kinase (plot reads raw values; p here)
    for sp_name in ["solani", "falciformis"]:
        e = ov["effector_genome"][sp_name]
        k = ov["kinase_genome"][sp_name]
        frow("effector_paired", sp_name, "mean_selected", e["mean_selected"], e["n"])
        frow("effector_paired", sp_name, "mean_all", e["mean_all"], e["n"])
        frow("effector_paired", sp_name, "p", e["p"], None)
        frow("kinase_paired", sp_name, "mean_selected", k["mean_selected"], k["n"])
        frow("kinase_paired", sp_name, "mean_all", k["mean_all"], k["n"])
        frow("kinase_paired", sp_name, "p", k["p"], None)

    # clinical-preference by functional category (LEfSe-style, replaces old
    # Pfam/KO proportion panels): one row per (set, category)
    panel_of_set = {"rf_selected": "pref_rf", "pan_0_20": "pref_pan"}
    for _, r in pref_cat.iterrows():
        frow(panel_of_set[r["Set"]], "combined", f'{r["Source"]}|{r["Category"]}',
             r["Effect_pp"], r["n_Genes"])
        # store p-value under a companion category key
        frow(panel_of_set[r["Set"]], "combined", f'{r["Source"]}|{r["Category"]}|p',
             r["pvalue"], r["n_Genes"])

    # preference per species
    for sp_name in ["solani", "falciformis"]:
        p = pref[sp_name]
        frow("preference", sp_name, "mean_clinical_rate", p["mean_clinical_rate"], p["n"])
        frow("preference", sp_name, "mean_nonclinical_rate", p["mean_nonclinical_rate"], p["n"])
        frow("preference", sp_name, "clinical_biased", p["clinical_biased"], None)
        frow("preference", sp_name, "nonclinical_biased", p["nonclinical_biased"], None)
        frow("preference", sp_name, "balanced", p["balanced"], None)
        frow("preference", sp_name, "signedrank_p", p["signedrank_p"], None)

    # significance per species
    for sp_name in ["solani", "falciformis"]:
        frow("significance", sp_name, "raw_p005", pref[sp_name]["n"], None)
        frow("significance", sp_name, "bonferroni_p005", pref[sp_name]["n_bonferroni"], None)

    fig_df = pd.DataFrame(fig_rows)
    fig_out = os.path.join(args.output_dir, "figure_data.rf.tsv")
    fig_df.to_csv(fig_out, sep="\t", index=False)
    print(f"figure data -> {fig_out}")

    # ---------- console report ----------
    print("\n=== headline (per species) ===")
    for sp_name in ["solani", "falciformis"]:
        p = pref[sp_name]
        g = genome_sp[sp_name]
        e = ov["effector_genome"][sp_name]
        k = ov["kinase_genome"][sp_name]
        print(f"  [{sp_name}] selected={p['n']} ({p['n']/n_total_hogs*100:.2f}%) | "
              f"per-genome median clin {g['median_clinical']:.0f} vs nonclin {g['median_nonclinical']:.0f} "
              f"(P={g['clin_vs_nonclin_p']:.2g}) | "
              f"effector per-genome {e['mean_selected']*100:.1f}% vs {e['mean_all']*100:.1f}% (P={e['p']:.2g}) | "
              f"kinase per-genome {k['mean_selected']*100:.1f}% vs {k['mean_all']*100:.1f}% (P={k['p']:.2g}) | "
              f"preference {p['clinical_biased']}/{p['nonclinical_biased']}/{p['balanced']} (P={p['signedrank_p']:.2g}) | "
              f"Bonferroni {p['n_bonferroni']}/{p['n']}")
    print("\n=== functional enrichment (top 10 by padj) ===")
    if len(fun_df):
        print(fun_df.head(10)[["Species", "Source", "Category", "n_selected",
                               "padj_BH"]].to_string(index=False))


if __name__ == "__main__":
    main()
