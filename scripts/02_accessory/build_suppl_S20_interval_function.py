#!/usr/bin/env python3
"""Supplementary Table S20 - functional composition of the pan-genome by presence interval.

The Results subsection on the pan-genome characterises only the 0-20% interval, although that
interval is a minority of the 29,080 orthogroups (39%) while the near-core 80-100% interval is
larger (44%). This table gives every interval the same treatment, for the complex as a whole and
for the two species that carry enough isolates to be analysed separately.

Interval definition: `05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv`
(presence over the 211 FSSC genomes) and, for the per-species sheets,
`05.stats_pan_by_species_shorter/per_species/<sp>.HOG.interval.tsv` (presence over that species'
own genomes). Collapsed groups 0-20 / 20-80 / 80-100 sit alongside the five individual bins.

Counting unit: one orthogroup, one vote. Each HOG is labelled by the majority class among its own
annotated member genes, and the intervals are then compared over HOGs. This is the only unit that
keeps the two axes comparable:

  * counting gene copies instead weights each interval by its copy number, and the near-core
    interval holds 80% of all gene copies, so its composition converges on the genome-wide
    average and every contrast is compressed;
  * labelling a HOG from the union of its members' domains lets the label track HOG size, so a
    large core orthogroup accumulates hundreds of domains and a 'first matching keyword' rule
    then reports 91.6% of core orthogroups as transporters, against 4.8% when each member gene
    is counted once.

Votes are cast only by genes that carry an annotation, and ties are broken alphabetically.

Writes 00.manuscript/supplementary/suppl_table.S20_pan_genome_by_interval.xlsx
  overview              per interval: HOGs, genes, effector and annotation-coverage rates
  Pfam_MF               GO molecular-function / EC class, complex-wide intervals
  KEGG_level2           KEGG level-2 category, complex-wide intervals
  Pfam_MF_species       the same axis on each species' own intervals
  KEGG_level2_species   the same axis on each species' own intervals
  candidates_intervals  where the random-forest candidates sit on the presence gradient
  candidates_function   candidate composition against the complex and the within-species baseline
  marker_search         genome-wide test of whether any orthogroup marks a clinical isolate
  marker_candidates_80_20  the orthogroups meeting the relaxed 80% clinical / 20% non-clinical
                        criterion, as a shortlist for experimental follow-up
  resistance_genes      antifungal-resistance orthogroups on the presence gradient
  resistance_source_assoc  the resistance orthogroups whose prevalence differs between clinical
                        and non-clinical isolates after BH correction
  definitions           interval definition and how the labels were assigned
"""

import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

IV = os.path.join(ROOT, "05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv")
IV_SP = os.path.join(ROOT, "05.stats_pan_by_species_shorter/per_species")
GENE_HOG = os.path.join(ROOT, "05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv")
GENE_PFAM = os.path.join(ROOT, "02.annotation/all_isolates/all.fusarium.gene_pfam.tsv")
GENE_EGGNOG = os.path.join(ROOT, "02.annotation/all_isolates/all.fusarium.gene_eggnog.tsv")
PFAM_LOOKUP = os.path.join(ROOT, "02.annotation/pfam_function_lookup.tsv")
KEGGCAT = os.path.join(ROOT, "08.pan_020_genes/kegg_pathway_categories.tsv")
RES_IV = os.path.join(ROOT, "08.pan_020_genes/resistance_genes_by_interval.tsv")
RES_TEST = os.path.join(ROOT, "08.pan_020_genes/resistance_genes_clinical_vs_env.tsv")
MARKER = os.path.join(ROOT, "08.pan_020_genes/marker_gene_search.tsv")
MARKER80 = os.path.join(ROOT, "08.pan_020_genes/marker_gene_search_80_20.tsv")
EFF_IV = os.path.join(ROOT, "08.pan_020_genes/effector_genes_by_interval.tsv")
EFF_TEST = os.path.join(ROOT, "08.pan_020_genes/effector_genes_per_genome.tsv")
OVERVIEW = os.path.join(ROOT, "08.pan_020_genes/suppl_table.020_overview_by_interval.tsv")
CAND = os.path.join(ROOT, "09.rf_clinical_genes/pav_rerun")
OUT = os.path.join(ROOT, "00.manuscript/supplementary",
                   "suppl_table.S20_pan_genome_by_interval.xlsx")
OUT_CAND = os.path.join(ROOT, "00.manuscript/supplementary",
                        "suppl_table.S21_rf_candidates.xlsx")
OUT_RES = os.path.join(ROOT, "00.manuscript/supplementary",
                       "suppl_table.S22_resistance.xlsx")

BINS = ["0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%"]
COLS = BINS + ["0-20", "20-80", "80-100"]
MID = ["20%-40%", "40%-60%", "60%-80%"]
SPECIES = ["solani", "falciformis"]
AXES = (("go_mf_class", "Pfam_MF"), ("KEGG_level2", "KEGG_level2"))


def collapse(b):
    return "0-20" if b == "0%-20%" else ("80-100" if b == "80%-100%" else "20-80")


def hog_labels(gene_sets, genes_of_hog):
    """majority class over the annotated member genes of one HOG; ties break alphabetically"""
    votes = Counter()
    for g in genes_of_hog:
        s = gene_sets.get(g)
        if s:
            votes.update(s)
    if not votes:
        return None
    top = max(votes.values())
    return sorted(k for k, v in votes.items() if v == top)[0]


def summarise(labels, iv, gap_label, key):
    d = iv.set_index("HOG").join(labels)
    d[key] = d[key].fillna(gap_label)
    t = pd.concat([pd.crosstab(d[key], d.bin), pd.crosstab(d[key], d.grp3)], axis=1)
    t = t.reindex(columns=COLS, fill_value=0).astype(int)
    n_all = t.sum()
    n_ann = n_all - t.loc[gap_label]
    p_all = (t / n_all * 100).round(2)
    p_ann = (t / n_ann * 100).round(2)
    p_ann.loc[gap_label] = np.nan
    core, rare = int(n_ann["80-100"]), int(n_ann["0-20"])
    rows = []
    for k in t.index:
        a, b = int(t.at[k, "80-100"]), int(t.at[k, "0-20"])
        if a + b == 0:
            continue
        try:
            chi = stats.chi2_contingency([[a, core - a], [b, rare - b]]).pvalue
        except ValueError:
            chi = float("nan")
        row = {key: k}
        row.update({c: p_all.at[k, c] for c in COLS})
        row.update({"pct_annotated_%s" % c: p_ann.at[k, c] for c in COLS})
        row["n_HOGs_annotated_0-20"] = b
        row["n_HOGs_annotated_80-100"] = a
        row["fold_core_vs_rare"] = round((a / core) / (b / rare), 2) if b else None
        row["chi2_p_core_vs_rare"] = chi
        rows.append(row)
    return pd.DataFrame(rows).sort_values("pct_annotated_80-100", ascending=False), n_ann


def bh(p):
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    if ok.sum():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[ok][order] if False else pv[order]
        n = len(pv)
        qv = np.empty(n)
        prev = 1.0
        for i in range(n - 1, -1, -1):
            prev = min(prev, ranked[i] * n / (i + 1))
            qv[order[i]] = prev
        q[ok] = np.clip(qv, 0, 1)
    return q


def baseline_share(tab, key):
    """class -> share of annotated HOGs, complex-wide (80-100% vs 0-20% columns are not used)"""
    n_ann = tab["n_HOGs_annotated_0-20"].sum()
    s = tab.set_index(key)["pct_annotated_20-80"].dropna()
    return s


def main():
    iv = pd.read_csv(IV, sep="\t")
    iv.columns = ["HOG", "rate", "bin"]
    iv = iv[iv.bin.isin(BINS)].copy()
    iv["grp3"] = iv.bin.map(collapse)
    print(f"orthogroups with an interval: {len(iv)}")

    gvh = pd.read_csv(GENE_HOG, sep="\t", usecols=["GeneID", "HOG_ID"])
    by_hog = gvh.groupby("HOG_ID").GeneID.apply(list).to_dict()

    # ---- Pfam: gene -> set of GO MF classes carried by its own domains ----
    lk = pd.read_csv(PFAM_LOOKUP, sep="\t", usecols=["pfam_acc", "go_mf_class"])
    gp = pd.read_csv(GENE_PFAM, sep="\t", header=None, usecols=[0, 1],
                     names=["GeneID", "pf"], dtype=str)
    gp["pfam_acc"] = gp.pf.str.split(".").str[0]
    gp = gp.merge(lk, on="pfam_acc", how="left")
    gp["mf"] = gp.go_mf_class.fillna("No assignment")
    pfam_sets = gp.drop_duplicates(["GeneID", "mf"]).groupby("GeneID").mf.apply(set).to_dict()

    # ---- KEGG: gene -> set of level-2 categories of its own pathways ----
    cat = pd.read_csv(KEGGCAT, sep="\t", dtype=str)
    l2 = {str(p).zfill(5): n for p, n in zip(cat.pathway_id, cat.level2)}

    def to_l2(v):
        if not isinstance(v, str) or v.strip() in ("", "-"):
            return set()
        out = set()
        for pid in v.split(","):
            m = re.match(r"^(?:ko|map|ec)?(\d{5})$", pid.strip())
            if m and m.group(1) in l2:
                out.add(l2[m.group(1)])
        return out

    eg = pd.read_csv(GENE_EGGNOG, sep="\t", header=None, usecols=[0, 5],
                     names=["GeneID", "path"], dtype=str).drop_duplicates("GeneID")
    kegg_sets = {g: s for g, s in zip(eg.GeneID, eg.path.map(to_l2)) if s}

    # ---- label every HOG once per axis ----
    pfam_lab, kegg_lab = {}, {}
    for hog, genes in by_hog.items():
        pfam_lab[hog] = hog_labels(pfam_sets, genes)
        kegg_lab[hog] = hog_labels(kegg_sets, genes)
    labels = {"go_mf_class": pd.Series(pfam_lab, name="go_mf_class"),
              "KEGG_level2": pd.Series(kegg_lab, name="KEGG_level2")}
    # one shared table, so that every other figure or table that needs a per-orthogroup
    # functional label reads the same one instead of re-deriving it with a different rule
    lab_out = pd.DataFrame({"HOG_ID": labels["go_mf_class"].index}).set_index("HOG_ID")
    lab_out["go_mf_class"] = labels["go_mf_class"]
    lab_out["KEGG_level2"] = labels["KEGG_level2"]
    lab_path = os.path.join(ROOT, "02.annotation/hog_function_labels.tsv")
    lab_out.reset_index().to_csv(lab_path, sep="\t", index=False)
    print(f"wrote {os.path.relpath(lab_path, ROOT)}  ({len(lab_out):,} orthogroups)")

    # ---- complex-wide tables ----
    tabs = {}
    for key, sheet in AXES:
        tab, n_ann = summarise(labels[key], iv, "No annotation", key)
        tabs[key] = tab
        print(f"{sheet}: annotated HOGs {n_ann['0-20']:,} / {n_ann['80-100']:,} "
              f"(0-20% / 80-100%)")

    # ---- per-species interval tables ----
    ivs = {}
    for sp in SPECIES:
        t = pd.read_csv(os.path.join(IV_SP, f"{sp}.HOG.interval.tsv"), sep="\t")
        t = t.rename(columns={"HOG_ID": "HOG", "Presence_Rate": "rate",
                              "Presence_Interval": "bin"})
        t["grp3"] = t.bin.map(collapse)
        ivs[sp] = t

    for key, sheet in AXES:
        wide = None
        for sp in SPECIES:
            tab, n_ann = summarise(labels[key], ivs[sp], "No annotation", key)
            tabs[f"{key}|{sp}"] = tab
            sub = tab[[key, "pct_annotated_0-20", "pct_annotated_20-80",
                       "pct_annotated_80-100", "fold_core_vs_rare"]].copy()
            sub.columns = [key] + [f"{sp}_{c}" for c in
                                   ("pct_0-20", "pct_20-80", "pct_80-100", "fold")]
            wide = sub if wide is None else wide.merge(sub, on=key, how="outer")
            print(f"{sheet}[{sp}]: annotated HOGs {n_ann['0-20']:,} / {n_ann['80-100']:,}")
        wide = wide.sort_values(f"{SPECIES[0]}_pct_80-100", ascending=False)
        tabs[sheet + "_species"] = wide

    # ---- candidates ----
    cand = {}
    for sp in SPECIES:
        f = pd.read_csv(os.path.join(CAND, f"rf_pav_{sp}_features.tsv"), sep="\t")
        cand[sp] = f.HOG_ID.tolist()
        print(f"candidates {sp}: {len(f)}")

    # interval placement of the candidates, on both presence definitions
    rows = []
    for sp in SPECIES:
        for scope, table, total_col in (
                ("complex (211 genomes)", iv, None),
                (f"{sp} only", ivs[sp], None)):
            base = table[["HOG", "bin"]]
            n_by_bin = base.bin.value_counts().reindex(BINS)
            n_tot = int(n_by_bin.sum())
            c_by_bin = (base[base.HOG.isin(cand[sp])].bin.value_counts()
                        .reindex(BINS).fillna(0).astype(int))
            n_c = int(c_by_bin.sum())
            for b in BINS:
                rows.append(dict(candidate_set=f"F. {sp}", presence_scope=scope,
                                 interval=b, n_candidates=int(c_by_bin[b]),
                                 pct_candidates=round(c_by_bin[b] / n_c * 100, 2),
                                 pct_background=round(n_by_bin[b] / n_tot * 100, 2),
                                 fold=round((c_by_bin[b] / n_c) / (n_by_bin[b] / n_tot), 2)
                                 if n_by_bin[b] else None))
    cand_int = pd.DataFrame(rows)

    # functional composition of the candidates against each baseline
    rows = []
    for key, sheet in AXES:
        broad = tabs[key]
        for sp in SPECIES:
            base_cx = broad.set_index(key)["pct_annotated_20-80"].dropna()
            base_sp = (tabs[f"{key}|{sp}"].set_index(key)["pct_annotated_20-80"].dropna())
            lab = labels[key]
            sel = lab.reindex(cand[sp]).dropna()
            n_c = len(sel)
            if not n_c:
                continue
            vc = sel.value_counts()
            ps, qs = [], []
            for k in vc.index:
                a = int(vc[k])
                b = int(round(base_sp.get(k, 0) / 100 * len(lab.dropna())))
                try:
                    ps.append(stats.chi2_contingency(
                        [[a, n_c - a], [b, int(len(lab.dropna())) - b]]).pvalue)
                except ValueError:
                    ps.append(np.nan)
            qs = bh(ps)
            for (k, a), p, q in zip(vc.items(), ps, qs):
                rows.append(dict(axis=sheet, candidate_set=f"F. {sp}", cls=k,
                                 n_candidates=int(a),
                                 pct_candidates=round(a / n_c * 100, 2),
                                 pct_complex_baseline=round(float(base_cx.get(k, np.nan)), 2),
                                 pct_species_baseline=round(float(base_sp.get(k, np.nan)), 2),
                                 fold_vs_complex=round((a / n_c * 100) / base_cx[k], 2)
                                 if k in base_cx and base_cx[k] else None,
                                 fold_vs_species=round((a / n_c * 100) / base_sp[k], 2)
                                 if k in base_sp and base_sp[k] else None,
                                 chi2_p_vs_species=p, bh_q_vs_species=q))
    cand_fun = pd.DataFrame(rows)

    overview = pd.read_csv(OVERVIEW, sep="\t")
    overview = overview[overview.Interval.isin(BINS)][
        ["Interval", "n_HOGs", "n_genes", "pct_effector_HOGs", "pct_effector_genes",
         "pct_HOGs_Pfam", "pct_HOGs_KEGG", "pct_HOGs_CAZy", "pct_HOGs_PHI"]]

    defs = pd.DataFrame([
        ("interval definition",
         "presence rate of each orthogroup over the 211 FSSC genomes, from "
         "05.stats_pan_by_species_shorter/all/N0.HOG.hog_count.interval.tsv; bins are "
         "0-20 / 20-40 / 40-60 / 60-80 / 80-100%, and the collapsed rows pool 20-80% because "
         "the three middle bins behave alike"),
        ("within-species intervals",
         "the *_species sheets use presence recomputed over each species' own genomes "
         "(!F. solani! 86 genomes, !F. falciformis! 66), from "
         "05.stats_pan_by_species_shorter/per_species/<sp>.HOG.interval.tsv. About one "
         "orthogroup in six changes bin when the genome set narrows to one species (83-84% stay "
         "put), and the near-core set grows from 12,853 complex-wide to 14,398 in !F. solani!, "
         "because a gene fixed within one species still scores low across all 211 genomes"),
        ("counting unit",
         "one orthogroup, one vote: each HOG is labelled by the majority class among its own "
         "annotated member genes and the intervals are compared over HOGs. Counting gene "
         "copies instead weights each interval by copy number, and the near-core interval "
         "holds 80% of all gene copies, so its composition converges on the genome-wide "
         "average and every contrast is compressed"),
        ("why the label is not taken from the HOG annotation field",
         "an orthogroup's Pfam field is the union over its members, so the number of domains "
         "it carries grows with orthogroup size (median 3 / 20 / 35 / 50 / 71 domains across "
         "the five bins) and a 'first matching keyword' rule then measures size rather than "
         "function: applied that way it labels 91.6% of core orthogroups as transporters, "
         "against 4.8% when each gene is counted once"),
        ("Pfam molecular-function axis",
         "each Pfam domain of each gene is mapped to a GO molecular-function / EC class through "
         "02.annotation/pfam_function_lookup.tsv (25,545 families), built by "
         "scripts/build_pfam_function_lookup.py from the Pfam 38.0 release files pfamA.txt and "
         "Pfam-A.clans.tsv: a family takes its clan's class where it has a clan and its own "
         "definition where it has none (13,524 families, 52.9%, have no clan). Accessions are "
         "matched on the 7-digit PF number"),
        ("KEGG axis",
         "each gene carries the level-2 categories of its own KEGG pathway IDs "
         "(02.annotation/all_isolates/all.fusarium.gene_eggnog.tsv, sixth column), mapped "
         "through 08.pan_020_genes/kegg_pathway_categories.tsv"),
        ("two denominators",
         "pct_* is the share of all HOGs of the interval; pct_annotated_* is the share of the "
         "annotated HOGs and is the one the Results text quotes. Both are given because "
         "annotation coverage itself tracks presence (fewer than half of the 0-20% HOGs carry a "
         "KEGG pathway, against all of the near-core HOGs), so dividing by all HOGs mixes "
         "composition with annotation. fold and chi-square use the annotated denominator"),
        ("candidates",
         "the clinical-associated orthogroups called by random forest within each species "
         "(390 for !F. solani!, 277 for !F. falciformis!, 12 shared; "
         "09.rf_clinical_genes/pav_rerun/rf_pav_<sp>_features.tsv). Their interval placement is "
         "reported under BOTH presence definitions because the two disagree for individual "
         "orthogroups, and the two must not be mixed. The 20-80% share of the matching species "
         "baseline is the comparator for the functional tests; BH correction is over classes "
         "within one axis and one species"),
        ("ties and missing votes",
         "a tie between classes is broken alphabetically; a HOG whose member genes carry no "
         "annotation at all is counted in 'No annotation' and excluded from pct_annotated_*"),
        ("resistance catalogue",
         "antifungal-resistance orthogroups were defined from the eggNOG preferred names, in four "
         "groups: azole core (ERG11, the lanosterol 14-alpha demethylase that azoles target, and "
         "UPC2, the sterol regulator that controls ERG expression); the ergosterol pathway around "
         "that target (ERG1-ERG13, ERG20, ERG24-ERG27, HMG1, IDI1, MVD1); efflux (CDR1 and SNQ2 "
         "as ABC transporters, MDR1 and ATR1 as MFS transporters, and the regulators PDR6 and "
         "PDR16); and other drug targets (FKS1 and its regulator RHO1 for the echinocandins, FUR1 "
         "for 5-fluorocytosine). No CARD or RGI run was used: CARD is centred on bacterial "
         "resistance and a negative result there would say nothing about a fungus, so the catalogue "
         "is the curated route"),
        ("resistance source cohort",
         "clinical versus non-clinical is taken from 09.rf_clinical_genes/rf_genes.per_genome.tsv, "
         "the same per-genome source and species table the random-forest screen uses, so that the "
         "two analyses run on one cohort: 60 clinical against 26 non-clinical genomes for "
         "!F. solani! and 48 against 18 for !F. falciformis!. Using the narrower source labels "
         "elsewhere in the project would leave 7 and 0 non-clinical genomes and no power"),
        ("resistance interval sheet",
         "resistance_genes gives, for each presence interval, the catalogue genes as a share of ALL "
         "genes in that interval, counted at gene level (one eggNOG-named gene counts once, and the "
         "denominator is every gene of the interval; the orthogroup-level share inflates with "
         "orthogroup size and is not used). The share is essentially flat: 0.15% of the genes in "
         "the 0-20% interval and 0.30% of those in the 80-100% interval, and 0.07% inside the "
         "random-forest candidate orthogroups against 0.29% elsewhere"),
        ("resistance association sheet",
         "resistance_source_assoc lists the catalogue orthogroups whose presence differs between "
         "clinical and non-clinical isolates of the same species (two-sided Mann-Whitney on the "
         "per-genome presence, BH across the 2,489 catalogue orthogroups within each species). "
         "Eighty-one orthogroups pass at q<0.05; clinical_rate and nonclinical_rate give the mean "
         "presence in each group and rate_delta their difference, so a negative delta marks an "
         "orthogroup more common outside hospitals"),
        ("effector interval sheet",
         "effector_genes gives, for each presence interval, the effector genes as a share of ALL "
         "genes in that interval, counted at GENE level from the EffectorP-3 consensus "
         "(08.pan_020_genes/effector_master.gene_level.tsv). EffectorP ran on the SignalP-positive "
         "secretome, so every scored protein carries an explicit call; the HOG-level share "
         "'orthogroups containing at least one effector gene' is not reported because it inflates "
         "with orthogroup size (82.1% of core orthogroups against 3.1% of core genes)"),
        ("effector per-genome sheet",
         "effector_per_genome gives the median number of effector genes per genome in clinical "
         "versus non-clinical isolates of each species (Mann-Whitney), on the random-forest cohort "
         "and excluding the four isolates with no EffectorP output (BMU12683-BMU12686), whose zero "
         "is a missing annotation rather than an absence of effectors"),    ], columns=["item", "value"])

    os.makedirs(os.path.dirname(OUT), exist_ok=True)

    # The workbook used to hold all of this at once, which mixed three different subjects: the
    # functional description of the intervals, the random-forest candidate sets, and the
    # antifungal-resistance catalogue. Each now has its own table and its own definitions page.
    s20_keys = ["interval definition", "counting unit",
                "why the label is not taken from the HOG annotation field",
                "Pfam molecular-function axis", "KEGG axis", "two denominators"]
    s21_keys = ["within-species intervals", "counting unit", "two denominators", "candidates",
                "ties and missing votes"]
    s22_keys = ["resistance catalogue", "interval definition", "resistance source cohort",
                "resistance interval sheet", "resistance association sheet"]

    def subset(keys):
        return defs[defs.item.isin(keys)].copy()

    with pd.ExcelWriter(OUT) as w:
        subset(s20_keys).to_excel(w, sheet_name="definitions", index=False)
        overview.to_excel(w, sheet_name="overview", index=False)
        tabs["go_mf_class"].to_excel(w, sheet_name="Pfam_MF", index=False)
        tabs["KEGG_level2"].to_excel(w, sheet_name="KEGG_level2", index=False)
        tabs["Pfam_MF_species"].to_excel(w, sheet_name="Pfam_MF_species", index=False)
        tabs["KEGG_level2_species"].to_excel(w, sheet_name="KEGG_level2_species", index=False)

    with pd.ExcelWriter(OUT_CAND) as w:
        subset(s21_keys).to_excel(w, sheet_name="definitions", index=False)
        cand_int.to_excel(w, sheet_name="candidates_intervals", index=False)
        cand_fun.to_excel(w, sheet_name="candidates_function", index=False)
        if os.path.exists(MARKER):
            pd.read_csv(MARKER, sep="\t").to_excel(w, sheet_name="marker_search", index=False)
        if os.path.exists(MARKER80):
            pd.read_csv(MARKER80, sep="\t").to_excel(w, sheet_name="marker_candidates_80_20", index=False)

    if os.path.exists(RES_IV) and os.path.exists(RES_TEST):
        res_iv = pd.read_csv(RES_IV, sep="\t")
        res_t = pd.read_csv(RES_TEST, sep="\t")
        res_sig = res_t[res_t.q_bh < 0.05].sort_values(["species", "q_bh"])
        with pd.ExcelWriter(OUT_RES) as w:
            subset(s22_keys).to_excel(w, sheet_name="definitions", index=False)
            res_iv.to_excel(w, sheet_name="resistance_genes", index=False)
            res_sig.to_excel(w, sheet_name="resistance_source_assoc", index=False)
            if os.path.exists(EFF_IV) and os.path.exists(EFF_TEST):
                pd.read_csv(EFF_IV, sep="\t").to_excel(w, sheet_name="effector_genes", index=False)
                pd.read_csv(EFF_TEST, sep="\t").to_excel(w, sheet_name="effector_per_genome", index=False)
        print(f"resistance table: {len(res_iv)} interval rows + "
              f"{len(res_sig)} source-associated orthogroups")
    print(f"\nwrote {os.path.relpath(OUT, ROOT)}")

    print("\ncandidate interval placement (complex-wide presence):")
    print(cand_int[cand_int.presence_scope == "complex (211 genomes)"].to_string(index=False))
    print("\ncandidate interval placement (own-species presence):")
    print(cand_int[(cand_int.presence_scope != "complex (211 genomes)")].to_string(index=False))


if __name__ == "__main__":
    main()
