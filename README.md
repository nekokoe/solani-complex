# FSSC Human Pathogenicity — Analysis Code

Analysis scripts for the manuscript:

> **The *Fusarium solani* species complex is an environmental opportunist whose clinical adaptation accompanies accessory-chromosome evolution**

Population-genomic and pan-genomic analysis of 211 *Fusarium solani* species complex (FSSC) isolates, including accessory-chromosome characterisation, clinical-association analysis, and experimental virulence testing in a *Galleria mellonella* model.

This repository contains only **analysis code** (scripts, R/Python/SHELL) for reproducibility and peer review. Input data are **not** included; see [DATA_ACCESS.md](DATA_ACCESS.md) for the data manifest and where to obtain the underlying datasets.

## Repository structure

```
code_repo/
├── scripts/
│   ├── 01_pan_genome/          # Pan-genome construction, clinical-association (RF), 0-20% gene characterisation
│   ├── 02_accessory/           # Accessory chromosome analysis: k-mer classification, gene transfer, TE mediation, repeat content, BMU02483 case study
│   ├── 03_virulence/           # Virulence phenotypes, literature gene screen, genome-wide association (GWA), temperature modules
│   ├── 04_population/          # Population genetics: AMOVA, PERMANOVA, ADMIXTURE plotting, PCoA
│   ├── 05_species_id/          # Species identification: barcode-gene ANI typing, multi-locus phylogenies (SHELL)
│   ├── 06_method_workflows/    # Method flowchart figures (S26/S27/S28) — grid-based R scripts
│   └── 07_suppl_tables_figures/# Supplementary-table assembly, enrichment statistics, figure utilities
├── README.md
├── DATA_ACCESS.md
└── requirements.txt
```

> Scripts synchronised 2026-10-27 (194 scripts; version 261027).

## Pipeline overview

```
External bioinformatics (not in this repo)
  OrthoFinder → N0.HOG.hog_count.tsv / N0.gene_vs_HOG.tsv
  mmseqs2     → all_solaniSC.cluster100_cluster.tsv
  RepeatMasker→ *.repeatmasker.out
  EffectorP-3 → per-isolate effectorp.tsv
  KMC/kmer-db → shared k-mer matrices (AC x isolates)
  plink       → strict core SNP genotypes

This repo: stats layer → plot layer → supplementary tables and figures
  stats_*.py / stats_*.R   → supplementary table TSVs
  plot_*.py / plot_*.R     → publication figures (PDF/PNG)
```

## Module pipelines

### 01 Pan-genome
```
stats_hog_count.py → HOG presence-rate intervals
stats_hog_curve.py → pan-genome/core-genome saturation curves → plot_hog_curve.py
build_effector_master.py → effector annotation master table
stats_020_genes.py → 0-20% presence-rate gene characterisation → plot_020_figure.py
binary_classification_feature_importance.py → random-forest clinical-association → feature_heatmap_clustering.R
stats_rf_genes.py / stats_rf_prevalence.py → RF-selected HOG statistics → plot_rf_v2.py / plot_rf_pies.py / plot_rf_heatmap.R
```

### 02 Accessory chromosomes
```
stats_ac_transfer_network.py → per-AC transfer profile
classify_kmer_type_v2.py → k-mer sharing type classification (block/mosaic)
stats_gene_transfer_211.py → 211-genome transfer detection + TE-flanking enrichment
stats_ac_hub_211.py → transfer hubs
stats_transfer_function.py → A-A vs C-A functional comparison
stats_repeat_terminal.py → chromosome/terminal repeat content → plot_repeat_terminal.py / ac_plot_repeat_panels_r.R
ac_bmu2483_gene_content.py → BMU02483 specific-HOG enrichment (Fisher)
ac_sliding_kmer_2483.py → 10-kb sliding-window k-mer sharing with BMU02483 (kmer-db)
plot_kmer_final.py → k-mer classification summary figures (F1/F2/F3)
ac_plot_transfer_mechanism_combined.R / ac_plot_transfer_function_combined.R / ac_plot_class_clinical_combined.R
  → combined supplementary figures S13/S16/S20
```

### 03 Virulence
```
phenotype_modeling.py → survival data → phenotype matrices (mortality/AUC/temperature preference)
literature_genes_phenotype.py → 16 literature pathogenicity genes (mmseqs easy-search)
genomewide_assoc.py → 3,630 variable HOGs x 9 traits, Spearman + permutation FDR (32,670 tests)
plot_phenotype_s25.py / plot_literature_s26.py / plot_gwa_s24.py / plot_hog_bias_s27v4.R
  → supplementary figures S25/S26/S24/S27
plot_virulence_fig8.R → main-text Figure 7 (virulence-genotype associations)
```

### 04 Population
```
prepare_amova.py → plink genotypes → 012 matrix → build_amova_dist.py → stats_amova.R (pegas)
prepare_permanova.py → PAV matrix → stats_permanova.R (vegan adonis2)
plot_admixture.R → ADMIXTURE Q-value plotting
stats_pcoa.py / plot_pcoa.py → PCoA of pan-genome features
```

### 05 Species identification
```
ani_typer/barcode_genes/.../run_ani_typer.*.sh   → barcode-gene ANI typing (BLAT, FastANI, ani-typer)
phylo_trees/.../build_tree.*.sh                 → single/multi-locus ML phylogenies (MAFFT, trimAl, IQ-TREE)
```
Note: these SHELL scripts reference data directories via environment variables (`FSSC_ROOT`, `FSSC_REFERENCE`, `FSSC_DOWNLOAD`, `FSSC_ANNOTATION`, `ANI_TYPER_BIN`); set them before running.

## Key external dependencies

| Tool | Version | Reference |
|---|---|---|
| Python | 3.13 | — |
| R | 4.5.1 | — |
| pandas / scipy / numpy / matplotlib / scikit-learn | — | McKinney 2010 / Virtanen 2020 / Hunter 2007 / Pedregosa 2011 |
| ggplot2 / patchwork / ComplexHeatmap / ggtree / ggalluvial / vegan / pegas | — | Wickham 2016 / Oksanen / Paradis 2010 |
| mmseqs2 | v14.7 | Steinegger & Söding 2017 |
| kmer-db | v1.11.1 | Deorowicz et al. 2019 |
| OrthoFinder | v2.5.4 | Emms & Kelly 2019 |
| EffectorP | 3.0 | Sperschneider & Dodds 2022 |
| BLAT / MAFFT / trimAl / IQ-TREE / FastANI / ADMIXTURE / PLINK / VCFtools / KMC | — | see manuscript Methods |

Full software citation list appears in the manuscript reference list.

## Reproducibility notes

- All Python scripts use `BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` and `os.chdir(BASE)`, so they run from any clone location.
- Figures are written next to the source data in the original project layout; adjust output paths if your layout differs.
- Supplementary data tables generated by the `stats_*` scripts are listed in [DATA_ACCESS.md](DATA_ACCESS.md).

## License / contact

Code is provided for peer-review purposes. For questions contact the corresponding authors listed in the manuscript.
