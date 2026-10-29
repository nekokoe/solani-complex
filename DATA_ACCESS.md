# Data Access Manifest

This repository contains **code only**. Underlying input data are large (tens of GB) and are hosted separately (unpublished, pending GenBank upload). This manifest lists the key inputs and outputs referenced by the scripts, with locations and how to obtain them.

## Input datasets (needed to reproduce analyses)

| Dataset | Location (in the original project) | Size | Notes |
|---|---|---|---|
| 211 isolate metadata | `sample.solaniSC.211.list`, `bmu_ref_groups.211.txt` | KB | clinical/environmental labels, species |
| HOG count matrix | `05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv` | ~150 MB | from OrthoFinder N0 output |
| Gene-to-HOG mapping | `05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv` | ~170 MB | |
| HOG functional annotation | `05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv` | ~330 MB | eggNOG/Pfam annotations |
| mmseqs protein clusters (100% identity) | `06.accessory-chrs/tgs_assemblies/all_solaniSC.cluster100_cluster.tsv` | ~178 MB | 211-isolate clusters |
| RepeatMasker outputs | `02.annotation/per_isolates/<strain>/.../fusarium.repeatmasker.out` | — | per-isolate |
| EffectorP predictions | `02.annotation/per_isolates/<strain>/signalp/effectorp.tsv` | — | |
| k-mer databases | `06.accessory-chrs/kmer_analysis/` | GB | KMC + kmer-db outputs |
| Core SNP genotypes | `04.population/*/...strict_core_snps.plink` | 10 GB+ | VCF/PLINK per species |
| TGS accessory chromosome definitions | `06.accessory-chrs/tgs_assemblies/*.accessory_chrs.id` | KB | 13 chromosome-level assemblies |
| Survival data (Galleria) | `07.lethal_groups/caterpillar_37C_narrow.txt` | KB | 15 strains x 6 conditions |
| Literature pathogenicity genes | `05.pan-genome/pathogenicity_genes/*.fasta` | KB | 16 genes |

**Large raw data** (genome assemblies, raw VCFs, cluster matrices >100 MB) are retained locally and are **not** intended for the GitHub repository. They are available from the corresponding authors on request, and will be deposited (BMU strain genomes are pending upload to GenBank; public reference genomes are listed in the manuscript Methods).

## Key output tables (supplementary tables)

These are the analysis results referenced by the manuscript; they are small and included in the manuscript's Supplementary Information:

| Supplementary table | Produced by | Content |
|---|---|---|
| Supplementary Table S1/S2 | (author) | isolate metadata / genome statistics |
| Supplementary Table S7 | genomewide_assoc.py + plot_gwa_s24.py | temperature modules + nominal candidates |
| Supplementary Table S10 | stats_gene_transfer_211.py | repeat-flanking enrichment |
| Supplementary Table S11 | ac_bmu2483_gene_content.py | BMU02483 gene-content evidence |
| kmer AC distribution | classify_kmer_type_v2.py / plot_kmer_final.py | AC k-mer classification + clinical association |
| transfer function tables | stats_transfer_function.py | A-A vs C-A functional enrichment |

## Environment variables for SHELL scripts (05_species_id)

The SHELL scripts expect these to point to your data tree:

```
export FSSC_ROOT=/path/to/project_root
export FSSC_REFERENCE=/path/to/fusarium.org_reference
export FSSC_DOWNLOAD=/path/to/persoonia_downloads
export FSSC_ANNOTATION=/path/to/per_isolate_assemblies
export ANI_TYPER_BIN=/path/to/ani-typer/bin
```

## Correspondence

Data requests: see corresponding authors in the manuscript.
