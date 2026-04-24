# FOP Mouse Microbiome Taxonomic Analysis

This folder contains code for the taxonomic microbiome analyses associated with the submitted manuscript examining the gut microbiome in a mouse model of Fibrodysplasia Ossificans Progressiva (FOP).

## Data Generation

Raw sequencing data were processed using a shotgun metagenomic pipeline executed on UCSF's Wynton high-performance computing cluster. The full pipeline is available at: [https://github.com/ethan-dinh/HPC-Metagenomics/tree/main](https://github.com/ethan-dinh/HPC-Metagenomics/tree/main)

The pipeline performs the following steps in order:
1. **Quality trimming** — adapter and low-quality base removal with Trimmomatic
2. **Host decontamination** — removal of mouse (*C57BL/6NJ*) reads using KneadData
3. **Taxonomic classification** — read-level classification with Kraken2 against a reference database
4. **Abundance re-estimation** — Bracken re-estimates species-level abundances from Kraken2 reports

## Main Workflow

Open and knit: [analysis/FOP_Mouse_Microbiome_Nature_Comm.rmd](analysis/FOP_Mouse_Microbiome_Nature_Comm.rmd)

## Required Software

- R (4.2+ recommended)
- RStudio
- Key R packages: `phyloseq`, `vegan`, `DESeq2`, `taxize`, `ggplot2`, `ggpubr`, `ComplexHeatmap`

