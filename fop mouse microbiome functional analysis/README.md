# FOP Mouse Microbiome Functional Analysis (HUMAnN)

This repository contains code for the functional metagenomic analyses associated with the submitted manuscript examining the gut microbiome in a mouse model of Fibrodysplasia Ossificans Progressiva (FOP).

## Data Generation

Raw sequencing data were processed using a shotgun metagenomic pipeline executed on UCSF's Wynton high-performance computing cluster. The full pipeline is available at: [https://github.com/ethan-dinh/HPC-Metagenomics/tree/main](https://github.com/ethan-dinh/HPC-Metagenomics/tree/main)

The pipeline performs the following steps in order:
1. **Quality trimming** — adapter and low-quality base removal with Trimmomatic
2. **Host decontamination** — removal of mouse (*C57BL/6NJ*) reads using KneadData
3. **Functional profiling** — gene family and pathway abundance estimation with HUMAnN3 against the UniRef90 protein database
4. **KO annotation** — KEGG Orthology (KO) mapping via HUMAnN3's regrouping utility

## Main Workflow

Open and knit: [analysis/FOP_Functional_analysis_NatureComm.rmd](analysis/FOP_Functional_analysis_NatureComm.rmd)

## Required Software

- R (4.2+ recommended)
- RStudio
- Key R packages: `DESeq2`, `ggplot2`, `ggpubr`, `ComplexHeatmap`, `dplyr`, `tidyr`

All required R packages are installed automatically in the setup section of the notebook.

## Notes

All file paths are relative and should run when the folder structure is preserved. Set the working directory to the `analysis/` folder before knitting.
