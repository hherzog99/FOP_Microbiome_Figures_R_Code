# FOP_MICROBIOME_PAPER_FIGURES
R scripts supporting the submitted manuscript on the gut microbiome Fibrodysplasia Ossificans Progressiva (FOP) paper.
# FOP_MICROBIOME_MOUSE_PAPER

R scripts supporting the submitted manuscript investigating the gut microbiome in a mouse model of Fibrodysplasia Ossificans Progressiva (FOP).

# Overview

This repository contains statistical code and figure-generation scripts used for analyses associated with the submitted manuscript examining how gut microbiome manipulation influences inflammation and heterotopic ossification in a mouse model of FOP.

All scripts are written in R and organized to allow reviewers to inspect the reported analyses.

# Repository Structure

## scripts/
R scripts used to perform statistical analyses and generate figures.

- `Anti-IL1 HO Quantification.R` --> Figure 4B
- `Cytokine Heatmap.R`--> Figure 3B
- `FOP GMA HO Quants.R` --> Figure 2D
- `FOP GMA KM Curve.R` --> Figure 2F

## figures/

Output figures and summary statistics generated from the scripts.

Includes:

- HO quantification plots
- Cytokine heatmaps
- Kaplan–Meier survival curves
- Statistical summary spreadsheets

# Reproducibility

All scripts use relative paths managed through the `here` package.

Recommended workflow:

1. Download or clone the repository  
2. Open `microbiome_paper.Rproj` in RStudio  
3. Open any script in `/scripts`  
4. Run the script

Figures and exported outputs will be written to the `/figures` folder.

# Required Software

- R (version 4.2 or newer recommended)
- RStudio (recommended)

Required packages are loaded and/or installed within scripts where applicable.

Common packages include:

- `here`
- `tidyverse`
- `ggplot2`
- `readxl`
- `survival`
- `survminer`
- `ComplexHeatmap`
- `pheatmap`
- `writexl`


# Notes

This repository contains R code only. No raw data or sequencing files are included.

