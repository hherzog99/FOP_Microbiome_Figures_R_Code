# FOP_MICROBIOME_PAPER_FIGURES
R scripts, processed data, and reproducible analysis workflows supporting the submitted manuscript on the gut microbiome Fibrodysplasia Ossificans Progressiva (FOP) paper.
# FOP_MICROBIOME_MOUSE_PAPER

R scripts, processed data, and reproducible analysis workflows supporting the submitted manuscript investigating the gut microbiome in a mouse model of Fibrodysplasia Ossificans Progressiva (FOP).

# Overview

This repository contains the processed datasets, statistical code, and figure-generation scripts used for analyses associated with the submitted manuscript examining how gut microbiome manipulation influences inflammation and heterotopic ossification in a mouse model of FOP.

All scripts are written in R and organized to allow reviewers to inspect and reproduce the reported analyses.

# Repository Structure

## data/

Processed input datasets used for analysis.

Files include:

- `anti_il_1_ho_quants.xlsx`  
  Anti–IL-1 treatment heterotopic ossification quantification data.

- `Cytokine.xlsx`  
  Serum cytokine dataset used for heatmap and statistical analyses.

- `FOP_ABX_KM_Curve.xlsx`  
  Survival data used for Kaplan–Meier analysis.

- `ho_quantification.xlsx`  
  Post-GMA heterotopic ossification volume dataset.

## scripts/
R scripts used to perform statistical analyses and generate figures.

- `Anti-IL1 HO Quantification.R`
- `Cytokine Heatmap.R`
- `FOP GMA HO Quants.R`
- `FOP GMA KM Curve.R`

## figures/

Output figures and summary statistics generated from the scripts.

Includes:

- HO quantification plots (`.pdf`, `.tiff`)
- Cytokine heatmaps
- Kaplan–Meier survival curves
- Statistical summary spreadsheets

## results/

Optional folder for additional outputs or intermediate processed files.

# Reproducibility

All scripts use relative paths managed through the `here` package.

When the repository folder structure is preserved, scripts should run without manually editing file paths.

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

This repository contains processed experimental data only. No raw sequencing files are included.


# Contact

For questions regarding analysis workflows or data organization, please contact the corresponding author listed in the submitted manuscript.
