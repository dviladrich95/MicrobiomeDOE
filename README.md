# Exploratory Microbiome Data Analysis for Design of Experiments Research

This project focuses on exploratory data analysis of microbiome datasets in the context of **Design of Experiments (DoE)**. The primary aim is to investigate how microbiome composition varies across different body mass index (BMI) groups and nationalities, using established diversity and ordination techniques.

## Dataset

The analysis is based on the `dietswap` dataset, available through the [**microbiome**](https://microbiome.github.io) and [**phyloseq**](https://joey711.github.io/phyloseq/) R packages. This dataset includes microbiome profiles of individuals categorized by:
- **BMI group**: lean, overweight, and obese
- **Nationality**: native African (AFR) and African American (AAM)

## Analyses Performed

- **Relative abundance heatmaps**  
  Visualizations of genus-level abundances across BMI groups and nationalities, aggregated by timepoint.
  
- **Alpha diversity (Shannon index)**  
  Diversity within samples was assessed to understand richness and evenness across BMI categories and nationalities.

- **Beta diversity (PCoA)**  
  Principal Coordinates Analysis (based on Bray-Curtis distances) was used to explore compositional differences between samples from native African and African American individuals.

## Goals

This analysis serves as a foundation for integrating microbiome data into experimental design research, highlighting key differences in microbial community structures that may inform future interventions or studies.
