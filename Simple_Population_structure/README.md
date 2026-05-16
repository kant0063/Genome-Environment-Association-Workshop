# Simple Population Structure Analysis

This directory contains workflows for analyzing and visualizing population structure using Principal Component Analysis (PCA) and spatial population structure methods. Understanding population structure is essential for correcting for background population stratification in genome-environment association analyses.

## Overview

Population structure—differences in allele frequencies among subpopulations—can create spurious associations if not properly accounted for. This section provides tools for:
- Detecting and quantifying population structure
- Visualizing genetic differentiation among populations
- Estimating ancestry and admixture proportions
- Incorporating structure in association analyses

## Contents

### Documentation

- **Population_Structure_for_EAA_model.md** - Guide explaining:
  - Why population structure matters in GEA
  - How to assess structure in your data
  - Approaches for correcting structure
  - Integration with environmental models

### Analysis Scripts

- **Generic_PCA_Analysis.Rmd** - Basic R Markdown workflow for:
  - Principal Component Analysis (PCA)
  - Scree plot generation
  - PCA score visualization
  - Interpreting principal components

- **Spatial_Population_Structure.Rmd** - Advanced R Markdown workflow for:
  - Spatial analysis of population structure
  - Geography-corrected population structure
  - Spatial autocorrelation analysis
  - Cline analysis along environmental gradients

- **Spatial_Population_Structure_README.md** - Detailed documentation for spatial analyses including:
  - Methodology explanations
  - Parameter selection
  - Result interpretation
  - Visualization strategies

## Key Features

### PCA Analysis
- Eigenvalue and scree plot analysis
- PC score calculation and visualization
- Variance explained summaries
- Outlier identification

### Spatial Population Structure
- Geographic weighted analysis
- Spatial autocorrelation (Moran's I)
- Cline analysis along environmental gradients
- Isolation-by-distance assessment
- Spatial variance partitioning

## Workflow Steps

### Standard PCA Workflow
1. Prepare genotype data
2. Perform PCA
3. Calculate variance explained
4. Generate PC plots
5. Identify population clusters
6. Extract PC scores for use as covariates

### Spatial Population Structure Workflow
1. Prepare geographic and genotype data
2. Calculate pairwise geographic distances
3. Assess geographic vs. genetic structure
4. Fit spatial models
5. Perform cline analysis
6. Visualize spatial patterns
7. Quantify geographic signal

## Input Requirements

- Cleaned genomic data (genotypes or allele frequencies)
- Sample metadata
- Geographic coordinates (for spatial analysis)
- Population assignments (optional)

## Output

Results include:
- PC eigenvalues and scree plots
- PC score plots and bioplots
- Population structure visualizations
- Spatial structure plots and maps
- Variance explained estimates
- Structure metrics and statistics

## Interpreting PCA

- **PC1, PC2**: Usually capture major population structure/differentiation
- **Eigenvalues**: Variance explained by each PC
- **Scree plot**: Visual assessment of significant PCs
- **PC scores**: Used as covariates in association models
- **Clustering**: Separate populations visible as distinct clusters

## Interpreting Spatial Structure

- **Geographic signal**: Significant correlation between genetic and geographic distance
- **Clines**: Gradual changes in allele frequency along geographic or environmental gradients
- **Spatial autocorrelation**: Nearby samples more similar than distant samples
- **Isolation-by-distance**: Classic pattern of genetic structure following geography

## Correcting for Structure in GEA

Approaches include:
- **Include PC scores as covariates** in association models
- **Mixed models** that explicitly model structure
- **LFMM and RDA** that account for latent population structure
- **Conditional analysis** adjusting for significant PCs

## Related Sections

See also:
- [Cleaning Genomic Data](../Cleaning%20Genomic%20Data/README.md) - Prepare data for structure analysis
- [Mapping Accessions](../Mapping%20Accessions/README.md) - Visualize geographic context
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Integrate structure correction
- [Latent Factor Mixed Model and RDA](../Latent_Factor_Mixed_Model_and_RDA/README.md) - Alternative structure-aware methods
