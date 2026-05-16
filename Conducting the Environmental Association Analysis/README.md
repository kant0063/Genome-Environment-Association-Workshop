# Conducting the Environmental Association Analysis

This directory contains comprehensive scripts and documentation for performing genome-environment association (GEA) analyses, including GWAS with environmental covariates, heritability estimation, and visualization of results.

## Overview

Genome-Environment Association (GEA) analysis tests for statistical associations between genetic variants and environmental variables. This workflow implements multiple approaches for identifying genomic regions under environmental selection.

## Contents

### Documentation

- **Conducting_the_EAA.md** - Overview and guide for conducting environmental association analyses

### Analysis Scripts

- **calculate_BLUEs_GWAS.Rmd** - R Markdown script for calculating Best Linear Unbiased Estimators (BLUEs) and conducting GWAS with environmental covariates

- **environmental_gwas.Rmd** - Main R Markdown workflow for conducting genome-wide association scans with environmental variables

- **environmental_gwas_from_vcf.R** - R script for performing environmental GWAS directly from VCF files

- **heritability_analysis.Rmd** - R Markdown script for estimating heritability and partitioning genetic variance explained by environmental factors

- **visualize_egwas_results.R** - Comprehensive R script for creating publication-quality visualizations of GEA results including:
  - Manhattan plots
  - Q-Q plots
  - Regional plots
  - Effect size summaries

## Key Features

- Multiple statistical frameworks for GEA (linear models, mixed models)
- Integration of environmental covariates
- Heritability and variance component estimation
- Correction for population structure
- Comprehensive result visualization
- Direct VCF input support

## Workflow Steps

1. Prepare cleaned genomic and environmental data
2. Calculate BLUEs or phenotypic values (if applicable)
3. Run environmental GWAS using preferred method
4. Estimate heritability and effect sizes
5. Visualize and interpret results

## Input Requirements

- Cleaned genomic data (VCF, PLINK format, or genotype matrix)
- Environmental variables (raster data, point data, or covariates)
- Sample metadata and geographic information
- Phenotypic data (if conducting phenotype-environment associations)

## Output

Results include:
- GWAS summary statistics
- Significant SNP associations
- Manhattan and Q-Q plots
- Regional association plots
- Heritability estimates

## Related Sections

See also:
- [Accessing Climate and Soil Data](../Accessing%20Climate%20and%20Soil%20Data/README.md) - Obtain environmental variables
- [Cleaning Genomic Data](../Cleaning%20Genomic%20Data/README.md) - Prepare genomic data
- [Simple Population Structure](../Simple_Population_structure/README.md) - Assess population stratification
- [Exploring Genomic Regions](../Exploring_genomic_regions_associated_with_environment/README.md) - Interpret results
