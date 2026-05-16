# Latent Factor Mixed Model and Redundancy Analysis (RDA)

This directory contains implementations of advanced multivariate statistical approaches for genome-environment association analysis, including Latent Factor Mixed Models (LFMM) and Redundancy Analysis (RDA).

## Overview

Latent Factor Mixed Models and Redundancy Analysis are powerful multivariate techniques for detecting associations between genetic variation and environmental variables while accounting for population structure and latent factors. These approaches complement univariate GEA analyses.

## Contents

### Documentation

- **LFMM_README.md** - Guide to Latent Factor Mixed Models including:
  - Theoretical background
  - Parameter selection
  - Result interpretation
  - Comparison with other approaches

- **RDA.md** - Overview of Redundancy Analysis framework

- **RDA_tutorial.md** - Detailed RDA tutorial with:
  - Step-by-step workflow
  - Interpretation guidelines
  - Example analyses
  - Visualization approaches

### Analysis Scripts

- **LFMM_Analysis.Rmd** - R Markdown workflow for:
  - LFMM model fitting
  - Latent factor estimation
  - Environmental association testing
  - Result visualization

- **RDA_tutorial.Rmd** - R Markdown implementation of:
  - RDA analysis pipeline
  - Canonical correspondence analysis
  - Variance partitioning
  - Significance testing

- **Generic_WZA_Analysis.Rmd** - Weighted Z-score analysis and advanced association methods

## Key Features

### LFMM
- Accounts for hidden population structure through latent factors
- Performs genome-wide environmental scans
- Provides association statistics (z-scores)
- Minimal parameter tuning required

### RDA
- Multivariate approach relating allele frequencies to environmental variables
- Variance partitioning between environmental and spatial components
- Visualization of environmental gradients
- Identifies principal environmental gradients

## Workflow Steps

### LFMM Workflow
1. Prepare genotype and environmental data
2. Determine number of latent factors
3. Fit LFMM model
4. Calculate z-scores for associations
5. Determine significance thresholds
6. Visualize and interpret results

### RDA Workflow
1. Format allele frequency or SNP data
2. Prepare environmental variables
3. Perform RDA analysis
4. Partition variance
5. Test significance of axes
6. Visualize ordination and associations

## Input Requirements

- Genotype data (SNP matrix or allele frequencies)
- Environmental variables (quantitative or categorical)
- Sample metadata
- Population information (optional, for stratification analysis)

## Output

Results include:
- Association z-scores (LFMM)
- Ordination plots (RDA)
- Variance partitioning tables
- Significantly associated SNPs
- Environmental gradient visualizations

## Advantages

- **LFMM**: Better handling of hidden structure, single model fitting
- **RDA**: Multivariate interpretation, variance partitioning, visual diagnostics

## Comparison with Univariate GEA

These multivariate approaches:
- Control for cryptic population structure more effectively
- Enable simultaneous analysis of multiple environmental variables
- Provide insights into redundancy and collinearity
- Often have greater statistical power for complex environmental scenarios

## Related Sections

See also:
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Univariate GEA approaches
- [Simple Population Structure](../Simple_Population_structure/README.md) - Population stratification assessment
- [Exploring Genomic Regions](../Exploring_genomic_regions_associated_with_environment/README.md) - Interpretation of results
