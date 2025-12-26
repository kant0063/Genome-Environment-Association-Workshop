# LFMM: Genotype-Environment Association Analysis

A comprehensive R workflow for identifying loci associated with environmental variables while controlling for population structure using Latent Factor Mixed Models (LFMM).

## Overview

**LFMM** is a powerful statistical method for genome-environment association (GEA) studies that identifies SNPs correlated with environmental variables while accounting for neutral population structure and confounding effects.

### What This Analysis Does

1. Reads genotype data from VCF files
2. Matches samples with environmental data
3. Controls for population structure using latent factors
4. Tests each SNP for association with environmental variables
5. Applies multiple testing corrections (Bonferroni, FDR)
6. Generates publication-quality Manhattan plots
7. Produces comprehensive summary statistics

### Key Applications

- **Identify adaptive loci**: Find SNPs under environmental selection
- **Discover candidate genes**: Locate genes associated with adaptation
- **Climate change vulnerability**: Predict genetic maladaptation
- **Conservation genomics**: Guide management decisions
- **Breeding programs**: Select for climate resilience

---

## Quick Start

### Minimal Example

```r
# Set your file paths
VCF_FILE <- "my_genotypes.vcf"
ENV_CSV <- "environmental_data.csv"

# Run analysis
rmarkdown::render("LFMM_Analysis.Rmd")
```

### Expected Runtime

- **Small dataset** (50 samples, 10K SNPs): ~5-10 minutes
- **Medium dataset** (100 samples, 50K SNPs): ~30-60 minutes
- **Large dataset** (200+ samples, 100K+ SNPs): ~2-4 hours

---

## Installation

### Required R Packages

```r
# Core packages
install.packages(c("vcfR", "dplyr", "ggplot2", "tidyr"))

# LFMM package
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("lfmm")

# Optional: For enhanced tables and colors
install.packages(c("knitr", "kableExtra", "RColorBrewer"))
```

---

## Input Data Requirements

### 1. VCF File

Standard VCF format with quality-filtered biallelic SNPs.

### 2. Environmental CSV

Required: Sample IDs matching VCF + environmental variables as continuous values.

**Example:**
```csv
SampleID,Latitude,Longitude,Temperature,Precipitation
Sample1,45.5,-122.7,15.3,800
Sample2,46.8,-121.8,16.1,750
```

---

## Key Parameters

```r
K_VALUE <- 6                  # Latent factors (use PCA to inform)
P_VALUE_THRESHOLD <- 0.001    # Nominal significance
FDR_THRESHOLD <- 0.01         # FDR threshold
```

---

## Output Files

- `lfmm_all_variables.csv` - Combined results
- `lfmm_summary.csv` - Summary statistics
- `lfmm_manhattan_faceted.png` - Multi-panel Manhattan plots
- Individual files per environmental variable

---

## Citation

Caye, K., et al. (2019). LFMM 2: Fast and accurate inference of gene-environment associations in genome-wide studies. *Molecular Biology and Evolution*, 36(4), 852-860.

**Last updated:** December 2024
