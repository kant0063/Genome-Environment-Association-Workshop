# Cleaning Genomic Data

This directory contains documentation and guidelines for quality control and cleaning of genomic datasets before conducting genome-environment association analyses.

## Overview

Genomic data quality is essential for valid association study results. This section outlines the standard data cleaning procedures for removing low-quality variants, poorly genotyped individuals, and other quality control measures.

## Contents

### Files

- **Steps_for_cleaning_data.md** - Comprehensive guide detailing the step-by-step data cleaning workflow, including:
  - Quality control thresholds and criteria
  - Filtering strategies for SNP and individual data
  - Handling missing data and rare variants
  - Recommended tools and software

## Key Quality Control Steps

The cleaning workflow includes:

1. **Individual-level QC**
   - Remove individuals with high missing data rates
   - Filter by sequencing depth or genotyping rate
   - Check for duplicates and sample identity issues

2. **Variant-level QC**
   - Filter variants by call rate
   - Remove variants failing Hardy-Weinberg equilibrium tests
   - Filter by minor allele frequency (MAF) thresholds
   - Remove redundant or low-quality variants

3. **Population-level QC**
   - Check for cryptic relatedness
   - Account for population stratification
   - Assess linkage disequilibrium patterns

## Input Data

This workflow typically works with:
- VCF (Variant Call Format) files
- PLINK binary format files (.bed, .bim, .fam)
- Other standard genomic data formats

## Output

Cleaned genomic datasets ready for downstream analyses including:
- Association studies
- Population structure analysis
- Linkage disequilibrium estimation

## Related Sections

See also:
- [Linkage Disequilibrium](../Linkage_Disequilibrium/README.md) - Analyze LD patterns in cleaned data
- [Simple Population Structure](../Simple_Population_structure/README.md) - Assess population structure
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Use cleaned data for GEA
