# Linkage Disequilibrium Analysis

This directory contains tools and workflows for analyzing linkage disequilibrium (LD) patterns in genomic datasets, including LD decay analysis and visualization.

## Overview

Linkage Disequilibrium (LD) describes the non-random association of alleles at different loci. Understanding LD structure is essential for:
- Interpreting GEA results and defining independent signals
- Choosing appropriate SNP pruning strategies
- Understanding recombination patterns
- Planning fine-mapping studies
- Conducting genome-wide association analyses

## Contents

### Documentation

- **sample_LD_calc.md** - Quick reference guide for LD calculation

### Analysis Scripts

- **LD_decay_analysis.Rmd** - Comprehensive R Markdown workflow for:
  - Calculating pairwise LD statistics (r², D', etc.)
  - Analyzing LD decay with physical distance
  - Estimating effective recombination rates
  - Visualizing LD patterns
  - Chromosome-specific and genome-wide analyses

## Key Features

- Efficient LD matrix calculation
- LD decay curve fitting and analysis
- Visualization of LD structure
- Genome-wide and regional LD analysis
- Multiple LD statistics supported
- Population-specific LD estimation

## LD Statistics

### Common LD Measures
- **r²** - Squared correlation coefficient (most commonly used)
- **D'** - Lewontin's standardized disequilibrium (normalized)
- **D** - Disequilibrium coefficient

### Interpretation
- r² = 1: Perfect LD (complete association)
- r² = 0: No LD (independent segregation)
- r² > 0.3-0.5: Moderate LD (commonly used thresholds)

## Workflow Steps

1. Load cleaned genomic data
2. Define LD analysis parameters (window size, thresholds)
3. Calculate pairwise LD statistics
4. Analyze LD decay with distance
5. Identify LD blocks
6. Visualize results
7. Interpret for downstream analyses

## Input Requirements

- Cleaned genomic data (VCF, PLINK format, or genotype matrix)
- SNP position information
- Sample metadata

## Output

Results include:
- LD matrices
- LD decay curves and parameters
- LD block definitions
- Visualization plots
- Summary statistics

## Applications

- **SNP pruning** - Select independent variants for analysis
- **Fine-mapping** - Determine independent causal variants
- **Population genetics** - Understand recombination and demographic history
- **GWAS interpretation** - Distinguish independent signals
- **Reference panel selection** - Choose tag SNPs

## LD Patterns and Interpretation

- **High LD**: Limited recombination or recent mutations
- **Low LD**: Ancient recombination or high mutation rates
- **LD blocks**: Regions with high internal LD
- **LD decay**: Rate of LD decline with distance reflects:
  - Population history (bottlenecks, admixture)
  - Recombination rates
  - Demographic processes

## Related Sections

See also:
- [Cleaning Genomic Data](../Cleaning%20Genomic%20Data/README.md) - Prepare data for LD analysis
- [Simple Population Structure](../Simple_Population_structure/README.md) - Population-specific LD patterns
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Use LD for result interpretation
