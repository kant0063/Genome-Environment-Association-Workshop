# Genetic Offset Analysis

This directory contains tools and workflows for assessing genetic offset—the mismatch between a population's current genetic adaptation and the environmental conditions it is projected to experience under future climate scenarios.

## Overview

Genetic offset quantifies the risk of maladaptation to future environments by measuring how much an organism's current genetic composition is mismatched with projected environmental conditions. This is critical for:
- Assessing climate change vulnerability
- Planning conservation and breeding strategies
- Identifying populations at risk
- Informing climate adaptation policy

## Contents

### Documentation

- **Genetic_Offset_README.md** - Comprehensive guide including:
  - Conceptual framework and definition of genetic offset
  - Methodology and statistical approaches
  - Input data requirements and preparation
  - Result interpretation and visualization
  - Case studies and applications

### Analysis Scripts

- **Genetic_Offset_Analysis.Rmd** - Main R Markdown workflow for:
  - Calculating genetic offsets
  - Assessing maladaptation risk
  - Visualizing offset scores
  - Comparing populations and scenarios
  - Generating summary statistics

### Sub-directory

- **koppen_geiger_future_climate_analysis/** - Specialized tools for analyzing climate classifications
  - See [koppen_geiger_future_climate_analysis/README.md](./koppen_geiger_future_climate_analysis/README.md)

## Key Features

- Quantifies mismatch between current and future environments
- Population-level risk assessment
- Multiple climate scenarios support
- Visualization of offset patterns
- Comparative population analysis
- Integration with GEA results

## Workflow Steps

1. Extract environmental association scores (from GEA analysis)
2. Obtain current environmental conditions
3. Obtain projected future climate scenarios
4. Calculate current genetic offset (validation)
5. Calculate projected genetic offset
6. Assess offset magnitude and uncertainty
7. Visualize and interpret results
8. Generate recommendations

## Input Requirements

- SNP allele frequencies or genotypes
- Environmental association statistics (β, z-scores, or effect sizes)
- Current environmental data
- Future climate projections (CMIP5, CMIP6, or other models)
- Population geographic coordinates

## Output

Results include:
- Genetic offset scores by population
- Comparison matrices between populations
- Risk classifications (low, medium, high)
- Uncertainty estimates
- Visualization plots and maps
- Summary tables and statistics

## Key Concepts

### Genetic Offset Calculation

Genetic offset typically calculated as:
- Distance between current allele frequencies and those predicted for future environment
- Weighted by environmental association effect sizes
- Measured in allele frequency units or phenotypic space

### Interpretation

- **Low offset**: Population well-adapted to future conditions
- **High offset**: Population potentially maladapted to future conditions
- **Moderate offset**: Variable expected (some loss of fitness)

## Climate Scenarios

Typically incorporates:
- **Multiple GCMs** (General Circulation Models)
- **Multiple RCPs/SSPs** (emission scenarios)
- **Multiple time periods** (2050s, 2100s, etc.)

## Applications

### Conservation
- Prioritize vulnerable populations for protection
- Plan assisted migration or translocation
- Inform seed sourcing and provenance strategies

### Agriculture/Breeding
- Identify germplasm needing climate adaptation
- Plan breeding for future climates
- Select environmental resilience traits

### Policy
- Climate change vulnerability assessments
- Adaptation planning
- Species/population-level risk evaluations

## Limitations and Considerations

- Based on current genetic associations (may change)
- Assumes linear relationships with environment
- Cannot predict new mutations or adaptation
- Depends on accuracy of climate projections
- Requires accurate environmental association estimates

## Climate Classification Analysis

The **koppen_geiger_future_climate_analysis** subdirectory provides additional tools for analyzing how Köppen-Geiger climate classifications may shift under future scenarios, offering complementary ecological interpretation.

## Related Sections

See also:
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Obtain environmental association scores
- [Exploring Genomic Regions](../Exploring_genomic_regions_associated_with_environment/README.md) - Understand associated genes
- [Environmental Genomic Selection](../Environmental%20Genomic%20Selection/README.md) - Breeding for adaptation
- [Mapping Accessions](../Mapping%20Accessions/README.md) - Geographic context for populations
