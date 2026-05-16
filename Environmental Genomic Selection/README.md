# Environmental Genomic Selection

This directory contains tools and workflows for genomic prediction and selection based on environmental adaptation traits. These methods enable predictive breeding for environmental resilience and adaptation.

## Overview

Environmental Genomic Selection (EGS) combines genomic selection with environmental genomic data to identify and select individuals or breeding candidates showing superior adaptation to specific environmental conditions. This is increasingly important for breeding climate-adapted crops and populations.

## Contents

### Documentation

- **EGS.md** - Overview and conceptual framework for environmental genomic selection

- **README_FINAL.md** - Comprehensive guide including:
  - Detailed workflow description
  - Input data requirements
  - Parameter specifications
  - Output interpretation
  - Advanced usage examples

### Analysis Scripts

- **genomic_pipeline_FINAL.R** - Complete R pipeline for environmental genomic selection including:
  - Genomic relationship matrix calculation
  - Breeding value prediction
  - Environmental adaptation indices
  - Selection recommendations

- **use_genomic_selection_pipeline.R** - Simple wrapper script demonstrating basic usage of the genomic selection pipeline

## Key Features

- Genomic relationship matrix (GRM) estimation
- Breeding value prediction under environmental scenarios
- Integration of environmental adaptation data
- Selection index calculations
- Candidate ranking and recommendation

## Workflow Steps

1. Prepare genomic data and calculate GRM
2. Define environmental scenarios for selection
3. Estimate breeding values for environmental traits
4. Calculate genomic selection indices
5. Rank and select breeding candidates
6. Implement selections in breeding program

## Input Requirements

- Cleaned genomic data (genotypes)
- Environmental adaptation scores (from GEA analysis)
- Pedigree information (if available)
- Breeding objectives and environmental scenarios

## Output

Selection outputs include:
- Breeding values for environmental traits
- Genomic selection indices
- Candidate rankings
- Predicted genetic gain
- Selection recommendations

## Applications

- Breeding for climate adaptation
- Crop improvement for new environments
- Population-level adaptation planning
- Sustainable agricultural production

## Related Sections

See also:
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Identify environmental adaptation loci
- [Exploring Genomic Regions](../Exploring_genomic_regions_associated_with_environment/README.md) - Characterize environmental associations
- [Genetic Offset](../genetic_offset/README.md) - Assess maladaptation risk
