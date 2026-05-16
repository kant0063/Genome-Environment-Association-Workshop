# Köppen-Geiger Future Climate Analysis

This directory contains specialized tools for analyzing how Köppen-Geiger climate classifications are projected to shift under future climate scenarios, providing ecological context for genetic offset assessments.

## Overview

The Köppen-Geiger climate classification system provides an intuitive ecological framework for understanding how environments may change. This workflow analyzes shifts in climate classifications under future scenarios to complement genetic offset analyses with ecological interpretation.

## Contents

### Documentation

- **README.md** - Overview of Köppen-Geiger classification system and analysis approach (this file)

- **VISUALIZATION_GUIDE.md** - Comprehensive guide for:
  - Creating climate classification maps
  - Visualizing classification shifts
  - Generating publication-quality figures
  - Color schemes and aesthetic choices
  - Interpretation of visual patterns

### Analysis Scripts

- **koppen_geiger_climate_change_analysis.R** - R script for:
  - Assigning Köppen-Geiger classifications
  - Comparing current vs. future classifications
  - Identifying climate shifts
  - Calculating areas of classification change
  - Visualizing climate scenario comparisons
  - Generating summary statistics

## Key Features

- Classification assignment based on bioclimatic variables
- Temporal comparison (current vs. future scenarios)
- Area-based analysis of climate shifts
- Multi-scenario support and comparison
- Publication-ready visualization functions
- Quantitative summary metrics

## Köppen-Geiger System Overview

The Köppen-Geiger classification uses:
- **Temperature** (annual, seasonal)
- **Precipitation** (annual, seasonal)
- Combines into hierarchical categories:
  - **Tropical** (A) - Hot climates
  - **Dry** (B) - Arid and semi-arid
  - **Temperate** (C) - Warm and cool temperate
  - **Continental** (D) - Cold with winter snow
  - **Polar** (E) - Extremely cold

## Workflow Steps

1. Obtain current and projected bioclimatic variables
2. Calculate Köppen-Geiger classifications
3. Compare current vs. future classifications
4. Identify areas of climate change
5. Quantify classification shifts by type
6. Analyze geographic patterns
7. Visualize changes and create maps
8. Generate summary reports

## Input Requirements

- Bioclimatic variables for current period
- Bioclimatic variables for future climate scenarios
- Climate projections (WORLDCLIM, downscaled GCMs, etc.)
- Geographic coordinate system information
- Raster or point-based environmental data

## Output

Results include:
- Köppen-Geiger classification maps (current and future)
- Climate change classification maps
- Statistical summaries of classification shifts
- Area affected by climate changes
- Transition matrices showing classification changes
- Publication-quality visualizations

## Interpreting Climate Shifts

### Types of Changes
- **Within-type variation** - Changes within category (e.g., wetter savanna to drier forest)
- **Category changes** - Shifts between major classes (e.g., temperate to continental)
- **Zone expansion** - Poleward or altitudinal shift of climate boundaries
- **Range contraction** - Loss of suitable climate space

### Ecological Implications

**Tropical to Dry**: 
- Decreased precipitation
- Potential for desertification
- Shift from forest to savanna ecosystems

**Temperate to Continental**:
- More extreme winters
- Greater seasonal variation
- Different frost and freeze patterns

**Dry to Temperate**:
- Increased precipitation and moisture
- Potential for vegetation expansion
- Greening and biomass increases

## Applications

### Conservation Biology
- Identify refugial climate space
- Plan species range shift corridors
- Assess habitat suitability changes

### Agriculture
- Determine crop suitability shifts
- Plan regional adaptation strategies
- Identify emerging growing zones

### Biodiversity Planning
- Assess ecosystem stability
- Identify high-priority conservation areas
- Plan for climate-responsive management

### Complementing Genetic Offset

This analysis provides:
- **Ecological context** for genetic offset scores
- **Climate analog** identification (finding current conditions like future climates)
- **Population-level interpretation** of offset in ecological terms
- **Cross-validation** of offset predictions

## Visualization Guidelines

See **VISUALIZATION_GUIDE.md** for:
- Color scheme recommendations
- Map projection choices
- Effective figure layouts
- Adding geographic context
- Creating publication figures
- Interpretation aids and legends

## Limitations and Considerations

- Köppen-Geiger is discrete classification (hides gradual changes)
- Sensitive to bioclimatic variable thresholds
- Does not capture within-classification variation
- Climate projections have inherent uncertainty
- Ecological changes may lag climate changes

## Integration with Genetic Analysis

Use climate classification shifts alongside genetic offset to:
- Validate offset predictions with ecological interpretation
- Identify populations in novel future climates (no current analog)
- Assess whether genetic adaptation matches ecological expectations
- Communicate findings to non-specialist audiences
- Develop management and policy recommendations

## Related Sections

See also:
- [Genetic Offset Analysis](../README.md) - Main genetic offset workflow
- [Accessing Climate and Soil Data](../../Accessing%20Climate%20and%20Soil%20Data/README.md) - Obtain bioclimatic data
- [Mapping Accessions](../../Mapping%20Accessions/README.md) - Geographic visualization
- [Environmental Genomic Selection](../../Environmental%20Genomic%20Selection/README.md) - Adaptation planning
