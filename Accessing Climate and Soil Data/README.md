# Accessing Climate and Soil Data

This directory contains resources and scripts for downloading and processing climate and soil environmental variables for use in genome-environment association analyses.

## Overview

Accessing quality environmental data is a critical first step in conducting genome-environment association (GEA) analyses. This section provides tools and guidance for obtaining bioclimatic variables, soil properties, and other environmental covariates.

## Contents

### Files

- **Get_bioclim_data_generic.Rmd** - R Markdown script for downloading bioclimatic data (WORLDCLIM2) for any geographic region. Provides a generic workflow that can be adapted to your study system.

- **climate_soil.md** - Documentation and references for understanding climate and soil data sources, including recommended databases and data formats.

## Key Features

- Download bioclimatic variables from WORLDCLIM2
- Extract environmental data for specific locations
- Format data for downstream GEA analyses
- Handle spatial data and environmental rasters

## Data Sources

This workflow utilizes publicly available environmental datasets including:
- **WORLDCLIM2** - Global bioclimatic variables (https://www.worldclim.org/)
- Soil databases (various regional and global sources)

## Usage

1. Review the `climate_soil.md` file for background information on available data sources
2. Open and adapt the `Get_bioclim_data_generic.Rmd` script for your study region
3. Modify geographic coordinates and data parameters as needed for your analysis
4. Run the script to download and process environmental variables

## Output

This workflow produces formatted environmental data files ready for use in subsequent GEA analyses.

## Related Sections

See also:
- [Mapping Accessions](../Mapping%20Accessions/README.md) - Geographic mapping of study samples
- [Cleaning Genomic Data](../Cleaning%20Genomic%20Data/README.md) - Preparing genomic data
