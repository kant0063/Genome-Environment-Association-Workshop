# Mapping Accessions

This directory contains workflows for geographic mapping and visualization of study samples/accessions, including tools for working with geographic coordinates, creating maps, and analyzing spatial distributions.

## Overview

Geographic information is fundamental to genome-environment association studies. This section provides tools for:
- Visualizing spatial distribution of study samples
- Mapping environmental gradients
- Assessing geographic bias and sampling coverage
- Creating publication-quality geographic maps

## Contents

### Documentation

- **Mapping Your Accessions lesson plan.md** - Educational lesson plan for:
  - Geographic data concepts
  - Mapping best practices
  - Interpretation guidelines

- **Geographic_Sample_Mapping_README.md** - Comprehensive guide including:
  - Data requirements and formats
  - Mapping methodology
  - Visualization options
  - Troubleshooting

- **geolocation databases.md** - Reference guide for:
  - Available geolocation databases
  - Data source descriptions
  - Accuracy and coverage information
  - Citation requirements

### Analysis Scripts

- **Geographic_Sample_Mapping.Rmd** - R Markdown workflow for:
  - Loading geographic and environmental data
  - Creating base maps
  - Plotting sample locations
  - Adding environmental layers
  - Customizing map aesthetics
  - Exporting publication-ready figures

## Key Features

- Interactive and static map creation
- Multiple mapping libraries and approaches
- Environmental overlay visualization
- Spatial analysis and clustering
- Publication-quality figure generation
- Multiple geographic projection support

## Workflow Steps

1. Prepare geographic coordinates
2. Load and validate sample metadata
3. Obtain environmental data (climate rasters, etc.)
4. Create base map
5. Plot sample locations
6. Add environmental information (coloring, contours)
7. Customize appearance and add map elements
8. Export for publication

## Input Requirements

- Sample geographic coordinates (latitude/longitude)
- Sample metadata and identifiers
- Environmental data (optional, for layering)
- Basemap data or tiles
- Geolocation database references

## Output

Map outputs include:
- Static map images (PNG, PDF, SVG formats)
- Interactive maps (HTML)
- Map legends and scale bars
- Spatial statistics and summaries

## Geolocation Data Sources

Common resources include:
- **Google Maps API** - Reverse geocoding
- **OpenStreetMap** - Geographic data and basemaps
- **GADM** - Administrative boundaries
- **Natural Earth** - Coastlines and geographic features
- **Species locality databases** - Plant collections, herbaria

## Spatial Considerations

- **Coordinate precision** - Appropriate decimal places for latitude/longitude
- **Projections** - Choosing appropriate map projections
- **Accuracy assessment** - Evaluating georeferencing quality
- **Sampling bias** - Recognizing geographic clustering or gaps

## Integrating with GEA Analysis

Maps created here help:
- Visualize spatial extent of environmental gradients
- Identify associations between geographic distribution and environmental variables
- Assess geographic bias in sampling
- Plan targeted sampling for coverage
- Interpret population structure in geographic context

## Related Sections

See also:
- [Accessing Climate and Soil Data](../Accessing%20Climate%20and%20Soil%20Data/README.md) - Obtain environmental data for overlays
- [Simple Population Structure](../Simple_Population_structure/README.md) - Integrate population structure with geography
- [Conducting the Environmental Association Analysis](../Conducting%20the%20Environmental%20Association%20Analysis/README.md) - Use geographic context for interpretation
