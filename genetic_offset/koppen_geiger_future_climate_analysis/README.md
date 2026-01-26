# Köppen-Geiger Climate Change Analysis

This R script calculates changes in Köppen-Geiger climate classifications between two time periods (1991-2020 to 2041-2070 under SSP585) for sample locations provided in a CSV file.

## Requirements

### R Packages
```r
install.packages(c("terra", "dplyr", "tidyr", "ggplot2", "readr", "viridis"))
```

### Data Files

**1. Köppen-Geiger Raster Data**
- Download from: https://figshare.com/articles/dataset/21789074
- File: `koppen_geiger_tif.zip` (125 MB)
- Extract the zip file to your working directory
- Reference: Beck et al. (2023) Scientific Data 10, 724

**2. Your Sample Data CSV**
- Must contain at minimum: longitude, latitude, and population columns
- See `sample_input_template.csv` for format example

## Input CSV Format

Your CSV file should have the following structure:

```csv
sample_id,longitude,latitude,population
Sample001,-122.4194,37.7749,Pop1
Sample002,-118.2437,34.0522,Pop1
Sample003,-87.6298,41.8781,Pop2
...
```

### Required Columns:
- **longitude**: Decimal degrees (-180 to 180)
- **latitude**: Decimal degrees (-90 to 90)
- **population**: Grouping variable (can be population, site, region, etc.)

### Optional Columns:
- **sample_id**: Unique identifier for each sample (will be auto-generated if not provided)
- Any other columns in your CSV will be preserved in the output

## Configuration

Edit these variables at the top of `koppen_geiger_climate_change_analysis.R`:

```r
# 1. Path to your input CSV file
input_csv <- "your_samples.csv"

# 2. Column names (adjust to match your CSV)
lon_column <- "longitude"
lat_column <- "latitude"
population_column <- "population"
id_column <- "sample_id"  # Set to NULL if you don't have this column

# 3. Paths to Köppen-Geiger rasters
kg_historical_path <- "koppen_geiger/1991_2020/koppen_geiger_0p00833333.tif"
kg_future_path <- "koppen_geiger/2041_2070/ssp585/koppen_geiger_0p00833333.tif"

# 4. Output file prefix
output_prefix <- "koppen_geiger_analysis"
```

## Usage

1. **Prepare your data:**
   - Create a CSV file with your sample locations
   - Ensure it has longitude, latitude, and population columns

2. **Download Köppen-Geiger data:**
   ```bash
   # Download koppen_geiger_tif.zip from Figshare
   # Extract it to your working directory
   ```

3. **Configure the script:**
   - Edit the configuration section with your file paths and column names

4. **Run the analysis:**
   ```r
   source("koppen_geiger_climate_change_analysis.R")
   ```

## Outputs

The script produces multiple output files:

### CSV Files:
1. **`*_full_results.csv`** - Complete dataset with all original columns plus:
   - `kg_historical_code` - Historical Köppen-Geiger numeric code
   - `kg_historical_class` - Historical climate class (e.g., "Cfa")
   - `kg_historical_description` - Human-readable description
   - `kg_historical_major` - Major climate class (A/B/C/D/E)
   - `kg_future_code` - Future Köppen-Geiger numeric code
   - `kg_future_class` - Future climate class
   - `kg_future_description` - Future description
   - `kg_future_major` - Future major class
   - `climate_changed` - Boolean: did climate class change?
   - `major_class_changed` - Boolean: did major class change?
   - `change_type` - Description of transition (e.g., "Cfa → Csa")
   - `major_change_type` - Major class transition

2. **`*_population_summary.csv`** - Summary statistics by population
   - Number of samples per population
   - Number and percentage that changed climate class
   - Number and percentage that changed major class

3. **`*_change_details.csv`** - Detailed breakdown of climate transitions by population

4. **`*_major_class_summary.csv`** - Summary of major climate class changes

5. **`*_overall_summary.csv`** - Overall statistics across all samples

6. **`*_top_transitions.csv`** - Most common climate transitions

### Plots (PNG, 300 DPI):

**Spatial Maps (4 maps):**
1. **`*_map_A_historical.png`** - Historical climate map (1991-2020) with sample points overlaid
2. **`*_map_B_future.png`** - Future climate map (2041-2070, SSP585) with sample points overlaid
3. **`*_map_change.png`** - Climate change map showing areas that changed (red) vs unchanged (gray)
4. **`*_map_AB_combined.png`** - Side-by-side comparison of historical and future climate (A/B figure)

**Statistical Plots (8-9 plots):**
5. **`*_by_population.png`** - Bar chart showing percentage of samples with climate change by population
6. **`*_stacked_bar.png`** - Stacked bar chart showing changed vs unchanged samples per population
7. **`*_top_transitions.png`** - Top 15 most common climate class transitions (major classes)
8. **`*_major_classes.png`** - Major climate class changes (A/B/C/D/E) by population
9. **`*_minor_transitions_top20.png`** - Top 20 detailed minor climate class transitions (e.g., Cfa → Csa)
10. **`*_minor_transitions_by_pop.png`** - Top 5 minor climate transitions for each population (faceted)
11. **`*_transition_heatmap.png`** - Transition matrix heatmap showing all historical → future class changes (only created if ≤15 unique classes)
12. **`*_class_distribution.png`** - Side-by-side comparison of climate class distribution in historical vs future periods
13. **`*_transition_flow.png`** - Stacked bar chart showing where the top 8 historical climate classes transition to

**Total: 6 CSV files + 12-13 plots (depending on class diversity)**

## Köppen-Geiger Classification

### Major Climate Classes:
- **A** - Tropical
- **B** - Arid
- **C** - Temperate
- **D** - Cold
- **E** - Polar

### Sub-classes (examples):
- **Af** - Tropical rainforest
- **Aw** - Tropical savannah
- **BSk** - Arid steppe, cold
- **Cfa** - Temperate, no dry season, hot summer
- **Dfb** - Cold, no dry season, warm summer
- **ET** - Polar tundra

See `kg_legend` in the script for complete classification system.

## Time Periods

- **Historical**: 1991-2020 (baseline, observational data)
- **Future**: 2041-2070 under SSP5-8.5 scenario

**Note:** The original request was for 1991-2000, but the Beck et al. (2023) dataset's closest baseline period is 1991-2020. This is actually better as it provides a more robust 30-year climatological baseline.

## Example Workflow

```r
# 1. Set working directory
setwd("/path/to/your/project")

# 2. Load packages (install if needed)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# 3. Run the analysis
source("koppen_geiger_climate_change_analysis.R")
```

## Troubleshooting

### Error: "Input CSV file not found"
- Check that `input_csv` path is correct
- Use absolute path if relative path doesn't work

### Error: "Missing required columns in CSV"
- Verify your CSV has the columns specified in the configuration
- Check for typos in column names (case-sensitive!)

### Error: "Köppen-Geiger raster not found"
- Ensure you've downloaded and extracted `koppen_geiger_tif.zip`
- Check that paths in configuration match your directory structure

### Warning: "Longitude/latitude values outside valid range"
- Check your coordinate data for errors
- Ensure longitude is -180 to 180, latitude is -90 to 90

### Warning: "Missing coordinates"
- Some rows have NA values in longitude or latitude
- These rows will be automatically removed from analysis

## Citation

If you use this analysis in publications, please cite:

Beck, H.E., T.R. McVicar, N. Vergopolan, A. Berg, N.J. Lutsko, A. Dufour, Z. Zeng, X. Jiang, A.I.J.M. van Dijk, D.G. Miralles (2023). High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections. Scientific Data 10, 724. https://doi.org/10.1038/s41597-023-02549-6

## Additional Notes

- The script preserves all original columns from your input CSV
- Missing coordinate values are automatically removed with a warning
- The script validates input data and provides helpful error messages
- All plots are saved as high-resolution PNG files (300 DPI)
- Results are sorted by percentage of climate change for easy interpretation

## Support

For questions about:
- **Köppen-Geiger data**: See Beck et al. (2023) documentation
- **This script**: Check the inline comments in the R script
- **Your specific data**: Validate your CSV format matches the template

## Version History

- v1.0 (2025-01-25): Initial release
  - CSV input support
  - Multiple output formats
  - Comprehensive visualizations
  - Detailed documentation
