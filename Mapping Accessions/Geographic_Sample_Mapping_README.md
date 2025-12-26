# Geographic Sample Mapping with Regional Zoom

A comprehensive R workflow for creating publication-quality maps showing the geographic distribution of samples, with automatic world overview and regional detail maps.

## Overview

This workflow creates professional maps for visualizing sample collection locations using modern R spatial packages. It automatically generates both broad overview maps and detailed regional zoom maps based on your data extent.

### What This Analysis Does

1. Reads sample coordinates from CSV
2. Validates and filters coordinate data
3. Creates world/continental overview map
4. Generates regional zoom map with data extent
5. Produces combined side-by-side visualization
6. Optional: Labeled maps and interactive maps
7. Exports high-resolution publication-ready figures

### Key Features

✓ **Automatic extent calculation**: No manual coordinate input needed  
✓ **Dual scale visualization**: Overview + regional detail  
✓ **High-resolution maps**: Country borders and coastlines  
✓ **Customizable styling**: Colors, sizes, themes  
✓ **Scale bars and north arrows**: Publication standards  
✓ **Group visualization**: Color by population/region  
✓ **Interactive maps**: Optional leaflet integration  

---

## Quick Start

### Minimal Example

```r
# Prepare your CSV with columns: SampleID, Latitude, Longitude
# Run analysis
rmarkdown::render("Geographic_Sample_Mapping.Rmd")
```

### Output

- World overview map (PNG, 300 dpi)
- Regional detail map (PNG, 300 dpi)
- Combined side-by-side map (PNG, 300 dpi)
- Optional: Labeled map, interactive HTML map

---

## Installation

### Required R Packages

```r
# Core packages
install.packages(c("tidyverse", "sf", "rnaturalearth", 
                   "rnaturalearthdata", "ggspatial", "patchwork",
                   "knitr", "kableExtra"))

# High-resolution map data (recommended)
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")

# Optional: Interactive maps
install.packages(c("leaflet", "htmlwidgets"))
```

### System Requirements

- **R version**: ≥ 4.0.0
- **RAM**: 4GB minimum
- **Storage**: ~500MB for map data

---

## Input Data Requirements

### CSV File Format

**Required columns:**
- Sample ID (any name)
- Latitude (decimal degrees, -90 to 90)
- Longitude (decimal degrees, -180 to 180)

**Optional columns:**
- Population/Region/Group (for color coding)
- Site names
- Collection dates
- Any other metadata

### Example CSV

```csv
SampleID,Latitude,Longitude,Population,Region
Sample001,45.523,-122.677,PopA,Northwest
Sample002,46.872,-121.769,PopA,Northwest
Sample003,42.360,-71.059,PopB,Northeast
Sample004,37.775,-122.419,PopC,California
```

### Coordinate Format

**IMPORTANT:**
- Use **decimal degrees** (not degrees/minutes/seconds)
- Use **WGS84** coordinate system
- Negative values:
  - Latitude: South = negative
  - Longitude: West = negative

**Convert from degrees/minutes/seconds:**
```r
# Example: 45°30'15"N, 122°40'30"W
lat_dd <- 45 + 30/60 + 15/3600  # = 45.50417
lon_dd <- -(122 + 40/60 + 30/3600)  # = -122.675
```

---

## Usage

### Basic Workflow

1. Prepare CSV with sample coordinates
2. Edit parameters in R Markdown
3. Run analysis
4. Examine and export maps

### Key Parameters

```r
# ==================== INPUT ====================
DATA_CSV <- "your_samples.csv"

# ==================== COLUMNS ====================
SAMPLE_ID_COLUMN <- "SampleID"
LATITUDE_COLUMN <- "Latitude"
LONGITUDE_COLUMN <- "Longitude"
GROUP_COLUMN <- "Population"        # Optional

# ==================== MAP EXTENT ====================
WORLD_XLIM <- c(-180, 180)          # NULL for auto
WORLD_YLIM <- c(-90, 90)            # NULL for auto
REGIONAL_BUFFER <- 2                # Degrees around data

# ==================== VISUAL ====================
POINT_SIZE <- 3
POINT_ALPHA <- 0.8
COLOR_PALETTE <- "Set1"             # RColorBrewer palette
MAP_FILL_COLOR <- "wheat"
OCEAN_COLOR <- "lightblue"
```

### Parameter Guide

**World Map Extent:**
- Set to `NULL` for automatic zoom to your data region
- Use `c(-180, 180)` and `c(-90, 90)` for full world view
- Customize for specific regions (e.g., Europe, Asia)

**Regional Buffer:**
- Degrees to add around data extent
- `2` = good default for most datasets
- Increase for more context (e.g., `5` or `10`)
- Decrease for tighter zoom (e.g., `1` or `0.5`)

**Color Palettes:**

Available RColorBrewer palettes:
- Qualitative: `"Set1"`, `"Set2"`, `"Set3"`, `"Dark2"`, `"Paired"`
- Sequential: `"Blues"`, `"Greens"`, `"Reds"`, `"Purples"`

**Point Size:**
- `1-2`: Very small (for dense datasets)
- `3-4`: Medium (default, good for most)
- `5-8`: Large (for sparse datasets or emphasis)

---

## Output Files

### File Structure

```
maps/
├── sample_map_world_overview.png       # World/continental view
├── sample_map_regional.png             # Regional zoom view
├── sample_map_combined.png             # Side-by-side comparison
├── sample_map_labeled.png              # With sample labels (if ≤50 samples)
├── sample_map_interactive.html         # Interactive leaflet map (optional)
├── sample_map_coordinates.csv          # Valid coordinate data
└── sample_map_summary.csv              # Summary statistics
```

### Map Descriptions

| File | Description | Dimensions |
|------|-------------|------------|
| `*_world_overview.png` | Broad view with regional extent box | 12" × 8" |
| `*_regional.png` | Zoomed view of data extent | 12" × 8" |
| `*_combined.png` | Both maps side-by-side | 16" × 8" |
| `*_labeled.png` | Regional map with sample IDs | 12" × 8" |
| `*_interactive.html` | Interactive web map | - |

All PNG files at 300 DPI (publication quality)

---

## Interpretation Guide

### Map Elements

**Scale Bar:**
- Shows actual distance
- Automatically sized for map extent
- Located bottom-left

**North Arrow:**
- Shows true north direction
- Fancy orienteering style
- Located top-right

**Red Dashed Box (on world map):**
- Shows extent of regional zoom map
- Helps orient between scales

### When to Use Each Map

**World Overview:**
- Show broad geographic context
- Display multi-continental sampling
- Presentations and proposals

**Regional Detail:**
- Show fine-scale distribution
- Emphasize local patterns
- Detailed analysis and publications

**Combined Map:**
- Best of both worlds
- Good for publications
- Shows context + detail

**Labeled Map:**
- Useful for ≤50 samples
- Reference specific samples
- Field site identification
- Only created if sample count allows

### Color Coding

**By Population/Region:**
- Different colors = different groups
- Legend shows group names
- Useful for comparing distributions

**Single Color:**
- All samples same color
- Simpler visualization
- Focus on geographic pattern

---

## Advanced Options

### Custom Map Extents

```r
# Europe only
WORLD_XLIM <- c(-10, 40)
WORLD_YLIM <- c(35, 70)

# North America
WORLD_XLIM <- c(-130, -60)
WORLD_YLIM <- c(25, 50)

# Southeast Asia
WORLD_XLIM <- c(90, 140)
WORLD_YLIM <- c(-10, 25)
```

### Multiple Regional Zoom Maps

```r
# Define specific regions of interest
region1_xlim <- c(lon_min1, lon_max1)
region1_ylim <- c(lat_min1, lat_max1)

# Create additional maps
region1_map <- ggplot() + 
  geom_sf(data = world_hires, ...) +
  coord_sf(xlim = region1_xlim, ylim = region1_ylim) +
  # ... rest of plot code
```

### Custom Colors

```r
# Manual color specification
library(RColorBrewer)
my_colors <- brewer.pal(5, "Set1")

# Or define your own
my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
```

### Add Additional Annotations

```r
# Add text labels
+ annotate("text", x = -120, y = 40, 
           label = "Region Name", size = 5)

# Add specific points
+ geom_point(aes(x = -122, y = 37.8), 
             color = "blue", size = 8, shape = 17)

# Add lines
+ geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "red", size = 1)
```

### Save in Different Formats

```r
# Save as PDF (vector format)
ggsave("map.pdf", plot = world_map, 
       width = 12, height = 8, device = "pdf")

# Save as SVG (vector format)
ggsave("map.svg", plot = world_map, 
       width = 12, height = 8, device = "svg")

# Save as TIFF (high quality raster)
ggsave("map.tiff", plot = world_map, 
       width = 12, height = 8, dpi = 600, device = "tiff")
```

---

## Troubleshooting

### Common Issues

#### 1. "Missing required columns"

**Problem:** Column names don't match parameters

**Solution:**
```r
# Check your column names
colnames(sample_data)

# Update parameters to match
SAMPLE_ID_COLUMN <- "ID"  # or whatever your column is named
LATITUDE_COLUMN <- "Lat"
LONGITUDE_COLUMN <- "Long"
```

#### 2. "All coordinates invalid/NA"

**Problem:** Coordinates not recognized

**Common causes:**
- Degrees/minutes/seconds format (not decimal)
- Comma instead of decimal point
- Missing values coded as text ("NA", "missing")

**Solutions:**
```r
# Convert DMS to decimal degrees
# Check for text values
table(sample_data$Latitude, useNA = "always")

# Replace text NA with actual NA
sample_data$Latitude[sample_data$Latitude == "NA"] <- NA
```

#### 3. "High-resolution maps not loading"

**Problem:** rnaturalearthhires not installed

**Solution:**
```r
# Install from r-universe
install.packages("rnaturalearthhires",
                 repos = "http://packages.ropensci.org",
                 type = "source")

# Or use medium resolution (already included)
world_hires <- world  # Falls back to medium resolution
```

#### 4. Maps look empty/blank

**Problem:** Coordinate system mismatch

**Possible causes:**
- Latitude and longitude swapped
- Wrong hemisphere (positive instead of negative)
- Coordinates far outside expected range

**Solutions:**
```r
# Check coordinate ranges
range(sample_data$Latitude, na.rm = TRUE)   # Should be -90 to 90
range(sample_data$Longitude, na.rm = TRUE)  # Should be -180 to 180

# Swap if needed
temp <- sample_data$Latitude
sample_data$Latitude <- sample_data$Longitude
sample_data$Longitude <- temp

# Fix hemisphere
sample_data$Longitude <- -sample_data$Longitude  # If in Western hemisphere
```

#### 5. Points don't appear on map

**Problem:** CRS mismatch or plotting order

**Solutions:**
```r
# Ensure WGS84 coordinates
# Check plotting order (points should come after map)
# Verify coordinate values are within map extent
print(range(sample_valid$Longitude))
print(range(sample_valid$Latitude))
```

#### 6. Labels overlap/unreadable

**Problem:** Too many samples or labels too large

**Solutions:**
```r
# Reduce label size
size = 2  # in geom_text

# Use check_overlap
check_overlap = TRUE

# Only label subset
sample_subset <- sample_valid %>% slice(1:20)
```

---

## Best Practices

### 1. Data Quality

**Before mapping:**
- Validate all coordinates in WGS84
- Check for outliers/errors
- Remove duplicates
- Verify hemisphere signs

### 2. Map Design

**For publications:**
- Use 300+ DPI resolution
- Include scale bar and north arrow
- Clear, readable labels
- Appropriate color scheme
- Consistent styling across figures

### 3. Resolution Selection

**When to use high-resolution:**
- Regional maps (country/state level)
- Small study areas
- Need detailed coastlines/borders

**When medium resolution is fine:**
- World/continental maps
- Large study areas
- Overview figures

### 4. Color Selection

**Best practices:**
- Maximum 7-8 colors for groups
- Use colorblind-friendly palettes
- Avoid red-green combinations
- Check contrast on printed versions

**Colorblind-safe palettes:**
```r
# RColorBrewer
COLOR_PALETTE <- "Set2"      # Good
COLOR_PALETTE <- "Dark2"     # Good
COLOR_PALETTE <- "Paired"    # Good

# Avoid
COLOR_PALETTE <- "Set1"      # Contains red-green
```

### 5. File Formats

**For different uses:**

| Format | Use Case | Command |
|--------|----------|---------|
| PNG | Web, presentations | Default |
| PDF | Publications, vector | `device = "pdf"` |
| TIFF | High-quality print | `device = "tiff"` |
| SVG | Editing in Illustrator | `device = "svg"` |

---

## Examples

### Example 1: Simple Map

```r
# Just show points, no grouping
DATA_CSV <- "samples.csv"
GROUP_COLUMN <- NULL  # No color coding
```

### Example 2: Regional Study

```r
# Focus on specific region
WORLD_XLIM <- c(10, 20)   # Italy
WORLD_YLIM <- c(36, 47)
REGIONAL_BUFFER <- 1      # Tight zoom
```

### Example 3: Multi-Continental

```r
# Show broad distribution
WORLD_XLIM <- c(-180, 180)  # Full world
WORLD_YLIM <- c(-90, 90)
REGIONAL_BUFFER <- 10       # Large buffer
```

### Example 4: Dense Sampling

```r
# Many samples in small area
POINT_SIZE <- 1.5           # Smaller points
POINT_ALPHA <- 0.6          # More transparent
REGIONAL_BUFFER <- 0.5      # Tight focus
```

---

## Citation

If you use this workflow in publications, please cite:

**For rnaturalearth:**
- South, A. (2017). rnaturalearth: World Map Data from Natural Earth. R package.

**For sf:**
- Pebesma, E. (2018). Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal*, 10(1), 439-446.

**For ggspatial:**
- Dunnington, D. (2021). ggspatial: Spatial Data Framework for ggplot2. R package.

**Example citation:**
> "Sample locations were visualized using R (R Core Team, 2023) with the sf (Pebesma, 2018), rnaturalearth (South, 2017), and ggspatial (Dunnington, 2021) packages. Maps show both continental overview and regional detail at 300 DPI resolution."

---

## Additional Resources

### Tutorials

- sf package: https://r-spatial.github.io/sf/
- rnaturalearth: https://docs.ropensci.org/rnaturalearth/
- ggspatial: https://paleolimbot.github.io/ggspatial/

### Related Packages

- **tmap**: Thematic maps
- **mapview**: Interactive viewing
- **ggmap**: Google Maps integration
- **leaflet**: Interactive web maps

### Map Data Sources

- Natural Earth: https://www.naturalearthdata.com/
- GADM: https://gadm.org/ (administrative boundaries)
- OpenStreetMap: https://www.openstreetmap.org/

---

## License

This workflow is provided under the MIT License.

## Contributing

Suggestions and improvements welcome!

---

**Last updated:** December 2024

**Version:** 1.0

**Tested with:** R 4.3.0, sf 1.0-14, rnaturalearth 0.3.4
