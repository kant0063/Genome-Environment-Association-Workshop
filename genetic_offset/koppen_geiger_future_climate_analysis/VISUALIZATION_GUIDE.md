# Visualization Guide - Köppen-Geiger Climate Change Analysis

This document describes all the plots produced by the analysis script.

## Overview

The script produces **12-13 high-resolution plots (300 DPI)** organized into two categories:
- **4 Spatial Maps** showing geographic distribution of climate classes
- **8-9 Statistical Plots** showing patterns and transitions

---

## PART 1: SPATIAL MAPS

These maps show the geographic distribution of Köppen-Geiger climate classes in your study region, with sample points overlaid.

### Map A: Historical Climate (`*_map_A_historical.png`)
**Type:** Spatial map with raster climate data and point overlay  
**Purpose:** Shows current (1991-2020) Köppen-Geiger climate classification  
**Features:**
- Full-color Köppen-Geiger classification using standard colors
- Study region automatically cropped to sample point extent + 10% buffer
- Sample points overlaid with distinct symbols/colors by population
- Black outline around points for visibility
- Complete legend for climate classes
- Background = light blue (ocean/missing data)

**Use case:** Understand the current climate context of your study area

---

### Map B: Future Climate (`*_map_B_future.png`)
**Type:** Spatial map with raster climate data and point overlay  
**Purpose:** Shows projected (2041-2070, SSP585) Köppen-Geiger climate classification  
**Features:**
- Same geographic extent as Map A for direct comparison
- Future climate classification colors
- Same sample points with population coding
- Allows side-by-side comparison with Map A

**Use case:** Visualize projected future climate conditions

---

### Map C: Climate Change Map (`*_map_change.png`)
**Type:** Binary change map  
**Purpose:** Highlights areas where climate class changed vs remained the same  
**Features:**
- **Red areas** = Climate class changed between time periods
- **Gray areas** = Climate class remained the same
- Sample points shown with different symbols:
  - **Triangles** = Samples that changed climate class
  - **Circles** = Samples with no climate change
- Clear visual of spatial patterns of change

**Use case:** Quickly identify where climate change is occurring geographically

---

### Map AB: Combined Comparison (`*_map_AB_combined.png`)
**Type:** Side-by-side faceted map  
**Purpose:** Direct visual comparison of historical vs future climate  
**Features:**
- Two panels: A) Historical | B) Future
- Same scale and extent for easy comparison
- Sample points shown as yellow dots with black outline
- Shared legend at bottom
- Publication-ready A/B figure format

**Use case:** Perfect for presentations and publications showing before/after

---

## PART 2: STATISTICAL PLOTS

### 1. Climate Change by Population (`*_by_population.png`)
**Type:** Horizontal bar chart  
**Purpose:** Shows the percentage of samples that changed climate class for each population  
**Features:**
- Populations sorted by percentage changed (highest to lowest)
- Percentage labels on each bar
- Quick overview of which populations are most affected

**Use case:** Identify which populations experience the most climate change

---

## 2. Stacked Bar - Changed vs Unchanged (`*_stacked_bar.png`)
**Type:** Stacked bar chart  
**Purpose:** Shows absolute counts of changed and unchanged samples per population  
**Features:**
- Red = Changed climate class
- Gray = Unchanged climate class
- Shows both relative and absolute impacts

**Use case:** Compare total sample sizes and change proportions across populations

---

## 3. Top 15 Climate Transitions (`*_top_transitions.png`)
**Type:** Horizontal bar chart  
**Purpose:** Shows the most common climate class transitions overall  
**Features:**
- Displays transitions like "Cfa → Csa" (detailed minor classes)
- Top 15 most common transitions
- Sample counts labeled on bars

**Use case:** Identify the most frequent types of climate change occurring

---

## 4. Major Climate Class Changes (`*_major_classes.png`)
**Type:** Horizontal bar chart  
**Purpose:** Shows percentage of samples changing major climate class by population  
**Features:**
- Major classes: A (Tropical), B (Arid), C (Temperate), D (Cold), E (Polar)
- Sorted by percentage changed
- More conservative metric than minor class changes

**Use case:** Assess fundamental climate zone shifts (e.g., Temperate → Arid)

---

## 5. Top 20 Minor Climate Transitions (`*_minor_transitions_top20.png`)
**Type:** Horizontal bar chart  
**Purpose:** Detailed view of the top 20 specific climate class transitions  
**Features:**
- Shows all 30 Köppen-Geiger sub-classes (Af, Am, Aw, BWh, etc.)
- Identifies precise climate transitions
- Larger plot to accommodate detailed labels

**Use case:** Understand specific climate changes (e.g., Cfa → Csa means "temperate no dry season hot summer" → "temperate dry summer hot summer")

---

## 6. Minor Transitions by Population (`*_minor_transitions_by_pop.png`)
**Type:** Faceted bar chart  
**Purpose:** Shows top 5 climate transitions for each population separately  
**Features:**
- Separate panel for each population
- Top 5 transitions per population
- Color-coded by population
- Identifies population-specific patterns

**Use case:** Compare which specific climate transitions affect different populations

---

## 7. Transition Matrix Heatmap (`*_transition_heatmap.png`)
**Type:** Heatmap  
**Purpose:** Complete matrix showing all historical → future class transitions  
**Features:**
- Rows = Historical climate classes (1991-2020)
- Columns = Future climate classes (2041-2070)
- Color intensity = number of samples
- Numbers in cells show exact counts
- **Only created if ≤15 unique classes** (to keep readable)

**Use case:** See complete transition patterns, identify cluster patterns in climate change

---

## 8. Climate Class Distribution (`*_class_distribution.png`)
**Type:** Grouped bar chart  
**Purpose:** Compare climate class frequencies between historical and future periods  
**Features:**
- Blue bars = Historical (1991-2020)
- Orange bars = Future (2041-2070)
- Percentage labels on bars
- Shows which classes are expanding vs contracting

**Use case:** Identify which climate classes are becoming more/less common

---

## 9. Transition Flow (`*_transition_flow.png`)
**Type:** Stacked bar chart  
**Purpose:** Shows where the top 8 historical climate classes transition to  
**Features:**
- X-axis = Historical climate classes (top 8 most common)
- Stacked bars show destination climate classes
- Uses viridis color palette for distinction
- Flow-like visualization showing transitions

**Use case:** Track where samples from each major historical climate class end up

---

## Map Specifications

**Spatial Maps (4 maps):**
- Automatically crop to study region (sample extent + 10% buffer)
- Köppen-Geiger standard color scheme
- Coordinate system: Geographic (WGS84)
- Resolution: Original raster resolution (~1km)
- Dimensions: 14" × 10"
- Sample points visible with population coding

**Statistical Plots (8-9 plots):**
- Dimensions vary by plot type (10-12" wide)
- Consistent color schemes
- Professional minimal theme

**All plots include:**
- Title and subtitle with time periods and scenario (SSP585)
- Axis labels with units
- High resolution (300 DPI) for publication
- Consistent color schemes
- Professional theme (minimal, clean)

**File format:** PNG  
**Resolution:** 300 DPI  

---

## Köppen-Geiger Color Scheme (Standard)

The spatial maps use the official Köppen-Geiger color scheme:

**Tropical (A) - Red tones:**
- Af (Rainforest) = Dark red
- Am (Monsoon) = Bright red
- Aw (Savannah) = Light red/pink

**Arid (B) - Yellow/Tan tones:**
- BWh (Hot desert) = Bright yellow
- BWk (Cold desert) = Yellow-tan
- BSh (Hot steppe) = Orange-tan
- BSk (Cold steppe) = Light tan

**Temperate (C) - Green/Yellow-green tones:**
- Cs* (Dry summer) = Yellow to olive
- Cw* (Dry winter) = Green tones
- Cf* (No dry season) = Bright green

**Cold (D) - Blue/Cyan/Magenta tones:**
- Ds* (Dry summer) = Magenta/purple
- Dw* (Dry winter) = Blue tones
- Df* (No dry season) = Cyan tones

**Polar (E) - Gray tones:**
- ET (Tundra) = Light gray
- EF (Ice cap) = Dark gray

---

## Color Schemes for Statistical Plots

**Consistent colors used:**
- **Changed** = Red (#e74c3c)
- **Unchanged** = Gray (#95a5a6)
- **Historical** = Blue (#3498db)
- **Future** = Orange (#e67e22)
- **Population bars** = Steel blue
- **Major classes** = Orange (#e67e22)
- **Minor transitions** = Purple (#9b59b6)
- **Transition flow** = Viridis (rainbow spectrum)

---

## Usage Tips

**For presentations:**
- Use Map AB (combined) for overview of study area
- Use Map C (change map) to highlight spatial patterns
- Use statistical plots 1, 2, and 8 for summary statistics
- Use plot 9 (flow) for storytelling

**For publications:**
- Map AB is publication-ready as Figure 1 (A/B panels)
- Maps A and B can be separate figures if needed
- Statistical plots 3, 5, and 7 provide most detail
- Plot 8 shows before/after comparison
- All plots are publication-ready at 300 DPI

**For reports:**
- Map A/B combined for geographic context
- Map C to show change hotspots
- Statistical plot 1 for executive summary
- Plot 4 for major climate zone shifts
- Plot 6 for population-specific details

**For posters:**
- Map AB works well as large background/context figure
- Use Map C to draw attention to change areas
- Combine with 2-3 key statistical plots

---

## Interpretation Guide

### Minor vs Major Class Changes

**Minor class changes** (30 classes):
- More sensitive to climate shifts
- Example: Cfa → Csa (change in precipitation seasonality)
- Tracks subtle climate modifications

**Major class changes** (5 classes):
- Fundamental climate zone shifts
- Example: C → B (Temperate → Arid)
- More dramatic transformations

### Reading Transitions

Transition notation: `Historical → Future`

**Examples:**
- `Cfa → Csa` = Lost year-round precipitation, gained dry summer
- `Cfb → Cfa` = Warmer summers (warm → hot)
- `BSk → BWk` = Steppe became desert (both cold, arid)
- `C → B` = Temperate became arid (major shift)

### Climate Class Codes

**First letter (Major class):**
- **A** = Tropical (hot year-round)
- **B** = Arid (dry)
- **C** = Temperate (mild)
- **D** = Cold (severe winters)
- **E** = Polar (very cold)

**Second letter (Precipitation):**
- **f** = No dry season
- **s** = Dry summer
- **w** = Dry winter
- **W** = Desert
- **S** = Steppe
- **T** = Tundra
- **F** = Frost

**Third letter (Temperature):**
- **a** = Hot summer
- **b** = Warm summer
- **c** = Cold summer
- **d** = Very cold winter
- **h** = Hot (for deserts)
- **k** = Cold (for deserts)

---

## Customization

You can modify plot appearance by editing the script:

**Change colors:**
```r
scale_fill_manual(values = c("Changed" = "#YOUR_COLOR"))
```

**Adjust dimensions:**
```r
ggsave(..., width = 12, height = 8)
```

**Change number of top items:**
```r
head(20)  # Change to show more/fewer
```

**Font sizes:**
```r
theme(axis.text = element_text(size = 12))
```

---

## File Naming

All plots use the prefix you specify in the script:

```r
output_prefix <- "koppen_geiger_analysis"
```

Results in files like:
- `koppen_geiger_analysis_by_population.png`
- `koppen_geiger_analysis_transition_heatmap.png`
- etc.

Change the prefix to organize multiple analyses:
```r
output_prefix <- "lentil_mediterranean_analysis"
output_prefix <- "pepper_tropics_analysis"
```
