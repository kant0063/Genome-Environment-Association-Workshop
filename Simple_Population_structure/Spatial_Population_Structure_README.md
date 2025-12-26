# Spatial Population Structure Analysis with TESS3

A comprehensive R workflow for analyzing spatial genetic structure using TESS3 (Terrain-based Estimation of Spatial Structure), which incorporates geographic information directly into ancestry estimation.

## Overview

**TESS3** is a spatially-explicit method for inferring population structure that uses geographic coordinates to model continuous variation in ancestry across landscapes. Unlike STRUCTURE or ADMIXTURE, TESS3 doesn't require discrete population assignments and produces beautiful spatial interpolation maps.

### What This Analysis Does

1. Reads genotype data from VCF files
2. Matches samples with geographic coordinates
3. Runs TESS3 for multiple K values
4. Determines optimal K via cross-validation
5. Generates STRUCTURE-like ancestry bar plots
6. Creates spatial interpolation maps
7. Exports ancestry coefficients (Q-matrices)

### Key Applications

- **Landscape genetics**: Understand gene flow patterns across space
- **Conservation**: Identify genetic management units
- **Biogeography**: Map genetic clines and transitions
- **Invasion biology**: Track introduction routes and spread
- **Agriculture**: Characterize germplasm geographic structure
- **Climate adaptation**: Identify locally adapted populations

### Advantages Over STRUCTURE/ADMIXTURE

✓ **Spatially explicit**: Directly uses geographic coordinates  
✓ **Continuous variation**: Models smooth transitions, not discrete pops  
✓ **Beautiful maps**: Generates spatial interpolation visualizations  
✓ **Faster**: More efficient computation  
✓ **No priors**: Doesn't require population assignments

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Data Requirements](#input-data-requirements)
- [Usage](#usage)
- [Output Files](#output-files)
- [Interpretation Guide](#interpretation-guide)
- [Advanced Options](#advanced-options)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Installation

### Required R Packages

```r
# Core packages
install.packages(c("vcfR", "maps", "rworldmap", "ggplot2", "dplyr", 
                   "tidyr", "patchwork"))

# TESS3 package
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("bcm-uga/TESS3_encho_sen")

# Alternative installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("tess3r")
```

### System Requirements

- **R version**: ≥ 4.0.0
- **RAM**: 8GB minimum (16GB+ recommended)
- **Storage**: ~2GB free space
- **CPU**: Multi-core recommended (parallel processing)

---

## Quick Start

### Minimal Example

```r
# Set your file paths
VCF_FILE <- "my_genotypes.vcf"
COORD_CSV <- "sample_coordinates.csv"

# Run analysis
rmarkdown::render("Spatial_Population_Structure.Rmd")
```

### Expected Runtime

- **Small dataset** (50 samples, 10K SNPs): ~5-10 minutes
- **Medium dataset** (100 samples, 50K SNPs): ~20-40 minutes
- **Large dataset** (200+ samples, 100K+ SNPs): ~1-3 hours

*Runtime depends on number of K values tested and CPU cores*

---

## Input Data Requirements

### 1. VCF File (Genotype Data)

**Format:** Standard VCF format (`.vcf` or `.vcf.gz`)

**Requirements:**
- SNP genotype calls
- Quality-filtered
- Biallelic SNPs only

**Example VCF:**
```
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
Chr1    1000    .       A    G    99    PASS    .     GT      0/0      0/1
Chr1    2000    .       C    T    99    PASS    .     GT      0/1      1/1
```

### 2. Coordinate CSV

**Required columns:**
- **Sample IDs**: Must match VCF exactly
- **Latitude**: Decimal degrees (WGS84)
- **Longitude**: Decimal degrees (WGS84)

**Example CSV:**
```csv
SampleID,Latitude,Longitude,Population,Site
Sample1,45.523,-122.677,PopA,Site1
Sample2,46.872,-121.769,PopA,Site2
Sample3,42.360,-71.059,PopB,Site3
```

**Important:**
- Use decimal degrees (not degrees/minutes/seconds)
- Latitude: -90 to 90 (negative = South)
- Longitude: -180 to 180 (negative = West)
- Missing coordinates not allowed
- WGS84 coordinate system

### Recommended Sampling Design

**For best results:**
- **Minimum 30 samples** (50+ preferred)
- Samples distributed across geographic range
- Multiple samples per location helpful
- Avoid extreme sampling gaps
- Cover environmental gradients

---

## Usage

### Basic Workflow

1. **Prepare data** (see Input Data Requirements)
2. **Edit parameters** at top of R Markdown
3. **Run analysis**
4. **Examine outputs** and interpret

### Key Parameters

```r
# ==================== INPUT FILES ====================
VCF_FILE <- "your_genotypes.vcf"
COORD_CSV <- "your_coordinates.csv"

# ==================== COLUMN NAMES ====================
SAMPLE_ID_COLUMN <- "SampleID"
LATITUDE_COLUMN <- "Latitude"
LONGITUDE_COLUMN <- "Longitude"

# ==================== TESS3 PARAMETERS ====================
K_VALUES <- 1:7                # Range of K to test
PLOIDY <- 2                    # 2 for diploid organisms
N_CORES <- 4                   # CPU cores to use

# ==================== MAP EXTENT ====================
MAP_XLIM <- c(-180, 180)       # Longitude range
MAP_YLIM <- c(-90, 90)         # Latitude range
```

### Parameter Selection Guide

**K Values:**
- Test range from 1 to 10 (or more)
- Optimal K determined by cross-validation
- Lower K = broader patterns
- Higher K = finer details
- May test multiple ranges

**Ploidy:**
- 2 for most diploid organisms
- 4 for autotetraploids
- 6 for autohexaploids

**CPU Cores:**
- Use all but 1-2 cores
- More cores = faster computation
- Check: `parallel::detectCores()`

**Map Extent:**
- Adjust based on your data range
- Add buffer around sample locations
- Use appropriate projection if needed

### Running from Command Line

```bash
# Generate HTML report
Rscript -e "rmarkdown::render('Spatial_Population_Structure.Rmd')"

# Generate PDF (requires LaTeX)
Rscript -e "rmarkdown::render('Spatial_Population_Structure.Rmd', output_format='pdf_document')"
```

---

## Output Files

### Directory Structure

```
tess3_results/
├── tess3_cross_validation.csv          # CV scores for all K values
├── tess3_summary.csv                   # Analysis summary
├── tess3_barplots_combined.png         # Ancestry bar plots (all K)
├── tess3_spatial_maps_baseR.png        # Spatial maps (base R)
├── tess3_spatial_maps_ggplot.png       # Spatial maps (ggplot2)
├── tess3_qmatrix_K2.csv                # Q-matrix for K=2
├── tess3_qmatrix_K3.csv                # Q-matrix for K=3
└── ...                                 # Q-matrices for other K values
```

### File Descriptions

| File | Description |
|------|-------------|
| `*_cross_validation.csv` | Cross-validation scores for each K value |
| `*_summary.csv` | Sample size, SNP count, optimal K |
| `*_barplots_combined.png` | Multi-panel STRUCTURE-like plots |
| `*_spatial_maps_baseR.png` | Spatial interpolation maps (base R) |
| `*_spatial_maps_ggplot.png` | Spatial maps with country borders (ggplot2) |
| `*_qmatrix_K*.csv` | Ancestry coefficients for each K |

### Q-Matrix CSV Structure

| Column | Description |
|--------|-------------|
| `Sample` | Sample ID |
| `Cluster_1` to `Cluster_K` | Ancestry proportion for each cluster (sum to 1) |
| `Longitude` | Sample longitude |
| `Latitude` | Sample latitude |

---

## Interpretation Guide

### Understanding Cross-Validation

**CV Score interpretation:**
- Lower CV score = better fit
- Look for "elbow" in CV plot
- Minimum CV typically = optimal K
- Plateau suggests additional K not informative

**Example:**
```
K=1: CV = 0.45
K=2: CV = 0.32  ← Major drop
K=3: CV = 0.28  ← Moderate drop
K=4: CV = 0.27  ← Small drop
K=5: CV = 0.26  ← Plateau
```
Optimal K likely = 3 or 4

### Interpreting Ancestry Bar Plots

**STRUCTURE-style plots show:**
- Each vertical bar = one individual
- Colors = ancestry clusters
- Height = ancestry proportion
- Individuals sorted by similarity

**Patterns:**

**Clear structure:**
- Distinct color blocks
- Minimal mixing
- Corresponds to geography

**Admixture:**
- Mixed colors within individuals
- Gradual transitions
- Gene flow between clusters

**Weak structure:**
- Similar colors across all individuals
- Minimal differentiation
- Possible isolation by distance only

### Interpreting Spatial Maps

**Color patterns indicate:**

1. **Distinct regions** (sharp boundaries)
   - Strong population differentiation
   - Geographic barriers to gene flow
   - Possible adaptation to local conditions

2. **Gradual transitions** (smooth clines)
   - Isolation by distance
   - Continuous gene flow
   - No major barriers

3. **Patchwork patterns**
   - Complex admixture
   - Multiple colonization events
   - Recent disturbance

**Geographic correspondence:**
- Boundaries at mountains, rivers
- Clusters match ecoregions
- Transitions along corridors

### Choosing Optimal K

**Multiple criteria:**

1. **Statistical**: Minimum CV score
2. **Biological**: Matches known subspecies/populations
3. **Geographic**: Corresponds to landscape features
4. **Practical**: Interpretable number of clusters

**Common K values:**

| K | Typical Pattern |
|---|-----------------|
| 2 | Major geographic split (e.g., east vs west) |
| 3-5 | Regional structure (e.g., north, central, south) |
| 6-10 | Fine-scale or family structure |
| >10 | May be overfitting |

### Management Implications

**For conservation:**

**Distinct clusters (high differentiation):**
- Manage as separate units
- Avoid mixing in translocations
- Prioritize unique genetic diversity
- Maintain connectivity between units

**Admixture zones:**
- Preserve as gene flow corridors
- May have high genetic diversity
- Important for adaptive potential
- Buffer against climate change

**For breeding/restoration:**

**Source populations:**
- High ancestry in target cluster
- Locally adapted
- Use for seed source

**Outcrossing:**
- Cross moderately differentiated clusters
- Avoid extreme admixture
- Consider local adaptation

---

## Advanced Options

### Testing Different K Ranges

```r
# Fine-scale analysis
K_VALUES <- 1:15

# Broad-scale only
K_VALUES <- 2:5

# Specific values
K_VALUES <- c(2, 4, 6, 8)
```

### Adjusting Spatial Resolution

```r
# Higher resolution maps (slower)
interpol = FieldsKrigModel(20)

# Lower resolution (faster)
interpol = FieldsKrigModel(5)

# Change map resolution
resolution = c(500, 500)  # Higher detail
resolution = c(100, 100)  # Lower detail
```

### Custom Color Palettes

```r
# Define your own colors
my.colors <- c("red", "blue", "green", "yellow", "purple")
my.palette <- CreatePalette(my.colors, 9)

# Use RColorBrewer palettes
library(RColorBrewer)
my.palette <- brewer.pal(7, "Set1")
```

### Filtering Samples by Geography

```r
# Example: Keep only samples in specific region
lon_range <- c(-130, -110)
lat_range <- c(30, 45)

keep <- coord_matrix[,1] >= lon_range[1] & 
        coord_matrix[,1] <= lon_range[2] &
        coord_matrix[,2] >= lat_range[1] & 
        coord_matrix[,2] <= lat_range[2]

genotype_subset <- genotype_final[keep, ]
coord_subset <- coord_matrix[keep, ]
```

### Comparing with STRUCTURE

```r
# Export Q-matrices for comparison
q.matrix <- qmatrix(tess3_obj, K = 3)

# Convert to STRUCTURE format
# Compare clustering patterns
# Calculate correlation between methods
```

---

## Troubleshooting

### Common Issues

#### 1. "No matching samples found"

**Problem:** Sample IDs don't match between VCF and CSV

**Solutions:**
```r
# Check sample names
vcf_samples <- colnames(gt_numeric)
coord_samples <- coordinates$SampleID

# Find mismatches
setdiff(vcf_samples, coord_samples)
setdiff(coord_samples, vcf_samples)

# Common issues:
# - Whitespace differences
# - Case sensitivity
# - Underscores vs dashes
```

#### 2. "Missing coordinates"

**Problem:** Some samples lack lat/long data

**Solutions:**
- Remove samples with missing coordinates
- Obtain coordinates from collection records
- Use approximate locations if appropriate
- Check coordinate column names match parameters

#### 3. TESS3 fails or crashes

**Possible causes:**

**Insufficient memory:**
```r
# Reduce SNP number
# Filter more stringently
# Use LD-pruned SNPs
```

**Too many NAs:**
```r
# Impute before running
genotype <- apply(genotype, 2, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
})
```

**Wrong matrix orientation:**
```r
# TESS3 requires: samples as rows, SNPs as columns
dim(genotype)  # Should be (samples x SNPs)
genotype <- t(gt_numeric)  # Transpose if needed
```

#### 4. "CV scores all similar"

**Problem:** No clear optimal K

**Possible causes:**
- Weak population structure
- Isolation by distance only
- Need larger K range
- Sample size too small

**Solutions:**
- Test higher K values
- Check with PCA
- Consider IBD models
- Increase sample size

#### 5. Spatial maps look wrong

**Problems:**

**Map extent too large/small:**
```r
# Adjust map limits
MAP_XLIM <- c(min(lon)-5, max(lon)+5)
MAP_YLIM <- c(min(lat)-5, max(lat)+5)
```

**Missing map features:**
```r
# Check rworldmap installation
library(rworldmap)
map.polygon <- getMap(resolution = "low")
```

**Interpolation issues:**
```r
# Adjust kriging model
interpol = FieldsKrigModel(k = 10)  # Change k value
```

---

## Best Practices

### 1. Study Design

**Optimal sampling:**
- **50-200 samples** ideal
- Distributed across range
- 5-10 samples per location
- Cover environmental gradients
- Avoid large geographic gaps

### 2. Data Quality

**Before analysis:**
- Filter VCF (MAF > 0.05, missing < 20%)
- Remove monomorphic SNPs
- Check coordinate accuracy
- Verify coordinate system (WGS84)
- Remove duplicate samples

### 3. K Selection

**Best practices:**
- Test K from 1 to 10 (minimum)
- Use cross-validation
- Consider biological knowledge
- Check consistency across runs
- Compare with other methods (PCA, STRUCTURE)

### 4. Validation

**Recommended:**
- Compare TESS3 with STRUCTURE/ADMIXTURE
- Run PCA for comparison
- Calculate pairwise FST
- Test for isolation by distance
- Validate with independent data

### 5. Visualization

**For publication:**
- Use high resolution (300+ dpi)
- Clear color schemes
- Appropriate map extent
- Add scale bars and north arrows
- Include sample points on maps

---

## Computational Performance

### Memory Requirements

| Dataset Size | RAM Required |
|-------------|--------------|
| < 50 samples, < 10K SNPs | 4GB |
| 50-100 samples, 10-50K SNPs | 8GB |
| 100-200 samples, 50-100K SNPs | 16GB |
| > 200 samples, > 100K SNPs | 32GB+ |

### Speed Optimization

**To speed up:**

1. Use multiple CPU cores
2. Reduce SNP number (LD-prune)
3. Test fewer K values
4. Lower spatial resolution
5. Filter samples geographically

---

## Citation

If you use this workflow, please cite:

**For TESS3:**
- Caye, K., Deist, T. M., Martins, H., Michel, O., & François, O. (2016). TESS3: Fast inference of spatial population structure and genome scans for selection. *Molecular Ecology Resources*, 16(2), 540-548.

**For tess3r package:**
- Frichot, E., Mathieu, F., Trouillon, T., Bouchard, G., & François, O. (2014). Fast and efficient estimation of individual ancestry coefficients. *Genetics*, 196(4), 973-983.

**Example citation:**
> "We analyzed spatial population structure using TESS3 (Caye et al., 2016), which incorporates geographic coordinates to model continuous variation in ancestry. We tested K from 1 to 7 and selected optimal K based on cross-validation scores."

---

## Additional Resources

### Related Methods

- **STRUCTURE**: Bayesian clustering without coordinates
- **ADMIXTURE**: Fast maximum likelihood clustering
- **sNMF**: Sparse non-negative matrix factorization
- **GENELAND**: Spatial clustering with genetic discontinuities
- **EEMS**: Effective migration surfaces

### Tutorials

- TESS3 vignette: In R package documentation
- tess3r GitHub: https://github.com/bcm-uga/TESS3_encho_sen
- Landscape genomics: https://popgen.nescent.org/

### Recommended Reading

1. **Spatial analysis:**
   - Guillot, G., et al. (2005). GENELAND: a computer package for landscape genetics. *Molecular Ecology Notes*, 5(3), 712-715.

2. **Methods comparison:**
   - Jombart, T., et al. (2010). Discriminant analysis of principal components. *BMC Genetics*, 11(1), 94.

3. **Landscape genetics:**
   - Storfer, A., et al. (2010). Putting the 'landscape' in landscape genetics. *Heredity*, 105(2), 128-133.

---

## License

This workflow is provided under the MIT License.

## Contact

For questions about this workflow:
- Review troubleshooting section
- Check TESS3 documentation: `?tess3`
- Post on GitHub issues or Biostars

---

**Last updated:** December 2024

**Version:** 1.0

**Tested with:** R 4.3.0, tess3r 1.1.0
