# Genetic Offset Analysis for Climate Change Vulnerability

A comprehensive R workflow for predicting population vulnerability to climate change using genomic data and future climate projections.

## Overview

**Genetic offset** quantifies the mismatch between a population's current genetic composition and the optimal genetic composition under future climate conditions. This metric serves as a powerful indicator of climate change vulnerability and helps identify populations that may require conservation intervention.

### What This Analysis Does

1. Identifies genotype-environment associations using LFMM (Latent Factor Mixed Models)
2. Downloads future climate projections from CMIP6 models
3. Calculates genetic offset across geographic space
4. Produces global and regional vulnerability maps
5. Generates high-resolution GeoTIFF outputs for GIS integration

### Key Applications

- **Conservation planning**: Identify vulnerable populations requiring intervention
- **Assisted gene flow**: Guide translocation and restoration decisions
- **Climate adaptation**: Predict which populations need support
- **Breeding programs**: Select climate-resilient genotypes
- **Protected area design**: Identify climate refugia

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
install.packages(c("vcfR", "terra", "dplyr", "fields", "maps"))

# LEA and LFMM (may require Bioconductor)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("LEA", "lfmm"))

# geodata for climate data download
install.packages("geodata")

# Optional: For R Markdown
install.packages(c("knitr", "kableExtra", "rmarkdown"))
```

### System Requirements

- **R version**: ≥ 4.0.0
- **RAM**: At least 16GB recommended (32GB for large datasets)
- **Storage**: Minimum 10GB free space for climate data downloads
- **Internet**: Required for downloading WorldClim and CMIP6 data

### CMIP6 Data Access

Climate data is downloaded automatically from WorldClim servers. No special accounts or authentication required.

---

## Quick Start

### Minimal Example

```r
# Set your file paths
VCF_FILE <- "my_genotypes.vcf"
ENV_CSV <- "sample_coordinates.csv"

# Run analysis with default parameters
source("Genetic_Offset_Analysis.R")
```

### Using R Markdown

```r
# Open in RStudio and click "Knit"
# Or from command line:
rmarkdown::render("Genetic_Offset_Analysis.Rmd")
```

### Expected Runtime

- **Small dataset** (50 samples, 10K SNPs): ~30 minutes
- **Medium dataset** (100 samples, 50K SNPs): ~2-4 hours
- **Large dataset** (200+ samples, 100K+ SNPs): ~6-12 hours

*Most time is spent downloading climate data and running LFMM*

---

## Input Data Requirements

### 1. VCF File (Genotype Data)

**Format:** Standard VCF format (`.vcf` or `.vcf.gz`)

**Requirements:**
- SNP genotype calls for all samples
- Quality-filtered SNPs
- Biallelic SNPs only

**Recommended filtering:**
```bash
# Using VCFtools
vcftools --vcf input.vcf \
  --maf 0.05 \
  --max-missing 0.8 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --out filtered
```

**Example VCF:**
```
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1  Sample2
Chr1    1000    .       A    G    99    PASS    .     GT      0/0      0/1
Chr1    2000    .       C    T    99    PASS    .     GT      0/1      1/1
```

### 2. Sample Coordinates CSV

**Required columns:**
- **First column**: Sample IDs (must match VCF exactly)
- **Latitude column**: Named "latitude", "lat", or "Latitude"
- **Longitude column**: Named "longitude", "lon", "long", or "Longitude"

**Optional columns:**
- Any additional metadata (population, location, etc.)

**Example CSV:**
```csv
SampleID,Latitude,Longitude,Population,Site
Sample1,45.5231,-122.6765,PopA,Site1
Sample2,46.8721,-121.7689,PopA,Site2
Sample3,42.3601,-71.0589,PopB,Site3
Sample4,37.7749,-122.4194,PopB,Site4
```

**Important notes:**
- Sample IDs must match VCF exactly (case-sensitive)
- Coordinates in decimal degrees (WGS84)
- Missing coordinates not allowed
- Latitude range: -90 to 90
- Longitude range: -180 to 180

### Data Preparation Checklist

Before running the analysis:

- [ ] VCF filtered for quality (MAF > 0.05, missing < 20%)
- [ ] Sample IDs match between VCF and CSV
- [ ] Coordinates are in decimal degrees (not degrees/minutes/seconds)
- [ ] At least 30 samples recommended
- [ ] Geographic coverage spans environmental gradients
- [ ] No duplicate sample IDs

---

## Usage

### Basic Workflow

1. **Prepare your data** (see Input Data Requirements)
2. **Edit parameters** at the top of the script
3. **Run the analysis**
4. **Examine outputs** and interpret results

### Customizing Parameters

Edit these parameters in the script:

```r
# ==================== INPUT FILES ====================
VCF_FILE <- "your_genotypes.vcf"
ENV_CSV <- "your_coordinates.csv"

# ==================== WORLDCLIM RESOLUTION ====================
RES_MINUTES <- 10    # Options: 2.5, 5, or 10 minutes
                     # Higher resolution = more detail but slower

# ==================== CLIMATE MODELS ====================
CMIP6_MODELS <- c(
  "ACCESS-ESM1-5",   # Australia
  "BCC-CSM2-MR",     # China
  "CanESM5",         # Canada
  "CNRM-CM6-1",      # France
  "IPSL-CM6A-LR",    # France
  "MIROC6"           # Japan
)

# ==================== SCENARIOS AND PERIODS ====================
SSP_SCENARIOS <- c("245", "585")
# "245" = SSP2-4.5 (moderate emissions)
# "585" = SSP5-8.5 (high emissions)

FUTURE_PERIODS <- c("2041-2060", "2081-2100")
# Mid-century and end-of-century

# ==================== ANALYSIS PARAMETERS ====================
GRID_RESOLUTION <- 200       # Grid cells for projection (50-500)
K_LFMM <- 7                  # Latent factors (typically 5-10)
PVALUE_CUTOFF <- 1e-5        # Significance threshold (1e-4 to 1e-6)
MAX_PROJ_CLIP <- NA          # Clip extreme values (optional)

# ==================== REGIONAL PLOTTING ====================
REGIONAL_BUFFER_DEGREES <- 2  # Buffer around data for regional maps
```

### Parameter Selection Guide

**Grid Resolution:**
- 50-100: Coarse, fast, global overview
- 200: Balanced (default), good for most analyses
- 300-500: Fine detail, slow, high memory use

**K (Latent Factors):**
- Use genetic structure analysis (PCA, ADMIXTURE) to inform
- Typically K = number of genetic clusters
- Too low: Poor correction for structure
- Too high: Loss of power
- Default K=7 works for most datasets

**P-value Cutoff:**
- 1e-4: Liberal, more candidates
- 1e-5: Balanced (default)
- 1e-6: Conservative, fewer candidates

### Running from Command Line

```bash
# Run R Markdown and generate HTML report
Rscript -e "rmarkdown::render('Genetic_Offset_Analysis.Rmd')"

# Generate PDF (requires LaTeX)
Rscript -e "rmarkdown::render('Genetic_Offset_Analysis.Rmd', output_format='pdf_document')"

# Run R script directly
Rscript Genetic_Offset_Analysis.R
```

---

## Output Files

### Directory Structure

```
genetic_offset_results/
├── candidate_snps_index.csv
├── LFMM_B_matrix.csv
├── per_model_debug_info.rds
├── data_extent.csv
├── offset_GLOBAL_MEAN.tif
├── offset_ssp245_2041-2060_MEAN.tif
├── offset_ssp245_2081-2100_MEAN.tif
├── offset_ssp585_2041-2060_MEAN.tif
├── offset_ssp585_2081-2100_MEAN.tif
└── regional_plots/
    ├── REGIONAL_offset_GRAND_MEAN.tif
    ├── REGIONAL_offset_GRAND_MEAN.png
    └── [additional regional files]

cmip6_cache/
└── [downloaded climate data - cached for reuse]
```

### File Descriptions

| File | Description |
|------|-------------|
| `candidate_snps_index.csv` | Indices of SNPs with significant genotype-environment associations |
| `LFMM_B_matrix.csv` | Effect size matrix from LFMM (SNPs × environmental variables) |
| `per_model_debug_info.rds` | Metadata for each climate model run |
| `data_extent.csv` | Geographic extent of sample data |
| `offset_GLOBAL_MEAN.tif` | **Main result**: Grand mean genetic offset across all scenarios |
| `offset_ssp*_MEAN.tif` | Mean genetic offset for specific SSP scenario and time period |
| `REGIONAL_*.tif` | Regional subsets focused on data extent |

### GeoTIFF Files

All raster outputs are in **GeoTIFF format** with:
- Coordinate Reference System: WGS84 (EPSG:4326)
- Projection: Geographic (latitude/longitude)
- Units: Decimal degrees
- Values: Genetic offset (unitless, relative measure)

**Loading in R:**
```r
library(terra)
offset <- rast("genetic_offset_results/offset_GLOBAL_MEAN.tif")
plot(offset)
```

**Loading in QGIS:**
1. Open QGIS
2. Layer → Add Raster Layer
3. Select `.tif` file
4. Apply color ramp (e.g., Spectral, RdYlBu)

**Loading in ArcGIS:**
1. Add Data → navigate to `.tif` file
2. Right-click → Properties → Symbology
3. Choose color scheme

---

## Interpretation Guide

### Understanding Genetic Offset

**What is genetic offset?**

Genetic offset measures the magnitude of environmental change relative to current genotype-environment relationships. It represents how "mismatched" a population's genetics will be to future environmental conditions.

**Units:** Genetic offset is a relative, unitless metric. Values are meaningful in comparison to each other, not as absolute measures.

### Interpreting Values

**General guidelines:**

| Offset Value | Interpretation | Management Implication |
|--------------|----------------|------------------------|
| < 0.5 | Low vulnerability | Monitor, natural adaptation likely |
| 0.5 - 1.0 | Moderate vulnerability | Enhanced monitoring, consider habitat protection |
| 1.0 - 1.5 | High vulnerability | Priority for intervention, consider assisted migration |
| > 1.5 | Very high vulnerability | Urgent action needed, high extinction risk |

**Important caveats:**

- Values are **relative** to your dataset
- Compare **within** your analysis, not across studies
- Higher values = greater environmental change
- Does not directly predict extinction

### Spatial Patterns

**What to look for:**

1. **Hotspots (High Offset Areas)**
   - Regions with greatest climate change impact
   - Populations most vulnerable
   - Priority areas for conservation intervention
   - Candidates for assisted gene flow

2. **Refugia (Low Offset Areas)**
   - Regions where current genetics match future climate
   - Potential climate refugia
   - Source populations for translocation
   - High conservation value

3. **Gradients**
   - Directional patterns in vulnerability
   - May guide translocation strategies
   - Indicate direction of climate velocity

### Comparing Scenarios

**SSP2-4.5 vs SSP5-8.5:**
- SSP2-4.5: Moderate emissions scenario
- SSP5-8.5: High emissions scenario
- Higher emissions = generally higher offset

**2041-2060 vs 2081-2100:**
- Mid-century vs end-of-century
- Later period = generally higher offset
- Shows temporal trajectory

### Management Recommendations

**For HIGH offset populations:**

✓ **Intensive monitoring**: Track population responses  
✓ **Habitat protection**: Secure quality habitat  
✓ **Assisted migration**: Consider translocation  
✓ **Ex-situ conservation**: Establish seed banks  
✓ **Breeding programs**: Select for climate resilience

**For LOW offset populations:**

✓ **Protect as refugia**: High conservation priority  
✓ **Source for restoration**: Use in translocation  
✓ **Study adaptation**: Investigate mechanisms  
✓ **Connectivity**: Ensure gene flow pathways

### Validation Approaches

**How to validate predictions:**

1. **Monitor populations**: Track fitness in response to climate change
2. **Common gardens**: Test populations under different climates
3. **Transplant experiments**: Move populations to predicted future climates
4. **Independent datasets**: Compare with other populations/species
5. **Functional validation**: Test candidate genes in lab

---

## Advanced Options

### Using Different Climate Models

Add or modify CMIP6 models:

```r
CMIP6_MODELS <- c(
  "ACCESS-CM2",
  "GFDL-ESM4",
  "GISS-E2-1-G",
  "HadGEM3-GC31-LL",
  "MRI-ESM2-0",
  "UKESM1-0-LL"
)
```

See available models: https://worldclim.org/data/cmip6/cmip6climate.html

### Adjusting Spatial Resolution

**For faster computation:**
```r
GRID_RESOLUTION <- 100      # Coarser grid
RES_MINUTES <- 10           # Coarser climate data
```

**For higher detail:**
```r
GRID_RESOLUTION <- 400      # Finer grid
RES_MINUTES <- 2.5          # Finer climate data (WARNING: large files!)
```

### Clipping Extreme Values

To prevent outliers from dominating color scales:

```r
MAX_PROJ_CLIP <- 3.0  # Clip offset values above 3.0
```

### Custom Environmental Variables

To use different environmental predictors:

```r
# After extracting climate data, subset to specific variables
X.env <- X.env[, c("bio1", "bio12", "bio15")]  # Temperature, precip, seasonality
```

### Regional Analysis Only

To skip global projections and focus on regional:

```r
# Set a smaller grid focused on your region
long.vec <- seq(min_lon - 5, max_lon + 5, length = 100)
lat.vec <- seq(min_lat - 5, max_lat + 5, length = 100)
```

### Multiple K Values

Test different numbers of latent factors:

```r
K_values <- c(3, 5, 7, 10)

for (K in K_values) {
  # Run analysis with current K
  # Save results with K in filename
}
```

---

## Troubleshooting

### Common Issues

#### 1. "No candidate SNPs found"

**Problem:** LFMM doesn't find significant genotype-environment associations

**Solutions:**
```r
# Loosen p-value threshold
PVALUE_CUTOFF <- 1e-4  # Instead of 1e-5

# Reduce K (may be over-correcting)
K_LFMM <- 5

# Check data quality:
# - Enough environmental variation?
# - Samples span environmental gradients?
# - SNPs polymorphic and not overly filtered?
```

#### 2. "Error downloading climate data"

**Problem:** Network issues or server unavailable

**Solutions:**
- Check internet connection
- Try again later (servers may be temporarily down)
- Download manually from https://worldclim.org
- Use cached data from previous runs

#### 3. Memory errors

**Problem:** R runs out of memory

**Solutions:**
```r
# Reduce grid resolution
GRID_RESOLUTION <- 100

# Use fewer climate models
CMIP6_MODELS <- CMIP6_MODELS[1:3]

# Increase R memory limit (Windows)
memory.limit(size = 32000)  # 32GB

# Process scenarios sequentially instead of storing all
```

#### 4. Sample name mismatch

**Problem:** Samples in VCF don't match CSV

**Solutions:**
```r
# Check VCF sample names
vcf <- read.vcfR(VCF_FILE)
colnames(vcf@gt)[-1]

# Check CSV sample names
coords <- read.csv(ENV_CSV)
coords[,1]

# Common issues:
# - Extra spaces
# - Case sensitivity
# - Special characters
# - Duplicates
```

#### 5. LFMM convergence issues

**Problem:** LFMM fails to converge or produces errors

**Solutions:**
```r
# Reduce K
K_LFMM <- 5

# Remove SNPs with high missing data
# Filter VCF before analysis:
vcftools --vcf input.vcf --max-missing 0.9 --recode --out clean

# Check for monomorphic SNPs (should be removed automatically)

# Ensure environmental variables have variation
```

#### 6. "All offset values are NA"

**Problem:** No valid predictions generated

**Causes and solutions:**
- **No common climate layers**: Check that WorldClim versions match
- **Coordinate mismatch**: Verify coordinates are in decimal degrees
- **Environmental data extraction failed**: Check coordinate validity

```r
# Debug by checking intermediate steps:
print(dim(X.env))           # Should have rows and columns
print(summary(X.env))       # Should have values
print(head(coords))         # Check coordinates look correct
```

---

## Best Practices

### 1. Data Quality

**Before running analysis:**
- Filter VCF for quality (MAF > 0.05, missing < 20%, HWE)
- Remove closely related individuals if doing population-level analysis
- Ensure samples span environmental gradients
- Check coordinates are accurate (use GPS data when possible)

### 2. Sample Design

**Optimal sampling:**
- **Minimum 30 samples** (more is better)
- Span full environmental range of species
- Include multiple populations
- Balance across environmental gradients
- Consider temporal replication if possible

### 3. Parameter Selection

**General recommendations:**
- Start with default parameters
- Use population structure analysis to inform K
- Be conservative with p-value thresholds initially
- Compare multiple scenarios and time periods

### 4. Validation

**Always validate predictions:**
- Compare with empirical fitness data if available
- Cross-validate with independent datasets
- Test predictions with transplant experiments
- Compare with other modeling approaches

### 5. Reporting Results

**Include in methods:**
- VCF filtering criteria
- Number of samples and SNPs
- K value and how it was selected
- P-value threshold
- Climate models and scenarios used
- Grid resolution

**Report in results:**
- Number of candidate SNPs
- Environmental variables used
- Mean and range of offset values
- Spatial patterns and hotspots
- Uncertainty across models

---

## Computational Performance

### Memory Requirements

| Dataset Size | Recommended RAM |
|-------------|-----------------|
| Small (< 50 samples, < 10K SNPs) | 8GB |
| Medium (50-100 samples, 10-50K SNPs) | 16GB |
| Large (100-200 samples, 50-100K SNPs) | 32GB |
| Very Large (> 200 samples, > 100K SNPs) | 64GB+ |

### Speed Optimization

**To speed up analysis:**

1. Use fewer climate models (3-4 sufficient)
2. Reduce grid resolution (100-150 adequate)
3. Use higher WorldClim resolution (10 minutes instead of 2.5)
4. Filter SNPs more stringently (MAF > 0.1)
5. Use parallel processing if available

**Parallel processing example:**
```r
library(parallel)
cl <- makeCluster(detectCores() - 1)
# Use parallel apply functions for model loop
stopCluster(cl)
```

### Disk Space

**Climate data cache:**
- 2.5 minute resolution: ~5-10GB per model
- 5 minute resolution: ~2-3GB per model  
- 10 minute resolution: ~500MB per model

Cache is reused across analyses, so only downloaded once.

---

## Citation

If you use this workflow, please cite the following:

**For genetic offset concept:**
- Fitzpatrick, M. C., & Keller, S. R. (2015). Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. *Ecology Letters*, 18(1), 1-16.

**For LFMM method:**
- Caye, K., Deist, T. M., Martins, H., Michel, O., & François, O. (2016). LFMM 2: Fast and accurate inference of gene-environment associations in genome-wide studies. *Molecular Biology and Evolution*, 36(4), 852-860.

**For LEA package:**
- Frichot, E., & François, O. (2015). LEA: An R package for landscape and ecological association studies. *Methods in Ecology and Evolution*, 6(8), 925-929.

**For WorldClim data:**
- Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1‐km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302-4315.

**For CMIP6 climate models:**
- Eyring, V., et al. (2016). Overview of the Coupled Model Intercomparison Project Phase 6 (CMIP6) experimental design and organization. *Geoscientific Model Development*, 9(5), 1937-1958.

**Example citation:**
> "We calculated genetic offset following the approach of Fitzpatrick & Keller (2015), using LFMM2 (Caye et al., 2016) to identify genotype-environment associations and projecting these associations onto future climate conditions from CMIP6 models (Eyring et al., 2016). Climate data were obtained from WorldClim 2.1 (Fick & Hijmans, 2017)."

---

## Additional Resources

### Related Methods

- **GDM** (Generalized Dissimilarity Modeling): Alternative approach to genetic offset
- **GF** (Gradient Forests): Machine learning approach to GEA
- **RDA** (Redundancy Analysis): Constrained ordination for GEA
- **BayPass**: Bayesian approach to environmental association

### Tutorials and Documentation

- LFMM tutorial: http://membres-timc.imag.fr/Olivier.Francois/lfmm/tutorial.html
- LEA vignette: https://bioconductor.org/packages/release/bioc/vignettes/LEA/inst/doc/LEA.pdf
- WorldClim: https://worldclim.org/data/index.html
- CMIP6: https://www.worldclim.org/data/cmip6/cmip6climate.html

### Recommended Reading

1. **Genomic vulnerability:**
   - Bay, R. A., et al. (2018). Genomic signals of selection predict climate-driven population declines in a migratory bird. *Science*, 359(6371), 83-86.

2. **Landscape genomics:**
   - Rellstab, C., et al. (2015). A practical guide to environmental association analysis in landscape genomics. *Molecular Ecology*, 24(17), 4348-4370.

3. **Conservation applications:**
   - Razgour, O., et al. (2019). Considering adaptive genetic variation in climate change vulnerability assessment reduces species range loss projections. *PNAS*, 116(21), 10418-10423.

4. **Methods comparison:**
   - Capblancq, T., & Forester, B. R. (2021). Redundancy analysis: A Swiss Army Knife for landscape genomics. *Methods in Ecology and Evolution*, 12(12), 2298-2309.

---

## License

This workflow is provided under the MIT License. Feel free to modify and distribute as needed.

## Contributing

Suggestions, improvements, and bug reports are welcome! Please submit issues or pull requests on the project repository.

## Contact

For questions about this workflow:
- Review the troubleshooting section
- Consult package documentation
- Post questions on Biostars or r-sig-geo mailing list

For issues with specific packages:
- LEA/LFMM: https://github.com/bcm-uga/lfmm
- geodata: https://github.com/rspatial/geodata
- terra: https://github.com/rspatial/terra

---

**Last updated:** December 2024

**Version:** 1.0

**Tested on:** R 4.3.0, Ubuntu 22.04, macOS 13.0, Windows 11
