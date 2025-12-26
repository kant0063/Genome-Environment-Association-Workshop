# Genomic Selection Pipeline v1.9

**Complete genomic prediction pipeline with environmental adaptation analysis and global geographic visualization**

---

## 🚀 Features

### Core Functionality
- ✅ **VCF Processing**: Automated conversion and quality control
- ✅ **Environmental Data Extraction**: WorldClim bioclimatic variables + custom GeoTIFFs
- ✅ **Smart Sample Matching**: Automatic alignment of genotypes and environmental data
- ✅ **SNP Filtering**: Flexible missing data tolerance
- ✅ **Genotype Imputation**: Mean, median, or mode imputation
- ✅ **Core Selection**: Greedy algorithm for diverse training set
- ✅ **Parallel Processing**: Multi-core support for speed

### Genomic Prediction Methods
1. **rrBLUP** - Ridge regression BLUP (fast, reliable)
2. **GBLUP** - Genomic BLUP with additive relationship matrix
3. **GBLUP-Gaussian** - Gaussian kernel for non-linear relationships
4. **GBLUP-Exponential** - Exponential kernel for distance decay
5. **BayesCπ** - Bayesian variable selection (requires BGLR package)

### Visualization & Analysis
- ✅ **Cross-Validation Plot**: Comparative prediction accuracy across methods
- ✅ **Regional Heatmap**: Geographic patterns (15 broad regions)
- ✅ **Country-Specific Heatmap**: Fine-scale patterns (160+ countries/regions)
- ✅ **Global Coverage**: All continents including comprehensive African coverage
- ✅ **Publication-Ready**: 300 DPI, 18×12" figures

---

## 📦 Installation

### Required R Packages

```r
# Core packages (required)
install.packages(c("rrBLUP", "SNPRelate", "dplyr", "ggplot2", 
                   "raster", "parallel", "foreach", "doParallel", 
                   "reshape2", "ggnewscale"))

# Optional (for BayesCπ method)
install.packages("BGLR")
```

### Download Pipeline

```bash
# Download and extract
tar -xzf genomic_pipeline_v1.9_GLOBAL_COVERAGE.tar.gz
# or
unzip genomic_pipeline_v1.9_GLOBAL_COVERAGE.zip
```

---

## 🎯 Quick Start

### Minimal Example

```r
source("genomic_pipeline_FINAL.R")

results <- run_genomic_pipeline_optimized(
  vcf_file = "genotypes.vcf.gz",
  coords_file = "coordinates.csv",
  worldclim_path = "wc2.0_30s_bio/",
  n_core = 50,
  output_dir = "results"
)
```

### Full Example with All Options

```r
results <- run_genomic_pipeline_optimized(
  # Input files
  vcf_file = "genotypes.vcf.gz",
  coords_file = "coordinates.csv",
  worldclim_path = "wc2.0_30s_bio/",
  
  # Core selection
  n_core = 50,
  max_snps_core = 10000,
  
  # SNP filtering
  max_missing = 0.05,  # Allow 5% missing per SNP
  
  # Imputation
  impute_method = "mean",  # "mean", "median", or "mode"
  
  # Prediction methods
  methods = c("rrBLUP", "GBLUP-Gaussian", "BayesCpi"),
  
  # Cross-validation
  k_fold = 10,
  n_reps = 25,
  
  # Performance
  n_cores = 10,
  
  # Environmental data
  use_all_tifs = TRUE,  # Use all GeoTIFF files, not just bio_01-19
  
  # Output
  output_dir = "lentil_results"
)
```

---

## 📁 Input File Formats

### 1. VCF File
- **Format**: Standard VCF or VCF.gz
- **Requirements**: SNP genotypes for all samples
- **Example**: `genotypes.vcf.gz`

### 2. Coordinates File (CSV)
- **Required columns**: `ID`, `Latitude`, `Longitude`
- **Example**:
```csv
ID,Latitude,Longitude
Sample001,41.9028,12.4964
Sample002,36.7783,119.4179
Sample003,-33.9249,18.4241
```

### 3. WorldClim Data
- **Format**: GeoTIFF (.tif) files
- **Variables**: bio_01 through bio_19 (bioclimatic variables)
- **Plus**: Any additional environmental rasters
- **Download**: https://worldclim.org/data/worldclim21.html
- **Resolution**: 30 seconds (~1km) recommended

---

## ⚙️ Parameters Explained

### SNP Missing Data Tolerance (`max_missing`)

Controls how many missing genotypes are allowed per SNP:

```r
max_missing = 0      # Strict: only complete SNPs (default)
max_missing = 0.05   # Allow SNPs with ≤5% missing samples
max_missing = 0.1    # Allow SNPs with ≤10% missing samples  
max_missing = 5      # Allow SNPs missing in ≤5 samples (absolute count)
```

**Recommendation**: 
- High coverage data: `0` (strict)
- Medium coverage: `0.05` (5%)
- Low coverage: `0.1-0.2` (10-20%)

### Imputation Method (`impute_method`)

How to fill in missing genotypes:

```r
impute_method = "mean"    # Mean genotype value (default, fastest)
impute_method = "median"  # Median value (robust to outliers)
impute_method = "mode"    # Most common value (for discrete data)
```

### Prediction Methods (`methods`)

Which genomic prediction models to use:

```r
# Fast and reliable
methods = c("rrBLUP")

# Multiple methods for comparison
methods = c("rrBLUP", "GBLUP-Gaussian", "BayesCpi")

# All methods (long runtime)
methods = c("rrBLUP", "GBLUP", "GBLUP-Gaussian", 
            "GBLUP-Exponential", "BayesCpi")
```

**Runtime estimates** (50 training samples, 103 total, 36 traits, 10 cores):
- rrBLUP: ~5 min CV + 10 sec prediction
- GBLUP-Gaussian: ~15 min CV + 30 sec prediction
- GBLUP-Exponential: ~15 min CV + 30 sec prediction
- BayesCπ: ~45 min CV + 2 min prediction

---

## 📊 Output Files

### For Each Method

```
results/
├── cv_plot.png                              # Cross-validation accuracy
│
├── gebv_heatmap_regional_rrBLUP.png        # Regional grouping (15 groups)
├── gebv_heatmap_country_rrBLUP.png         # Country-specific (160+ groups)
│
├── GEBVs_rrBLUP_all.csv                    # Predictions for all samples
├── GEBVs_rrBLUP_training.csv               # Predictions for training only
├── cv_rrBLUP.rds                           # Cross-validation results
│
└── [same files for each method]
```

### CSV File Format

**GEBVs_[method]_all.csv**:
```csv
Sample,bio_01,bio_02,bio_03,...,bio_19
Sample001,12.34,-1.23,0.45,...,234.56
Sample002,15.67,2.34,-0.12,...,198.23
...
```

---

## 🌍 Geographic Coverage

### Regional Heatmap (~15 groups)
Broad continental regions for overview:
- Turkey, Italy, Spain, Greece, France, Europe-Other
- Iran, Syria-Iraq, Egypt, North-Africa-Middle-East
- India, Pakistan, Bangladesh, South-Asia
- China, East-Asia
- USA-Canada, Americas
- Sub-Saharan-Africa

### Country-Specific Heatmap (160+ groups)
Fine-scale country/regional resolution:

**Africa (40+ regions)**:
Morocco, Algeria, Tunisia, Libya, Egypt, Sudan, Ethiopia, Nigeria, Kenya, Tanzania, South-Africa, Madagascar, and 30+ more

**Europe (40+ regions)**:
Individual countries including Turkey, Italy, Spain, France, Germany, UK, Poland, Ukraine, Russia, and Balkan states

**Asia (40+ regions)**:
Afghanistan, Pakistan, India, Bangladesh, China, Korea, Japan, Thailand, Vietnam, Indonesia, Philippines, and more

**Americas (20+ regions)**:
USA, Canada, Mexico, Brazil, Argentina, Chile, Colombia, and Central American countries

**Oceania (5+ regions)**:
Australia, New-Zealand, Papua-New-Guinea, Pacific-Islands

---

## 🔬 Workflow

```
1. Load VCF → All samples
2. Extract environmental data → For samples with coordinates
3. Filter samples → Keep those with complete environmental data
4. SNP filtering → Based on max_missing tolerance
5. Core selection → Diverse training set
6. Genotype imputation → Fill remaining missing values
7. Cross-validation → Test prediction accuracy (training set)
8. Genomic selection → Predict all samples
9. Generate plots → CV plot + 2 heatmaps per method
```

---

## 💡 Tips & Best Practices

### For Best Results

1. **Sample Size**: 
   - Minimum 50 training samples
   - More diverse = better predictions

2. **SNP Filtering**:
   - Start with `max_missing = 0.05`
   - Adjust based on SNP retention

3. **Method Selection**:
   - rrBLUP: Fast, good baseline
   - GBLUP-Gaussian: Better for non-linear relationships
   - BayesCπ: Best for complex traits (but slowest)

4. **Computational Resources**:
   - Use `n_cores = [available cores - 1]`
   - More cores = faster CV
   - BayesCπ is RAM-intensive

### Troubleshooting

**Error: "No SNPs after filtering"**
- Increase `max_missing` parameter
- Check VCF file quality

**Error: "Not enough core samples"**
- Decrease `n_core` parameter
- Check if samples have complete environmental data

**Error: "Kernel computation failed"**
- Likely numerical issues with small sample size
- Try rrBLUP or GBLUP instead of kernel methods

---

## 🏃 HPC Usage (Mana Cluster)

Edit paths in `run_genomic_pipeline_mana.sh` then submit:

```bash
sbatch run_genomic_pipeline_mana.sh
```

---

## 📖 Citation

If you use this pipeline in your research, please cite:

```
Genomic Selection Pipeline v1.9
University of Hawaii at Manoa
https://github.com/[your-repo]
```

---

## 🤝 Contributing

Issues, suggestions, and contributions welcome!

---

## 📧 Contact

**Pipeline Support**: mbkantar@hawaii.edu  
**HPC Support**: mana-help@hawaii.edu

---

## 🔄 Version History

### v1.9 (Current - 2024)
- ✅ Complete global coverage (160+ countries/regions)
- ✅ Comprehensive African country coverage (40+ regions)
- ✅ All traits included in heatmaps
- ✅ Dual heatmap system (regional + country-specific)

### v1.8
- ✅ Added country-specific heatmaps
- ✅ Genotype imputation (mean/median/mode)
- ✅ All traits in visualizations

### v1.7
- ✅ GEBV heatmaps with geographic clustering
- ✅ Improved kernel methods with regularization

### v1.6
- ✅ Genotype imputation step
- ✅ Multiple imputation methods

### v1.5
- ✅ Multiple prediction methods
- ✅ SNP missing data tolerance
- ✅ Comparative visualizations

### v1.3-1.4
- ✅ Fixed VCF conversion (transpose-first method)
- ✅ Automatic sample matching
- ✅ 13x speedup with parallel processing

### Earlier Versions
- Initial pipeline development
- Core functionality implementation

---

**Last Updated**: December 2024
