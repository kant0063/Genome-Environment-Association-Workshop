# Quick Start Guide: Environmental Association Analysis
## Adapting the Workflow to Your Species

---

## Prerequisites

✅ **Data Requirements:**
- [ ] Population metadata with coordinates
- [ ] Filtered SNP genotype data (hapmap and matrix format)
- [ ] Environmental/climate data matched to samples
- [ ] eGWAS results (from previous analysis)
- [ ] Population genetics statistics (F~ST~, SPA scores)
- [ ] GFF annotation file (optional but recommended)

✅ **Software Requirements:**
- [ ] R (≥ 4.0)
- [ ] RStudio (recommended)
- [ ] Required R packages installed

---

## Quick Setup (5 minutes)

### 1. Install Required Packages

```r
# Core packages
install.packages(c("tidyverse", "sommer", "raster", "readxl", 
                   "qvalue", "slider", "patchwork", "cowplot"))

# Optional packages
install.packages(c("neyhart", "snps", "LDheatmap", "ape"))
```

### 2. Organize Your Files

Create this folder structure:

```
your_project/
├── data/
│   ├── population_metadata_and_genotypes.RData
│   ├── germplasm_origin_bioclim_data.RData
│   ├── species_annotations.gff (optional)
├── results/
│   ├── eaa_gwas_results.RData
│   ├── population_genetics_stats.RData
├── figures/
├── environmental_association_analysis.Rmd
└── DOCUMENTATION_environmental_association_analysis.md
```

### 3. Configure for Your Species

Open `environmental_association_analysis.Rmd` and update these parameters:

```r
# === EDIT THESE ===
species_name <- "YourSpecies"           # e.g., "Pepper", "Lentil"
ld_decay_distance <- 17000              # YOUR LD decay distance (bp)
wild_category <- "Wild"                 # Match your metadata category
fdr_threshold <- 0.20                   # Adjust as needed
gff_file_name <- "your_species.gff"     # or NULL if unavailable
```

---

## Essential Parameters to Determine

### 1. LD Decay Distance 🔍

**What it is:** Distance (bp) at which linkage disequilibrium (r²) drops to ~0.2

**How to find it:**
- Run LD decay analysis on your data
- Check published papers on your species or close relatives
- Conservative estimates by crop type:
  - Selfing annuals: 50-200 kb
  - Outcrossing annuals: 10-50 kb  
  - Perennials: 1-10 kb
  - Clonal species: >500 kb

**Default if unknown:** 20,000 bp (20 kb)

### 2. Wild Population Category 🌿

**What it is:** The label in your `pop_metadata$category` for wild/native/landrace samples

**How to find it:**
```r
load("data/population_metadata_and_genotypes.RData")
table(pop_metadata$category)
```

**Common values:** "Wild", "Landrace", "Native", "Unimproved", "Primitive"

### 3. Model Selection Strategy 📊

**Options:**

1. **"best_lambda"** (recommended) - Uses model with best genomic inflation factor
   - Requires `p_lambda` object with `best_model` column
   
2. **"manual"** - Specify which model to use
   - Set `manual_model <- "model2"` (or whichever you prefer)
   
3. **"all_models"** - Use all models
   - Increases multiple testing burden

**If unsure:** Use `"best_lambda"` if available, otherwise `"manual"` with your most conservative model

---

## Pre-Flight Checklist

Before running the full analysis:

### Step 1: Load and Check Data

```r
# Load main data
load("data/population_metadata_and_genotypes.RData")

# Verify required objects exist
required <- c("pop_metadata", "snp_info", "K_wild")
missing <- required[!sapply(required, exists)]
if(length(missing) > 0) {
  stop("Missing objects: ", paste(missing, collapse = ", "))
}

# Check dimensions
dim(pop_metadata)      # Should have your sample size
dim(snp_info)          # Should match number of SNPs
table(pop_metadata$category)  # Verify wild_category exists
```

### Step 2: Verify Column Names

```r
# Check critical columns
names(pop_metadata)  # Should have: individual, category, latitude, longitude
names(snp_info)      # Should have: marker, chrom, pos

# Check genotype objects
ls(pattern = "geno_mat")  # Find your genotype matrices
ls(pattern = "geno_hmp")  # Find your hapmap objects
```

### Step 3: Test Data Matching

```r
# Load environmental data
load("data/germplasm_origin_bioclim_data.RData")

# Check overlap
sum(pop_metadata$individual %in% names(eaa_environmental_data))

# Should have reasonable overlap (>50% ideally)
```

---

## Running the Analysis

### Option 1: Full Analysis (Recommended)

1. Open `environmental_association_analysis.Rmd` in RStudio
2. Click "Knit" button (or Ctrl/Cmd + Shift + K)
3. Wait for completion (5-30 minutes depending on data size)
4. Review HTML output

### Option 2: Section-by-Section

Run sections in order:
1. ✅ Configuration
2. ✅ Data Loading  
3. ✅ Data Preparation
4. ✅ Model Selection
5. ✅ Population Genetic Comparisons
6. ✅ Gene Annotation (if GFF available)
7. ✅ Save Results

Check output after each section before proceeding.

---

## Expected Runtime

| Dataset Size | Runtime |
|--------------|---------|
| <10,000 SNPs, <100 samples | 2-5 min |
| 10-50K SNPs, 100-500 samples | 5-15 min |
| >50K SNPs, >500 samples | 15-60 min |

**Note:** Gene annotation adds 5-20 minutes

---

## Interpreting Your Results

### Key Questions to Ask

**1. How many significant associations?**
```r
# Check in output
n_total_assoc      # Total marker-variable pairs
n_unique_markers   # Unique SNPs
n_unique_vars      # Environmental variables
```

**Expected:** Depends on species and environment diversity
- Well-powered study: 50-500 associations
- Underpowered: <20 associations
- Very strong signal: >1000 associations

**2. Are eGWAS markers enriched for high F~ST~/SPA?**
```r
p_value_fst   # Permutation test p-value
p_value_spa   # Permutation test p-value
```

**Interpretation:**
- **p < 0.01:** Strong evidence of environmental selection
- **p = 0.01-0.05:** Moderate evidence
- **p > 0.05:** No clear enrichment (not necessarily bad!)

**3. How many eGWAS markers are also outliers?**
```r
n_egwas_fst_overlap   # Overlap with FST outliers
n_egwas_spa_overlap   # Overlap with SPA outliers
```

**Expected:** Usually 5-20% overlap

### Red Flags 🚩

- **No significant markers:** FDR too stringent or weak signal
- **All markers significant:** FDR too liberal or model misspecification  
- **Negative p-values from permutation:** Bug in code (check permutation setup)
- **p = 0 exactly:** Need more permutations or very strong signal

---

## Troubleshooting Guide

### Problem: "Object not found"

**Fix:**
```r
# Find actual object names
ls()  # List all objects
ls(pattern = "geno")  # Find genotype objects

# Manually assign
geno_mat_wild <- your_actual_name
```

### Problem: No significant markers

**Try:**
```r
# More liberal threshold
fdr_threshold <- 0.30

# Check p-value distribution
hist(gwas_out$pvalue)
min(gwas_out$pvalue)  # Best p-value
```

### Problem: GFF chromosomes don't match

**Fix:**
```r
# Check formats
unique(snp_info$chrom)       # e.g., "1", "2", "3"
unique(species_gff$seqid)    # e.g., "Chr1", "Chr2"

# Standardize
species_gff$seqid <- gsub("Chr", "", species_gff$seqid)
```

### Problem: Analysis runs out of memory

**Fix:**
```r
# Reduce permutations
n_permutations <- 1000

# Filter to significant markers only
# Clean up large objects
rm(large_object)
gc()
```

---

## Output Files

After successful completion:

✅ **`{species_name}_egwas_analysis_results.RData`**
- Contains all results for downstream analysis
- Load with: `load("results/YourSpecies_egwas_analysis_results.RData")`

✅ **`{species_name}_eaaSPA_versus_randomSPA.png`**  
- Figure comparing eGWAS markers to random background (SPA)

✅ **`{species_name}_eaaFST_versus_randomFST.png`**
- Figure comparing eGWAS markers to random background (F~ST~)

✅ **`environmental_association_analysis.html`**
- Complete report with all analysis and figures

---

## Next Steps

After running this analysis:

1. **Identify candidate genes** - Look at markers with nearby gene annotations
2. **Functional annotation** - GO enrichment, pathway analysis
3. **Validate top hits** - Independent population or experimental validation
4. **Geographic visualization** - Map allele frequencies across environment
5. **Genomic selection** - Use eGWAS markers for climate-adapted breeding

---

## Need More Help?

📖 **See full documentation:** `DOCUMENTATION_environmental_association_analysis.md`

🔧 **Check your data structure:**
```r
str(your_object)
summary(your_data)
```

💬 **Get R session info:**
```r
sessionInfo()
```

---

## Common Variations by Species

### Self-Pollinating Crops (Rice, Wheat, Barley, Lentil)
```r
ld_decay_distance <- 100000  # Extended LD
fdr_threshold <- 0.15        # Can be more stringent
```

### Outcrossing Crops (Maize, Sunflower, Alfalfa)
```r
ld_decay_distance <- 5000    # Rapid LD decay
fdr_threshold <- 0.20        # More liberal
outlier_percentile <- 0.995  # Top 0.5%
```

### Clonal Perennials (Potato, Sweet Potato, Yams)
```r
ld_decay_distance <- 500000  # Very extended LD
fdr_threshold <- 0.10        # Stringent
# May need separate analysis by clone group
```

### Tree Crops (Apple, Citrus, Coffee)
```r
ld_decay_distance <- 2000    # Rapid LD decay
wild_category <- "Wild"      # Or "Landrace" or "Traditional"
# Consider population structure carefully
```

---

## Quick Reference Card

| Parameter | Typical Range | Your Value |
|-----------|---------------|------------|
| Species name | Any | _____________ |
| LD decay (bp) | 2,000 - 500,000 | _____________ |
| FDR threshold | 0.05 - 0.30 | _____________ |
| Wild category | Species-specific | _____________ |
| Sample size | 50 - 1000+ | _____________ |
| SNP count | 5,000 - 500,000+ | _____________ |

---

**Ready to start?** Update the configuration section and run the analysis!

**Questions?** Consult the full documentation for detailed explanations.
