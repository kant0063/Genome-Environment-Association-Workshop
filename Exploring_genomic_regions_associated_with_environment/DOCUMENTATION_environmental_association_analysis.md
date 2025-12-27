# Environmental Association Analysis Workflow
## Documentation for Species Adaptation

---

## Table of Contents

1. [Overview](#overview)
2. [Required Data Inputs](#required-data-inputs)
3. [Species-Specific Configuration](#species-specific-configuration)
4. [Data Structure Requirements](#data-structure-requirements)
5. [Step-by-Step Adaptation Guide](#step-by-step-adaptation-guide)
6. [Troubleshooting](#troubleshooting)
7. [Expected Outputs](#expected-outputs)
8. [References and Resources](#references-and-resources)

---

## Overview

This workflow analyzes environmental association study (eGWAS) results to:

1. **Identify significant marker-environment associations** using FDR correction
2. **Compare eGWAS markers with population genetic outliers** (F~ST~ and spatial analysis [SPA])
3. **Test for enrichment** using permutation-based approaches
4. **Annotate markers** with nearby gene information (if GFF available)
5. **Generate publication-quality figures**

**Key Features:**
- Generic framework applicable to any plant species
- Flexible model selection strategies
- Comprehensive statistical testing
- Optional gene annotation analysis
- Reproducible R Markdown workflow

---

## Required Data Inputs

### 1. Population Metadata and Genotypes (`RData` file)

**Required Objects:**

| Object | Type | Description |
|--------|------|-------------|
| `pop_metadata` | data.frame | Individual metadata with columns: `individual`, `category`, `latitude`, `longitude` |
| `geno_mat_*_filter*` | matrix | Filtered genotype matrix (individuals × markers) |
| `geno_hmp_*_filter*` | data.frame | Hapmap-format genotypes with columns: `marker`, `chrom`, `pos`, `alleles`, individual columns |
| `snp_info` | data.frame | Marker information with columns: `marker`, `chrom`, `pos` |
| `K_wild` | matrix | Kinship/relationship matrix for wild individuals |

**Example Structure:**

```r
# pop_metadata
#   individual  category  latitude  longitude  population
#   IND001      Wild      45.2      -93.1      POP1
#   IND002      Cultivar  42.8      -88.5      POP2

# geno_mat_wild (numeric matrix, -1/0/1 coding)
#        SNP001  SNP002  SNP003
# IND001    -1       0       1
# IND002     0       1      -1

# geno_hmp_wild
#   marker  chrom  pos    alleles  IND001  IND002
#   SNP001  1      1250   A/G      -1      0
#   SNP002  1      5820   C/T      0       1

# snp_info
#   marker  chrom  pos
#   SNP001  1      1250
#   SNP002  1      5820
```

### 2. Environmental/Climate Data (`RData` file)

**Required Objects:**

| Object | Type | Description |
|--------|------|-------------|
| `eaa_environmental_data` | data.frame | Environmental variables for each location/individual |
| `eaa_environmental_vars` | data.frame | Metadata about environmental variables with columns: `variable`, `full_name`, `class` (optional) |

**Example Structure:**

```r
# eaa_environmental_data
#   location_abbr  bio1   bio2   bio12  elevation  ...
#   LOC001         15.2   8.3    850    245
#   LOC002         18.5   9.1    620    180

# eaa_environmental_vars
#   variable  full_name                      class
#   bio1      Annual Mean Temperature        temperature
#   bio2      Mean Diurnal Range            temperature  
#   bio12     Annual Precipitation          precipitation
#   elevation Elevation                     geography
```

### 3. eGWAS Results (`RData` file)

**Required Objects:**

| Object | Type | Description |
|--------|------|-------------|
| `gwas_out` | data.frame | GWAS results with columns: `marker`, `variable`, `model`, `score`, `pvalue` |
| `p_lambda` | data.frame | Genomic inflation factors with columns: `model`, `variable`, `lambda`, `best_model` (optional) |
| `eaa_gwas_sigmar` | data.frame | Pre-computed significant markers (optional; can be calculated from `gwas_out`) |

**Example Structure:**

```r
# gwas_out
#   marker  variable  model   score     pvalue
#   SNP001  bio1      model1  15.32     0.0001
#   SNP001  bio1      model2  18.45     0.00002
#   SNP002  bio2      model1  8.21      0.004

# p_lambda (optional)
#   model   variable  lambda  best_model
#   model1  bio1      1.15    
#   model2  bio1      1.02    *
#   model1  bio2      1.08    *
```

### 4. Population Genetics Results (`RData` file)

**Required Objects:**

| Object | Type | Description |
|--------|------|-------------|
| `global_marker_fst` | data.frame | F~ST~ values with columns: `marker`, `Fst` |
| `spa_out` | data.frame | Spatial analysis scores with columns: `marker`, `spa_score` |

**Example Structure:**

```r
# global_marker_fst
#   marker  Fst
#   SNP001  0.125
#   SNP002  0.089
#   SNP003  0.234

# spa_out
#   marker  spa_score
#   SNP001  0.456
#   SNP002  0.312
#   SNP003  0.789
```

### 5. GFF Annotation File (Optional)

- Standard GFF3 format
- Must contain genes and other features
- Chromosome IDs must match those in `snp_info$chrom`

**Example:**

```
chr1  source  gene  1000  5000  .  +  .  ID=gene001;Name=ABC1
chr1  source  mRNA  1000  5000  .  +  .  ID=mRNA001;Parent=gene001
```

---

## Species-Specific Configuration

### Essential Parameters to Update

```r
# 1. SPECIES IDENTIFIER
species_name <- "YourSpecies"  
# Examples: "Pepper", "Lentil", "Barley", "Yam"

# 2. LD DECAY DISTANCE
ld_decay_distance <- 17000  # bp
# HOW TO DETERMINE:
# - Run LD decay analysis on your population
# - Find distance where r² drops to ~0.2
# - Literature values for similar species
# - Conservative estimate: 10-50 kb for most crops

# 3. FDR THRESHOLD
fdr_threshold <- 0.20
# Adjust based on:
# - Study design and statistical power
# - Number of environmental variables tested
# - Typical values: 0.05 (strict), 0.10, 0.20 (liberal)

# 4. WILD POPULATION CATEGORY
wild_category <- "Wild"
# Must match your metadata exactly
# Common alternatives: "Landrace", "Native", "Unimproved"

# 5. MODEL SELECTION
model_selection <- "best_lambda"
# Options:
# - "best_lambda": Use model with best genomic inflation factor
# - "manual": Specify a single model (e.g., "model2")
# - "all_models": Use all models (increases multiple testing)

# 6. GFF FILE
gff_file_name <- "species_annotations.gff"
# Set to NULL if not available
```

### Optional Parameters

```r
# Outlier threshold (for SPA and FST)
outlier_percentile <- 0.999  # Top 0.1%

# Number of permutations
n_permutations <- 10000  # More = better p-value precision

# Seed for reproducibility
set.seed(209)  # Use consistent seed across analyses
```

---

## Data Structure Requirements

### Critical Column Names

The workflow expects specific column names. If your data uses different names, rename them before running:

```r
# Example: Renaming columns to match expected format
pop_metadata <- your_metadata %>%
  rename(
    individual = sample_id,
    category = germplasm_type,
    latitude = lat,
    longitude = long
  )

snp_info <- your_snp_info %>%
  rename(
    marker = snp_id,
    chrom = chromosome,
    pos = position
  )
```

### Genotype Coding

**Expected Format:** Additive coding (-1, 0, 1) where:
- `-1` = Homozygous reference (0 copies of alternate allele)
- `0` = Heterozygous (1 copy of alternate allele)  
- `1` = Homozygous alternate (2 copies of alternate allele)

**If your data uses different coding:**

```r
# Example: Convert 0/1/2 to -1/0/1
geno_mat <- geno_mat - 1

# Example: Convert AA/AB/BB to -1/0/1
geno_mat <- case_when(
  geno_mat == "AA" ~ -1,
  geno_mat == "AB" ~ 0,
  geno_mat == "BB" ~ 1
)
```

### Missing Data

- Missing genotypes should be `NA`
- The workflow handles missing data automatically
- High missingness (>20%) may affect results

---

## Step-by-Step Adaptation Guide

### Step 1: Prepare Your Data

1. **Organize data files:**
   ```
   project/
   ├── data/
   │   ├── population_metadata_and_genotypes.RData
   │   ├── germplasm_origin_bioclim_data.RData
   │   ├── species_annotations.gff (optional)
   ├── results/
   │   ├── eaa_gwas_results.RData
   │   ├── population_genetics_stats.RData
   ├── figures/
   └── environmental_association_analysis.Rmd
   ```

2. **Verify required objects exist:**
   ```r
   # Load each RData file and check
   load("data/population_metadata_and_genotypes.RData")
   ls()  # Should show: pop_metadata, geno_mat_*, snp_info, K_wild
   
   # Check dimensions
   dim(pop_metadata)
   dim(snp_info)
   table(pop_metadata$category)  # Verify your wild_category exists
   ```

3. **Standardize column names:**
   ```r
   # Use the exact column names expected by workflow
   # See "Critical Column Names" section above
   ```

### Step 2: Configure Parameters

1. Open `environmental_association_analysis.Rmd`

2. Update the **Configuration** section:
   ```r
   species_name <- "YourSpecies"
   fdr_threshold <- 0.20  # Adjust as needed
   ld_decay_distance <- 15000  # YOUR VALUE HERE
   wild_category <- "Wild"  # Match your data
   model_selection <- "best_lambda"
   gff_file_name <- "your_species.gff"  # or NULL
   ```

3. Update **File Paths** if different from defaults:
   ```r
   data_dir <- "path/to/your/data"
   result_dir <- "path/to/results"
   fig_dir <- "path/to/figures"
   ```

### Step 3: Run Data Checks

Run the **Data Loading** and **Data Structure Check** sections first:

```r
# Should complete without errors
# Check output for:
# - Correct number of SNPs and individuals
# - All required objects present
# - Reasonable dimensions
```

### Step 4: Execute Analysis

Run sections sequentially:

1. **Data Preparation** - Prepares genotype matrices
2. **Model Selection** - Identifies significant markers
3. **Population Genetic Comparisons** - SPA and FST analyses  
4. **Gene Annotation** - If GFF available
5. **Save Results** - Outputs RData and figures

### Step 5: Interpret Results

Key results to examine:

1. **Number of significant associations:**
   - How many marker-variable pairs?
   - How many unique markers?
   - Distribution across environmental variables

2. **Permutation test p-values:**
   - Are eGWAS markers enriched for high SPA scores?
   - Are eGWAS markers enriched for high F~ST~ values?
   - p < 0.05 suggests true signal vs. random background

3. **Gene proximity:**
   - Are significant markers near genes?
   - Any enrichment in genic vs. intergenic regions?

---

## Troubleshooting

### Common Issues and Solutions

#### Issue: "Object not found" errors

**Cause:** Variable naming doesn't match expected patterns

**Solution:**
```r
# Check what genotype objects exist
ls(pattern = "geno_mat")
ls(pattern = "geno_hmp")

# Manually assign if needed
geno_mat_wild <- geno_mat_wild_final  # Use your actual name
geno_hmp_wild <- geno_hmp_wild_final
```

#### Issue: No significant markers found

**Possible causes:**
- FDR threshold too stringent
- Weak environmental associations
- Small sample size

**Solutions:**
```r
# Try more liberal FDR
fdr_threshold <- 0.30

# Check distribution of p-values
hist(gwas_out$pvalue, breaks = 50)

# Check if any markers approach significance
gwas_out %>%
  group_by(variable) %>%
  summarize(min_p = min(pvalue))
```

#### Issue: GFF chromosome IDs don't match

**Cause:** Different chromosome naming conventions

**Solution:**
```r
# Check formats
unique(snp_info$chrom)  # e.g., "1", "2", "3"
unique(species_gff$seqid)  # e.g., "chr1", "chr2"

# Convert if needed
species_gff <- species_gff %>%
  mutate(seqid = str_remove(seqid, "chr"))

# Or convert SNP info
snp_info <- snp_info %>%
  mutate(chrom = paste0("chr", chrom))
```

#### Issue: Permutation test runs very slowly

**Cause:** Large number of markers and/or permutations

**Solutions:**
```r
# Reduce permutations for testing
n_permutations <- 1000  # Instead of 10000

# Use parallel processing (advanced)
library(parallel)
cl <- makeCluster(4)  # Use 4 cores
# Modify permutation code to use parReplicate()
```

#### Issue: Memory errors with large datasets

**Solutions:**
```r
# Process chromosomes separately
# Filter to relevant markers only
geno_mat_wild1 <- geno_mat_wild[, colnames(geno_mat_wild) %in% egwas_sigmar$marker]

# Use data.table for large data
library(data.table)
gwas_out <- as.data.table(gwas_out)
```

---

## Expected Outputs

### 1. RData File

**Filename:** `{species_name}_egwas_analysis_results.RData`

**Contains:**
- `egwas_sigmar`: Significant marker-variable associations
- `eaa_gwas_spa_scores`: SPA scores for eGWAS markers
- `eaa_gwas_fst`: F~ST~ values for eGWAS markers
- `spa_outliers`: Top SPA score outliers
- `fst_outliers`: Top F~ST~ outliers
- `p_value_spa`: Permutation test p-value (SPA)
- `p_value_fst`: Permutation test p-value (F~ST~)
- `all_sigmar_nearby_annotation`: Gene annotations (if GFF provided)
- `all_markers_nearby_annotation`: Background gene annotations

### 2. Figures

**Figure 1:** `{species_name}_eaaSPA_versus_randomSPA.png`
- Histogram of mean SPA scores from random marker sets
- Vertical line showing mean SPA for eGWAS markers
- Demonstrates if eGWAS markers have elevated SPA scores

**Figure 2:** `{species_name}_eaaFST_versus_randomFST.png`
- Histogram of mean F~ST~ from random marker sets
- Vertical line showing mean F~ST~ for eGWAS markers
- Demonstrates if eGWAS markers have elevated F~ST~

**In-document figures:**
- Significant markers per environmental variable (bar plot)
- Significant markers by environmental class (bar plot)
- Distribution of associations per marker (histogram)

### 3. HTML Report

The R Markdown generates a comprehensive HTML report with:
- Interactive table of contents
- All code and output
- Embedded figures
- Summary statistics
- Session information for reproducibility

---

## Interpretation Guidelines

### Understanding Results

**1. Permutation Test P-values**

- **p < 0.01:** Strong evidence that eGWAS markers are enriched for SPA/F~ST~
- **0.01 < p < 0.05:** Moderate evidence of enrichment
- **p > 0.05:** No evidence of enrichment (could still be real biological signal)

**What this means:**
- Low p-values suggest environmental selection is driving population differentiation
- High p-values don't invalidate eGWAS results—just indicate selection isn't the primary driver

**2. Overlap with Outliers**

- **Few overlaps:** Environmental adaptation may involve subtle allele frequency shifts
- **Many overlaps:** Strong selective sweeps in response to environment
- **No overlap:** Doesn't invalidate eGWAS—different biological processes

**3. Gene Annotation Results**

- **Enrichment in genic regions:** Suggests functional variants
- **No enrichment:** Could indicate regulatory variants or incomplete annotation
- **Specific genes:** Candidates for functional validation

### Biological Interpretation

**Elevated SPA scores suggest:**
- Spatially structured selection
- Geographic isolation with gene flow
- Local adaptation to environment

**Elevated F~ST~ suggests:**
- Strong divergent selection
- Ancient population splits
- Limited gene flow between populations

**Both elevated:**
- Strong evidence for environmental adaptation
- Likely candidates for functional studies

---

## Species-Specific Considerations

### Clonal/Asexual Species
- LD may extend much further (use larger `ld_decay_distance`)
- Population structure may be extreme
- Consider additional corrections for clonal structure

### Polyploid Species
- Genotype coding may need adjustment
- Allele dosage effects important
- Consider ploidy level in interpretation

### Mixed Mating Systems
- Expect intermediate LD patterns
- May need separate analyses for selfing vs. outcrossing lineages

### Domesticated vs. Wild
- Different LD patterns
- May need separate analyses by germplasm type
- Selection signatures may differ

---

## Advanced Customization

### Adding Custom Environmental Variables

```r
# Add your own environmental data
my_env_data <- read.csv("custom_environmental_data.csv")

# Merge with existing data
eaa_environmental_data <- eaa_environmental_data %>%
  left_join(my_env_data, by = "location_abbr")

# Update variable list
new_vars <- tibble(
  variable = c("my_var1", "my_var2"),
  full_name = c("My Variable 1", "My Variable 2"),
  class = c("custom", "custom")
)

eaa_environmental_vars <- bind_rows(eaa_environmental_vars, new_vars)
```

### Custom Plotting Themes

```r
# Define your own theme
theme_myspecies <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
}

# Use in plots
ggplot(...) + theme_myspecies()
```

### Subsetting to Specific Chromosomes

```r
# Analyze only certain chromosomes
chromosomes_to_analyze <- c("1", "2", "3")

egwas_sigmar_subset <- egwas_sigmar %>%
  filter(chrom %in% chromosomes_to_analyze)

# Continue analysis with subset
```

---

## Quality Control Checklist

Before finalizing results:

- [ ] All required data loaded successfully
- [ ] Genotype matrices have correct dimensions
- [ ] Environmental data matches genetic samples
- [ ] Column names standardized
- [ ] No excessive missing data (>20%)
- [ ] Significant markers found (if expected)
- [ ] Permutation tests completed
- [ ] Figures generated correctly
- [ ] Results saved to output directory
- [ ] HTML report renders without errors

---

## References and Resources

### Key Papers

1. **Environmental Association Analysis:**
   - Rellstab et al. (2015). "A practical guide to environmental association analysis in landscape genomics." Mol Ecol.
   - Forester et al. (2018). "Comparing methods for detecting multilocus adaptation." Mol Ecol.

2. **Spatial Analysis:**
   - Duforet-Frebourg et al. (2016). "Detecting genomic signatures of natural selection with principal component analysis." Mol Biol Evol.

3. **F~ST~ Outlier Detection:**
   - Whitlock & Lotterhos (2015). "Reliable detection of loci responsible for local adaptation." Genetics.

### Software Documentation

- **sommer:** https://cran.r-project.org/package=sommer
- **qvalue:** http://bioconductor.org/packages/qvalue/
- **ape (GFF reading):** https://cran.r-project.org/package=ape

### Related Workflows

If using this workflow, you likely also need:
1. **Population structure analysis** (PCA, ADMIXTURE)
2. **LD decay estimation**
3. **Environmental GWAS** (generates input for this workflow)
4. **F~ST~ calculation**
5. **Spatial analysis** (SPA/PCAdapt)

---

## Getting Help

### Common Resources

1. **Check object structure:**
   ```r
   str(your_object)
   head(your_object)
   ```

2. **Verify data integrity:**
   ```r
   summary(your_data)
   table(is.na(your_data))
   ```

3. **Test with small subset:**
   ```r
   test_data <- your_data[1:100, ]
   # Run workflow on test_data first
   ```

### Reporting Issues

When asking for help, include:
- Species and study system
- Dimensions of your data
- Specific error messages
- Output from `sessionInfo()`
- Relevant code chunk

---

## Citation

If you use this workflow, please cite:

- The original cranberry eGWAS study (if applicable)
- R packages used (see Session Info)
- This workflow (if publicly available)

---

## Version History

- **v1.0** (2024): Initial generic version adapted from cranberry-specific analysis
- Designed for maximum flexibility across plant species
- Comprehensive documentation and error handling

---

**Document prepared for:** Generic plant species environmental association analysis  
**Last updated:** December 2024  
**Maintained by:** Michael Price-Waldman, University of Hawaii
