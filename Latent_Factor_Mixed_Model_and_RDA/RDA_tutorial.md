# Redundancy Analysis (RDA) Tutorial for Genomic Data

## Overview

Redundancy Analysis (RDA) is a constrained ordination method that identifies associations between genomic variation (from SNP data) and environmental or climate variables. This tutorial provides a complete workflow for performing RDA with genomic data.

## Table of Contents

1. [Introduction](#introduction)
2. [Required R Packages](#required-r-packages)
3. [Input Data Requirements](#input-data-requirements)
4. [Analysis Workflow](#analysis-workflow)
5. [Interpreting Results](#interpreting-results)
6. [Complete R Script](#complete-r-script)

---

## Introduction

### What is RDA?

Redundancy Analysis is a multivariate method that:
- Identifies relationships between genomic variation and environmental predictors
- Detects candidate loci potentially under environmental selection
- Provides visualization of genotype-environment associations
- Tests statistical significance through permutation tests

### When to Use RDA

Use RDA when you have:
- Genomic data (SNPs) from multiple individuals/populations
- Environmental or climate data for sample collection locations
- Interest in identifying adaptive genetic variation
- Questions about genotype-environment associations

---

## Required R Packages

```r
# Install and load required packages
required_packages <- c("vcfR", "vegan", "ggplot2", "dplyr", "corrplot")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
```

**Package purposes:**
- `vcfR`: Reading and handling VCF files
- `vegan`: Performing RDA and multivariate analyses
- `ggplot2`: Creating publication-quality plots
- `dplyr`: Data manipulation
- `corrplot`: Visualizing correlation matrices

---

## Input Data Requirements

### 1. VCF File (Genotype Data)

Your VCF file should contain:
- SNP genotype calls for all samples
- Chromosome and position information
- Quality-filtered SNPs

**Example VCF structure:**
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3
Chr1    1000    .       A       G       99      PASS    .       GT      0/0     0/1     1/1
Chr1    2000    .       C       T       99      PASS    .       GT      0/1     0/1     0/0
```

### 2. Environmental/Climate CSV File

The CSV file must include:
- **First column**: Sample IDs (must match VCF sample names)
- **Optional columns**: Longitude, Latitude (will be excluded from analysis)
- **Remaining columns**: Environmental/climate variables (continuous)

**Example CSV structure:**
```
SampleID,Longitude,Latitude,Temperature,Precipitation,Elevation
Sample1,-120.5,45.2,15.3,800,1200
Sample2,-118.3,44.8,16.1,750,950
Sample3,-122.1,46.5,14.2,900,1450
```

---

## Analysis Workflow

### Step 1: Set Analysis Parameters

```r
# === FILE INPUTS ===
VCF_FILE <- "your_genotypes.vcf"
CLIMATE_CSV <- "your_environmental_data.csv"
SAMPLE_ID_COLUMN <- "SampleID"
EXCLUDE_COLUMNS <- c("SampleID", "Longitude", "Latitude")

# === ANALYSIS PARAMETERS ===
CORRELATION_CUTOFF <- 0.9        # Remove highly correlated variables
OUTLIER_QUANTILE <- 0.999        # 99.9th percentile for candidate SNPs
N_PERMUTATIONS <- 100            # Permutations for significance testing

# === OUTPUT OPTIONS ===
OUTPUT_PREFIX <- "rda"
OUTPUT_DIR <- "rda_results"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}
```

### Step 2: Read and Process VCF Data

```r
# Read VCF file
vcf <- read.vcfR(VCF_FILE)
cat("Loaded", nrow(vcf@fix), "SNPs from", ncol(vcf@gt) - 1, "samples\n")

# Extract numeric genotype matrix
gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
```

### Step 3: Filter SNPs by Missing Data

```r
# Calculate missing data proportion per SNP
missing_per_snp <- rowSums(is.na(gt_numeric))
missing_proportion <- missing_per_snp / ncol(gt_numeric)
max_missing_allowed <- 0.30  # 30% threshold

# Filter SNPs
snps_pass_filter <- missing_proportion <= max_missing_allowed
gt_numeric <- gt_numeric[snps_pass_filter, , drop = FALSE]
vcf@fix <- vcf@fix[snps_pass_filter, , drop = FALSE]
vcf@gt <- vcf@gt[snps_pass_filter, , drop = FALSE]

cat("Retained", nrow(gt_numeric), "SNPs after filtering\n")
```

### Step 4: Read and Match Environmental Data

```r
# Read environmental data
env_data <- read.csv(CLIMATE_CSV, header = TRUE, stringsAsFactors = FALSE)

# Match samples between VCF and CSV
vcf_samples <- colnames(gt_numeric)
csv_samples <- env_data[[SAMPLE_ID_COLUMN]]

# Find common samples
common_samples <- intersect(vcf_samples, csv_samples)
cat("Found", length(common_samples), "samples with both genetic and environmental data\n")

# Subset data to common samples
gt_matched <- gt_numeric[, common_samples, drop = FALSE]
env_matched <- env_data[match(common_samples, csv_samples), ]
```

### Step 5: Prepare Environmental Variables

```r
# Extract environmental variables
env_vars <- env_matched[, !colnames(env_matched) %in% EXCLUDE_COLUMNS, drop = FALSE]

# Convert to numeric and handle missing data
env_vars <- as.data.frame(lapply(env_vars, as.numeric))
rownames(env_vars) <- common_samples

# Remove variables with too much missing data
missing_threshold <- 0.1  # 10%
env_vars <- env_vars[, colSums(is.na(env_vars)) / nrow(env_vars) <= missing_threshold, drop = FALSE]
```

### Step 6: Handle Multicollinearity

Highly correlated environmental variables can cause problems in RDA. Remove redundant variables:

```r
# Function to identify highly correlated variables
identify_highly_correlated <- function(cor_matrix, cutoff = 0.9) {
  upper_tri <- upper.tri(cor_matrix)
  high_cor_pairs <- which(abs(cor_matrix) > cutoff & upper_tri, arr.ind = TRUE)
  
  if (nrow(high_cor_pairs) == 0) {
    return(integer(0))
  }
  
  vars_to_remove <- c()
  
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- high_cor_pairs[i, 1]
    var2 <- high_cor_pairs[i, 2]
    
    mean_cor_var1 <- mean(abs(cor_matrix[var1, -var1]), na.rm = TRUE)
    mean_cor_var2 <- mean(abs(cor_matrix[var2, -var2]), na.rm = TRUE)
    
    if (mean_cor_var1 > mean_cor_var2) {
      vars_to_remove <- c(vars_to_remove, var1)
    } else {
      vars_to_remove <- c(vars_to_remove, var2)
    }
  }
  
  return(unique(vars_to_remove))
}

# Calculate correlation matrix
env_cor <- cor(env_vars, use = "complete.obs")

# Identify and remove highly correlated variables
vars_to_remove <- identify_highly_correlated(env_cor, CORRELATION_CUTOFF)
if (length(vars_to_remove) > 0) {
  env_vars_final <- env_vars[, -vars_to_remove, drop = FALSE]
} else {
  env_vars_final <- env_vars
}

cat("Retained", ncol(env_vars_final), "environmental variables after removing correlated ones\n")
```

### Step 7: Scale and Prepare Data for RDA

```r
# Scale environmental variables (mean = 0, sd = 1)
env_vars_scaled <- scale(env_vars_final, center = TRUE, scale = TRUE)

# Prepare genotype data (transpose and impute missing values)
genotype_data <- t(gt_matched)
genotype_data[is.na(genotype_data)] <- mean(genotype_data, na.rm = TRUE)

# Ensure sample order matches
genotype_data <- genotype_data[match(rownames(env_vars_scaled), rownames(genotype_data)), ]
```

### Step 8: Run RDA

```r
# Perform RDA
rda_result <- rda(genotype_data ~ ., data = as.data.frame(env_vars_scaled))

# Print summary
cat("\n=== RDA Summary ===\n")
summary(rda_result)
```

### Step 9: Test Significance

```r
# Overall model significance
anova_result <- anova.cca(rda_result, permutations = N_PERMUTATIONS)
cat("\n=== Overall Model Significance ===\n")
print(anova_result)

# Individual axis significance
anova_axes <- anova.cca(rda_result, by = "axis", permutations = N_PERMUTATIONS)
cat("\n=== Axis Significance ===\n")
print(anova_axes)

# Individual variable significance
anova_terms <- anova.cca(rda_result, by = "terms", permutations = N_PERMUTATIONS)
cat("\n=== Environmental Variable Significance ===\n")
print(anova_terms)
```

### Step 10: Extract Variance Explained

```r
# Calculate variance components
total_var <- rda_result$tot.chi
constrained_var <- rda_result$CCA$tot.chi
unconstrained_var <- rda_result$CA$tot.chi
var_explained <- constrained_var / total_var

# Variance per axis
eigenvalues <- rda_result$CCA$eig
axis_var <- (eigenvalues / total_var) * 100

cat("\nVariance Explained:\n")
cat("  - Total variance:", round(total_var, 2), "\n")
cat("  - Constrained variance:", round(constrained_var, 2), "\n")
cat("  - Proportion explained:", round(var_explained * 100, 2), "%\n")
cat("  - RDA1:", round(axis_var[1], 2), "%\n")
cat("  - RDA2:", round(axis_var[2], 2), "%\n")
```

### Step 11: Create RDA Biplot

```r
# Extract sample scores
sample_scores <- scores(rda_result, choices = 1:2, display = "sites", scaling = 1)
rda_df <- data.frame(
  Sample = rownames(sample_scores),
  RDA1 = sample_scores[, 1],
  RDA2 = sample_scores[, 2]
)

# Extract environmental variable loadings
env_scores <- scores(rda_result, choices = 1:2, display = "bp", scaling = 1)
env_arrows <- data.frame(
  Variable = rownames(env_scores),
  RDA1 = env_scores[, 1],
  RDA2 = env_scores[, 2]
)

# Create biplot
library(ggplot2)

rda_plot <- ggplot(rda_df, aes(x = RDA1, y = RDA2)) +
  geom_point(size = 2.5, alpha = 0.7, color = "steelblue") +
  geom_segment(data = env_arrows, 
               aes(x = 0, y = 0, xend = RDA1 * 3, yend = RDA2 * 3),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "red", size = 1, alpha = 0.8) +
  geom_text(data = env_arrows, 
            aes(x = RDA1 * 3.2, y = RDA2 * 3.2, label = Variable),
            color = "black", size = 3, fontface = "bold") +
  labs(
    title = "RDA Biplot: Genomic Variation Constrained by Environment",
    subtitle = paste0("Total variance explained: ", round(var_explained * 100, 2), "%"),
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  coord_fixed()

print(rda_plot)
```

### Step 12: Identify Candidate SNPs

```r
# Extract SNP loadings
species_scores <- scores(rda_result, choices = 1:2, display = "species", scaling = 1)

# Calculate Mahalanobis distances
loadings_subset <- species_scores[, 1:2]
maha_dist <- mahalanobis(loadings_subset, 
                         center = colMeans(loadings_subset), 
                         cov = cov(loadings_subset))

# Identify outliers
outlier_threshold <- quantile(maha_dist, OUTLIER_QUANTILE)
outlier_snps <- which(maha_dist > outlier_threshold)

cat("Identified", length(outlier_snps), "candidate SNPs\n")

# Create outlier data frame
outlier_df <- data.frame(
  SNP_index = outlier_snps,
  RDA1_loading = species_scores[outlier_snps, 1],
  RDA2_loading = species_scores[outlier_snps, 2],
  Mahalanobis_distance = maha_dist[outlier_snps],
  chromosome = vcf@fix[outlier_snps, "CHROM"],
  position = as.numeric(vcf@fix[outlier_snps, "POS"])
)

# Plot SNP loadings
loading_df <- data.frame(
  SNP = 1:nrow(species_scores),
  RDA1 = species_scores[, 1],
  RDA2 = species_scores[, 2],
  Outlier = 1:nrow(species_scores) %in% outlier_snps
)

loading_plot <- ggplot(loading_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Outlier), size = 1, alpha = 0.6) +
  scale_color_manual(values = c("gray70", "red"), 
                     labels = c("Background SNPs", "Candidate SNPs")) +
  labs(
    title = "SNP Loadings on RDA Axes",
    subtitle = paste("Candidate SNPs (n =", length(outlier_snps), ") highlighted in red"),
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  coord_fixed()

print(loading_plot)
```

### Step 13: Save Results

```r
# Save RData object
save(rda_result, file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_results.RData")))

# Save plots
ggsave(file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_biplot.png")), 
       rda_plot, width = 10, height = 8, dpi = 300)

ggsave(file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_snp_loadings.png")), 
       loading_plot, width = 10, height = 8, dpi = 300)

# Save candidate SNPs
write.csv(outlier_df, 
          file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_outlier_snps.csv")), 
          row.names = FALSE)

# Create and save summary table
rda_summary_table <- data.frame(
  Metric = c("Total Variance", "Constrained Variance", "Unconstrained Variance", 
             "Proportion Explained", "Number of Environmental Variables", 
             "Number of Samples", "Number of SNPs", "Candidate SNPs"),
  Value = c(
    round(total_var, 2),
    round(constrained_var, 2), 
    round(unconstrained_var, 2),
    paste0(round(var_explained * 100, 2), "%"),
    ncol(env_vars_scaled),
    nrow(genotype_data),
    ncol(genotype_data),
    length(outlier_snps)
  )
)

write.csv(rda_summary_table, 
          file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_summary_table.csv")), 
          row.names = FALSE)
```

---

## Interpreting Results

### 1. Variance Explained

The **proportion of constrained variance** indicates how much of the total genomic variation is explained by environmental variables.

- **< 5%**: Weak environmental association
- **5-15%**: Moderate environmental association
- **> 15%**: Strong environmental association

### 2. Significance Testing

- **Overall model p-value**: Tests whether environmental variables significantly explain genomic variation
- **Axis p-values**: Tests significance of individual RDA axes
- **Variable p-values**: Identifies which environmental variables significantly contribute

**Interpretation:**
- p < 0.05: Statistically significant
- p < 0.01: Highly significant
- p > 0.05: Not significant

### 3. RDA Biplot

The biplot shows:
- **Points**: Individual samples positioned by genomic variation
- **Arrows**: Environmental variable gradients
  - Arrow length = strength of relationship
  - Arrow direction = gradient direction
  - Samples in arrow direction have higher values for that variable

**Reading the biplot:**
- Samples close together are genetically similar given environmental constraints
- Long arrows indicate strong environmental gradients
- Arrow angles show correlations between environmental variables
  - Acute angles (< 90°): Positive correlation
  - Obtuse angles (> 90°): Negative correlation
  - Right angles (90°): No correlation

### 4. Candidate SNPs

Candidate SNPs (outliers) are those with:
- Extreme loadings on RDA axes
- High Mahalanobis distances
- Potential adaptive significance

**Follow-up analyses:**
- Annotate candidate SNPs (genes, functional regions)
- Validate with independent datasets
- Test for signatures of selection
- Investigate biological functions

---

## Complete R Script

Below is the complete, ready-to-run script:

```r
#!/usr/bin/env Rscript

################################################################################
# Generic RDA (Redundancy Analysis) Script
# 
# This script performs RDA to identify associations between genomic variation
# (from VCF) and environmental/climate variables (from CSV).
################################################################################

# Load required libraries
required_packages <- c("vcfR", "vegan", "ggplot2", "dplyr", "corrplot")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

################################################################################
# USER-DEFINED PARAMETERS
################################################################################

# === FILE INPUTS ===
VCF_FILE <- "your_genotypes.vcf"
CLIMATE_CSV <- "your_environmental_data.csv"
SAMPLE_ID_COLUMN <- "SampleID"
EXCLUDE_COLUMNS <- c("SampleID", "Longitude", "Latitude")

# === ANALYSIS PARAMETERS ===
CORRELATION_CUTOFF <- 0.9
OUTLIER_QUANTILE <- 0.999
N_PERMUTATIONS <- 100
MAX_MISSING_DATA <- 0.30  # Maximum 30% missing data per SNP

# === OUTPUT OPTIONS ===
OUTPUT_PREFIX <- "rda"
OUTPUT_DIR <- "rda_results"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

################################################################################
# HELPER FUNCTIONS
################################################################################

identify_highly_correlated <- function(cor_matrix, cutoff = 0.9) {
  upper_tri <- upper.tri(cor_matrix)
  high_cor_pairs <- which(abs(cor_matrix) > cutoff & upper_tri, arr.ind = TRUE)
  
  if (nrow(high_cor_pairs) == 0) {
    return(integer(0))
  }
  
  vars_to_remove <- c()
  
  for (i in 1:nrow(high_cor_pairs)) {
    var1 <- high_cor_pairs[i, 1]
    var2 <- high_cor_pairs[i, 2]
    
    mean_cor_var1 <- mean(abs(cor_matrix[var1, -var1]), na.rm = TRUE)
    mean_cor_var2 <- mean(abs(cor_matrix[var2, -var2]), na.rm = TRUE)
    
    if (mean_cor_var1 > mean_cor_var2) {
      vars_to_remove <- c(vars_to_remove, var1)
    } else {
      vars_to_remove <- c(vars_to_remove, var2)
    }
  }
  
  return(unique(vars_to_remove))
}

################################################################################
# MAIN ANALYSIS
################################################################################

cat("=== Starting RDA Analysis ===\n\n")

# Step 1: Read VCF file
cat("Step 1: Reading VCF file...\n")
vcf <- read.vcfR(VCF_FILE)
cat("  - Loaded", nrow(vcf@fix), "SNPs from", ncol(vcf@gt) - 1, "samples\n")

# Step 2: Extract genotype matrix
cat("\nStep 2: Extracting genotype data...\n")
gt_numeric <- extract.gt(vcf, element = "GT", as.numeric = TRUE)

# Step 3: Filter SNPs by missing data
cat("\nStep 3: Filtering SNPs by missing data...\n")
missing_per_snp <- rowSums(is.na(gt_numeric))
missing_proportion <- missing_per_snp / ncol(gt_numeric)
snps_pass_filter <- missing_proportion <= MAX_MISSING_DATA

gt_numeric <- gt_numeric[snps_pass_filter, , drop = FALSE]
vcf@fix <- vcf@fix[snps_pass_filter, , drop = FALSE]
vcf@gt <- vcf@gt[snps_pass_filter, , drop = FALSE]

cat("  - Retained", nrow(gt_numeric), "SNPs after filtering\n")

# Step 4: Read environmental data
cat("\nStep 4: Reading environmental data...\n")
env_data <- read.csv(CLIMATE_CSV, header = TRUE, stringsAsFactors = FALSE)
cat("  - Loaded data for", nrow(env_data), "samples\n")

# Step 5: Match samples
cat("\nStep 5: Matching samples...\n")
vcf_samples <- colnames(gt_numeric)
csv_samples <- env_data[[SAMPLE_ID_COLUMN]]
common_samples <- intersect(vcf_samples, csv_samples)

cat("  - Found", length(common_samples), "common samples\n")

gt_matched <- gt_numeric[, common_samples, drop = FALSE]
env_matched <- env_data[match(common_samples, csv_samples), ]

# Step 6: Prepare environmental variables
cat("\nStep 6: Preparing environmental variables...\n")
env_vars <- env_matched[, !colnames(env_matched) %in% EXCLUDE_COLUMNS, drop = FALSE]
env_vars <- as.data.frame(lapply(env_vars, as.numeric))
rownames(env_vars) <- common_samples

# Remove variables with too much missing data
env_vars <- env_vars[, colSums(is.na(env_vars)) / nrow(env_vars) <= 0.1, drop = FALSE]

# Step 7: Handle multicollinearity
cat("\nStep 7: Checking for multicollinearity...\n")
env_cor <- cor(env_vars, use = "complete.obs")
vars_to_remove <- identify_highly_correlated(env_cor, CORRELATION_CUTOFF)

if (length(vars_to_remove) > 0) {
  env_vars_final <- env_vars[, -vars_to_remove, drop = FALSE]
  cat("  - Removed", length(vars_to_remove), "highly correlated variables\n")
} else {
  env_vars_final <- env_vars
}

cat("  - Retained", ncol(env_vars_final), "environmental variables\n")

# Step 8: Scale data
cat("\nStep 8: Scaling data...\n")
env_vars_scaled <- scale(env_vars_final, center = TRUE, scale = TRUE)
genotype_data <- t(gt_matched)
genotype_data[is.na(genotype_data)] <- mean(genotype_data, na.rm = TRUE)
genotype_data <- genotype_data[match(rownames(env_vars_scaled), rownames(genotype_data)), ]

# Step 9: Run RDA
cat("\nStep 9: Running RDA...\n")
rda_result <- rda(genotype_data ~ ., data = as.data.frame(env_vars_scaled))

# Step 10: Test significance
cat("\nStep 10: Testing significance...\n")
anova_overall <- anova.cca(rda_result, permutations = N_PERMUTATIONS)
anova_axes <- anova.cca(rda_result, by = "axis", permutations = N_PERMUTATIONS)
anova_terms <- anova.cca(rda_result, by = "terms", permutations = N_PERMUTATIONS)

cat("\n=== Overall Model Significance ===\n")
print(anova_overall)
cat("\n=== Axis Significance ===\n")
print(anova_axes)
cat("\n=== Variable Significance ===\n")
print(anova_terms)

# Step 11: Calculate variance explained
cat("\nStep 11: Calculating variance components...\n")
total_var <- rda_result$tot.chi
constrained_var <- rda_result$CCA$tot.chi
unconstrained_var <- rda_result$CA$tot.chi
var_explained <- constrained_var / total_var
eigenvalues <- rda_result$CCA$eig
axis_var <- (eigenvalues / total_var) * 100

cat("  - Total variance:", round(total_var, 2), "\n")
cat("  - Proportion explained:", round(var_explained * 100, 2), "%\n")

# Step 12: Create RDA biplot
cat("\nStep 12: Creating RDA biplot...\n")
sample_scores <- scores(rda_result, choices = 1:2, display = "sites", scaling = 1)
rda_df <- data.frame(
  Sample = rownames(sample_scores),
  RDA1 = sample_scores[, 1],
  RDA2 = sample_scores[, 2]
)

env_scores <- scores(rda_result, choices = 1:2, display = "bp", scaling = 1)
env_arrows <- data.frame(
  Variable = rownames(env_scores),
  RDA1 = env_scores[, 1],
  RDA2 = env_scores[, 2]
)

rda_plot <- ggplot(rda_df, aes(x = RDA1, y = RDA2)) +
  geom_point(size = 2.5, alpha = 0.7, color = "steelblue") +
  geom_segment(data = env_arrows, 
               aes(x = 0, y = 0, xend = RDA1 * 3, yend = RDA2 * 3),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "red", size = 1, alpha = 0.8) +
  geom_text(data = env_arrows, 
            aes(x = RDA1 * 3.2, y = RDA2 * 3.2, label = Variable),
            color = "black", size = 3, fontface = "bold") +
  labs(
    title = "RDA Biplot: Genomic Variation Constrained by Environment",
    subtitle = paste0("Total variance explained: ", round(var_explained * 100, 2), "%"),
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  coord_fixed()

# Step 13: Identify candidate SNPs
cat("\nStep 13: Identifying candidate SNPs...\n")
species_scores <- scores(rda_result, choices = 1:2, display = "species", scaling = 1)
loadings_subset <- species_scores[, 1:2]

maha_dist <- mahalanobis(loadings_subset, 
                         center = colMeans(loadings_subset), 
                         cov = cov(loadings_subset))

outlier_threshold <- quantile(maha_dist, OUTLIER_QUANTILE)
outlier_snps <- which(maha_dist > outlier_threshold)

cat("  - Identified", length(outlier_snps), "candidate SNPs\n")

outlier_df <- data.frame(
  SNP_index = outlier_snps,
  RDA1_loading = species_scores[outlier_snps, 1],
  RDA2_loading = species_scores[outlier_snps, 2],
  Mahalanobis_distance = maha_dist[outlier_snps],
  chromosome = vcf@fix[outlier_snps, "CHROM"],
  position = as.numeric(vcf@fix[outlier_snps, "POS"])
)

loading_df <- data.frame(
  SNP = 1:nrow(species_scores),
  RDA1 = species_scores[, 1],
  RDA2 = species_scores[, 2],
  Outlier = 1:nrow(species_scores) %in% outlier_snps
)

loading_plot <- ggplot(loading_df, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = Outlier), size = 1, alpha = 0.6) +
  scale_color_manual(values = c("gray70", "red"), 
                     labels = c("Background SNPs", "Candidate SNPs")) +
  labs(
    title = "SNP Loadings on RDA Axes",
    subtitle = paste("Candidate SNPs (n =", length(outlier_snps), ") highlighted in red"),
    x = paste0("RDA1 (", round(axis_var[1], 2), "%)"),
    y = paste0("RDA2 (", round(axis_var[2], 2), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  coord_fixed()

# Step 14: Save results
cat("\nStep 14: Saving results...\n")

save(rda_result, file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_results.RData")))

ggsave(file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_biplot.png")), 
       rda_plot, width = 10, height = 8, dpi = 300)

ggsave(file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_snp_loadings.png")), 
       loading_plot, width = 10, height = 8, dpi = 300)

write.csv(outlier_df, 
          file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_outlier_snps.csv")), 
          row.names = FALSE)

rda_summary_table <- data.frame(
  Metric = c("Total Variance", "Constrained Variance", "Unconstrained Variance", 
             "Proportion Explained", "Number of Environmental Variables", 
             "Number of Samples", "Number of SNPs", "Candidate SNPs"),
  Value = c(
    round(total_var, 2),
    round(constrained_var, 2), 
    round(unconstrained_var, 2),
    paste0(round(var_explained * 100, 2), "%"),
    ncol(env_vars_scaled),
    nrow(genotype_data),
    ncol(genotype_data),
    length(outlier_snps)
  )
)

write.csv(rda_summary_table, 
          file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_summary_table.csv")), 
          row.names = FALSE)

cat("\n=== RDA Analysis Complete! ===\n")
cat("Results saved to:", OUTPUT_DIR, "\n")
```

---

## Tips and Best Practices

### Data Preparation

1. **Quality control**: Filter SNPs by MAF, HWE, and missing data before RDA
2. **Sample size**: Aim for at least 50 samples for reliable results
3. **Environmental variables**: Use biologically meaningful predictors
4. **Multicollinearity**: Always check and remove correlated variables (r > 0.9)

### Analysis Decisions

1. **Scaling**: Always scale environmental variables (mean = 0, sd = 1)
2. **Missing data**: Impute with mean or remove samples/SNPs with >30% missing
3. **Permutations**: Use at least 100 permutations for significance testing (999 for publication)
4. **Outlier threshold**: 99.9th percentile is conservative; adjust based on your goals

### Validation

1. **Cross-validation**: Use independent datasets to validate candidate SNPs
2. **Population structure**: Run PCA to check if RDA axes correlate with neutral structure
3. **Functional annotation**: Prioritize candidate SNPs in genes or regulatory regions
4. **Alternative methods**: Compare with GEA methods (LFMM, Bayenv2, BayPass)

### Common Issues

**Problem**: Low variance explained (< 2%)
- **Solution**: Check environmental variables for relevance; increase sample size

**Problem**: Non-significant model
- **Solution**: Increase permutations; check for insufficient environmental gradients

**Problem**: Too many/few candidate SNPs
- **Solution**: Adjust `OUTLIER_QUANTILE` parameter; use biological knowledge

**Problem**: Computational slowness
- **Solution**: Reduce SNP number via LD pruning; use subset of SNPs for testing

---

## References

1. Forester, B. R., et al. (2018). Comparing methods for detecting multilocus adaptation with multivariate genotype–environment associations. *Molecular Ecology*, 27(9), 2215-2233.

2. Capblancq, T., & Forester, B. R. (2021). Redundancy analysis: A Swiss Army Knife for landscape genomics. *Methods in Ecology and Evolution*, 12(12), 2298-2309.

3. Legendre, P., & Legendre, L. (2012). *Numerical ecology* (3rd ed.). Elsevier.

4. Rellstab, C., et al. (2015). A practical guide to environmental association analysis in landscape genomics. *Molecular Ecology*, 24(17), 4348-4370.

---

## Contact and Support

For questions or issues with this tutorial:
- Check the vegan package documentation: `?rda`
- Review vcfR tutorials: https://knausb.github.io/vcfR_documentation/
- Post questions on Biostars or Stack Overflow with tags `rda`, `genomics`, `r`

---

**Last updated**: December 2024
