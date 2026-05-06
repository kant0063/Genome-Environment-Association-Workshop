################################################################################
# Environmental Genome-Wide Association Study (eGWAS)
# 
# This script performs environmental association analysis starting from:
# 1. A VCF file containing genotype data
# 2. A CSV file containing environmental data for each sample
#
# The script implements multiple linear mixed models (LMM) with various 
# population structure corrections to identify genetic variants associated
# with environmental variables.
#
################################################################################

# Load required packages ------------------------------------------------------

library(vcfR)           # For reading VCF files
library(sommer)         # For mixed model analysis and FDR correction
library(tidyverse)      # For data manipulation and visualization
library(rrBLUP)         # For GWAS and kinship matrix calculation
library(qvalue)         # For FDR correction
library(QCEWAS)         # For calculating genomic inflation factor
library(qqman)          # For QQ plots and Manhattan plots

# CRITICAL NOTE: Package conflicts
# Several packages mask dplyr functions (especially select, filter, rename)
# To avoid conflicts, these functions are explicitly called with dplyr:: prefix
# throughout this script. Do NOT remove the dplyr:: prefixes or you will get
# "unused arguments" errors.
# 
# Common masking packages:
# - MASS masks dplyr::select()
# - stats masks dplyr::filter()  
# - rrBLUP may mask various functions

# Configuration ---------------------------------------------------------------

# Set your file paths here
vcf_file <- "path/to/your/genotypes.vcf"        # Path to VCF file
env_file <- "path/to/your/environment.csv"      # Path to environmental CSV file
output_dir <- "output"                          # Directory for results
figures_dir <- "figures"                        # Directory for figures

# Create output directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis parameters
min_maf <- 0.05              # Minimum minor allele frequency
fdr_threshold <- 0.20        # FDR threshold for significance
missing_thresh <- 0.10       # Maximum proportion of missing data per marker


# Read and process VCF file ---------------------------------------------------

cat("Reading VCF file...\n")
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# Extract genotype matrix
geno_mat <- extract.gt(vcf, element = "GT")

# Convert VCF genotypes to numeric format (-1, 0, 1, 2)
# This converts genotypes like "0/0", "0/1", "1/1" to numeric values
geno_numeric <- matrix(NA, nrow = nrow(geno_mat), ncol = ncol(geno_mat))
rownames(geno_numeric) <- rownames(geno_mat)
colnames(geno_numeric) <- colnames(geno_mat)

for (i in 1:nrow(geno_mat)) {
  for (j in 1:ncol(geno_mat)) {
    gt <- geno_mat[i, j]
    if (is.na(gt)) {
      geno_numeric[i, j] <- NA
    } else if (gt == "0/0" || gt == "0|0") {
      geno_numeric[i, j] <- -1
    } else if (gt == "0/1" || gt == "0|1" || gt == "1|0") {
      geno_numeric[i, j] <- 0
    } else if (gt == "1/1" || gt == "1|1") {
      geno_numeric[i, j] <- 1
    } else {
      # Handle other cases (e.g., multiallelic)
      alleles <- as.numeric(unlist(strsplit(gt, "[/|]")))
      geno_numeric[i, j] <- sum(alleles) - 1
    }
  }
  if (i %% 1000 == 0) cat("  Processed", i, "of", nrow(geno_mat), "markers\n")
}

# Transpose so samples are rows and markers are columns
geno_mat_t <- t(geno_numeric)

# Remove monomorphic markers (markers with no variation)
cat("  Checking for monomorphic markers...\n")
is_polymorphic <- apply(geno_mat_t, 2, function(x) {
  # Remove missing values
  x_clean <- x[!is.na(x)]
  # Check if there's more than one unique value
  length(unique(x_clean)) > 1
})

n_monomorphic <- sum(!is_polymorphic)
if (n_monomorphic > 0) {
  cat("  Removing", n_monomorphic, "monomorphic markers\n")
  geno_mat_t <- geno_mat_t[, is_polymorphic]
  geno_numeric <- geno_numeric[is_polymorphic, ]
}

# Create SNP information data frame
snp_info <- data.frame(
  marker = rownames(geno_numeric),
  chrom = getCHROM(vcf)[is_polymorphic],
  pos = getPOS(vcf)[is_polymorphic],
  stringsAsFactors = FALSE
)

cat("VCF file loaded:", nrow(snp_info), "markers x", ncol(geno_mat_t), "samples\n")
cat("  (After removing monomorphic markers)\n")


# Read environmental data -----------------------------------------------------

cat("\nReading environmental data...\n")
env_data <- read.csv(env_file, stringsAsFactors = FALSE)

# The CSV file should have the following structure:
# - First column: sample IDs (matching VCF sample names)
# - Subsequent columns: environmental variables
# Example columns: sample_id, latitude, longitude, bio1, bio2, ..., elevation

# Check that sample column exists
if (!"sample_id" %in% colnames(env_data)) {
  # Assume first column is sample ID
  colnames(env_data)[1] <- "sample_id"
  cat("  Note: Assuming first column contains sample IDs\n")
}

# Identify environmental variables (exclude sample_id, latitude, longitude)
env_vars <- setdiff(colnames(env_data), c("sample_id", "latitude", "longitude"))
cat("  Found", length(env_vars), "environmental variables:", 
    paste(head(env_vars, 3), collapse = ", "), "...\n")


# Match samples between genotype and environmental data -----------------------

# Find common samples
common_samples <- intersect(rownames(geno_mat_t), env_data$sample_id)
cat("\nFound", length(common_samples), "samples with both genotype and environmental data\n")

if (length(common_samples) == 0) {
  stop("No common samples found between VCF and environmental data. Check sample ID matching.")
}

# Subset to common samples
geno_mat_t <- geno_mat_t[common_samples, ]
env_data <- env_data[match(common_samples, env_data$sample_id), ]


# Quality control and filtering -----------------------------------------------

cat("\nPerforming quality control...\n")

# Calculate missing data per marker
missing_per_marker <- colMeans(is.na(geno_mat_t))
markers_to_keep <- missing_per_marker <= missing_thresh

cat("  Removed", sum(!markers_to_keep), "markers with >", 
    missing_thresh * 100, "% missing data\n")

geno_mat_filtered <- geno_mat_t[, markers_to_keep]
snp_info_filtered <- snp_info[markers_to_keep, ]

# Calculate minor allele frequency
maf <- apply(geno_mat_filtered, 2, function(x) {
  x <- x[!is.na(x)]
  p <- (sum(x == 1) + 2 * sum(x == 0)) / (2 * length(x))
  return(min(p, 1 - p))
})

# Filter by MAF
markers_pass_maf <- maf >= min_maf
cat("  Removed", sum(!markers_pass_maf), "markers with MAF <", min_maf, "\n")

geno_mat_final <- geno_mat_filtered[, markers_pass_maf]
snp_info_final <- snp_info_filtered[markers_pass_maf, ]

cat("  Final dataset:", nrow(geno_mat_final), "samples x", 
    ncol(geno_mat_final), "markers\n")

# Impute missing genotypes with column means
for (j in 1:ncol(geno_mat_final)) {
  missing_idx <- is.na(geno_mat_final[, j])
  if (any(missing_idx)) {
    geno_mat_final[missing_idx, j] <- mean(geno_mat_final[, j], na.rm = TRUE)
  }
}


# Calculate kinship matrix ----------------------------------------------------

cat("\nCalculating kinship matrix...\n")
K <- A.mat(geno_mat_final, min.MAF = 0, max.missing = 1)
cat("  Kinship matrix calculated:", nrow(K), "x", ncol(K), "\n")


# Calculate principal components from kinship matrix -------------------------

cat("Calculating principal components...\n")
K_pca <- prcomp(K)
K_eigens <- K_pca$x[, 1:3] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  dplyr::rename(PC1 = PC1, PC2 = PC2, PC3 = PC3)

cat("  PC1 variance explained:", 
    round(summary(K_pca)$importance[2, 1] * 100, 2), "%\n")


# Prepare data for GWAS -------------------------------------------------------

cat("\nPreparing data for GWAS...\n")

# Format genotype data for rrBLUP GWAS function
geno_df <- t(geno_mat_final) %>%
  as.data.frame() %>%
  rownames_to_column("marker") %>%
  dplyr::inner_join(snp_info_final, ., by = "marker") %>%
  as.data.frame()

# Ensure the column order is correct: marker, chrom, pos, then individuals
# and that column names match the individuals in pheno and K
col_order <- c("marker", "chrom", "pos", rownames(geno_mat_final))
geno_df <- geno_df[, col_order]

cat("  Genotype data formatted for GWAS:\n")
cat("    Markers:", nrow(geno_df), "\n")
cat("    Individuals:", ncol(geno_df) - 3, "\n")
cat("    First few individual IDs:", paste(head(colnames(geno_df)[-(1:3)], 3), collapse = ", "), "\n")

# Prepare phenotype (environmental) data
pheno <- env_data %>%
  dplyr::select(sample_id, all_of(env_vars)) %>%
  dplyr::rename(individual = sample_id) %>%
  as.data.frame()

# Handle missing environmental data
for (col in env_vars) {
  if (any(is.na(pheno[[col]]))) {
    cat("  Warning: Missing values in", col, "- imputing with median\n")
    pheno[[col]][is.na(pheno[[col]])] <- median(pheno[[col]], na.rm = TRUE)
  }
}


# Define models for testing ---------------------------------------------------

# Test multiple models with different population structure corrections:
# 1. K only (kinship matrix)
# 2. K + 3 PCs (kinship + principal components)
# 3. K + latitude
# 4. K + 3 PCs + latitude
# 5. K + latitude + longitude
# 6. K + 3 PCs + latitude + longitude

egwas_models <- tribble(
  ~model, ~use_K, ~n_PCs, ~covariates,
  "model1", TRUE, 0, character(),
  "model2", TRUE, 3, character(),
  "model3", TRUE, 0, c("latitude"),
  "model4", TRUE, 3, c("latitude"),
  "model5", TRUE, 0, c("latitude", "longitude"),
  "model6", TRUE, 3, c("latitude", "longitude")
) %>%
  mutate(
    n_params = pmap_dbl(list(n_PCs, covariates), 
                        ~1 + 1 + ..1 + length(..2) + 1),
    results = list(NULL)
  )


# Run environmental GWAS ------------------------------------------------------

cat("\nRunning environmental GWAS...\n")
cat("  Testing", nrow(egwas_models), "different models\n")
cat("  Testing", length(env_vars), "environmental variables\n")

# Add geographic coordinates to phenotype data if they exist
if (all(c("latitude", "longitude") %in% colnames(env_data))) {
  pheno <- pheno %>%
    dplyr::left_join(
      env_data %>% dplyr::select(sample_id, latitude, longitude) %>%
        dplyr::rename(individual = sample_id),
      by = "individual"
    )
  cat("  Added latitude and longitude to phenotype data\n")
} else {
  cat("  Note: latitude/longitude not found in environmental data\n")
}

cat("  Phenotype data columns:", paste(colnames(pheno), collapse = ", "), "\n")

# Iterate over models
for (i in seq_len(nrow(egwas_models))) {
  
  model_name <- egwas_models$model[i]
  cat("\n  Running", model_name, "...\n")
  
  # Extract model parameters
  n_pcs <- egwas_models$n_PCs[i]
  covariates <- egwas_models$covariates[[i]]
  
  if (length(covariates) > 0) {
    cat("    Model covariates:", paste(covariates, collapse = ", "), "\n")
  }
  
  # Storage for results
  results_list <- list()
  
  # Iterate over environmental variables
  for (var_idx in seq_along(env_vars)) {
    env_var <- env_vars[var_idx]
    
    if (var_idx %% 5 == 0) {
      cat("    Progress:", var_idx, "/", length(env_vars), "variables\n")
    }
    
    # Remove the current variable from covariates if present
    covariates_use <- setdiff(covariates, env_var)
    if (length(covariates_use) == 0) covariates_use <- NULL
    
    # Prepare phenotype for this variable
    # Only include covariates that actually exist in pheno
    cols_to_use <- c("individual", env_var)
    if (!is.null(covariates_use)) {
      # Filter covariates to only those that exist in pheno
      covariates_use <- covariates_use[covariates_use %in% colnames(pheno)]
      if (length(covariates_use) > 0) {
        cols_to_use <- c(cols_to_use, covariates_use)
      } else {
        covariates_use <- NULL
      }
    }
    
    pheno_var <- pheno[, cols_to_use, drop = FALSE]
    
    # Skip if any required columns are missing
    if (!env_var %in% colnames(pheno)) {
      cat("      Warning: Skipping", env_var, "- variable not found in pheno\n")
      next
    }
    
    # If we have covariates, adjust the phenotype values
    if (!is.null(covariates_use) && length(covariates_use) > 0) {
      # Fit a linear model to adjust for covariates
      covariate_formula <- reformulate(covariates_use, response = env_var)
      covariate_model <- lm(covariate_formula, data = pheno_var)
      # Get residuals (phenotype adjusted for covariates)
      pheno_var[[env_var]] <- residuals(covariate_model)
    }
    
    # Prepare final phenotype data (just individual and the environmental variable)
    pheno_var_final <- pheno_var[, c("individual", env_var), drop = FALSE]
    
    # CRITICAL: Ensure pheno is a data frame with exactly 2 columns
    # and that the individual IDs match between pheno, geno, and K
    pheno_var_final <- as.data.frame(pheno_var_final)
    
    # Check that individuals match between pheno and geno
    geno_individuals <- colnames(geno_df)[-(1:3)]  # Skip marker, chrom, pos
    pheno_individuals <- pheno_var_final$individual
    
    if (!all(pheno_individuals %in% geno_individuals)) {
      cat("      Warning: Not all pheno individuals found in geno for", env_var, "\n")
      next
    }
    
    # Check that individuals match between pheno and K
    K_individuals <- rownames(K)
    if (!all(pheno_individuals %in% K_individuals)) {
      cat("      Warning: Not all pheno individuals found in K for", env_var, "\n")
      next
    }
    
    # Run GWAS using rrBLUP
    tryCatch({
      gwas_result <- rrBLUP::GWAS(
        pheno = pheno_var_final,
        geno = geno_df,
        n.PC = n_pcs,
        min.MAF = 0,  # Already filtered
        P3D = TRUE,   # Population parameters previously determined
        plot = FALSE,
        K = K
      )
      
      results_list[[env_var]] <- gwas_result
      
    }, error = function(e) {
      cat("      Error in", env_var, ":", e$message, "\n")
      cat("      Pheno dimensions:", nrow(pheno_var_final), "x", ncol(pheno_var_final), "\n")
      cat("      Pheno column names:", paste(colnames(pheno_var_final), collapse = ", "), "\n")
      cat("      First few pheno individuals:", paste(head(pheno_var_final$individual, 3), collapse = ", "), "\n")
      cat("      First few geno individuals:", paste(head(geno_individuals, 3), collapse = ", "), "\n")
      cat("      First few K individuals:", paste(head(K_individuals, 3), collapse = ", "), "\n")
    })
  }
  
  # Combine results for this model
  if (length(results_list) > 0) {
    results_df <- reduce(results_list, full_join, 
                        by = c("marker", "chrom", "pos"))
    egwas_models$results[[i]] <- results_df
  }
}


# Combine and analyze results -------------------------------------------------

cat("\nCombining results across models...\n")

# Combine all model results
gwas_combined <- egwas_models %>%
  dplyr::select(model, results) %>%
  dplyr::filter(!map_lgl(results, is.null)) %>%
  unnest(results) %>%
  pivot_longer(
    cols = all_of(env_vars),
    names_to = "variable",
    values_to = "score"
  ) %>%
  mutate(
    p_value = 10^-score,
    score = ifelse(is.na(score), 0, score)
  )

# Store as gwas_analyzed (no need to calculate uniform scores - qqman will do this)
gwas_analyzed <- gwas_combined

# Calculate genomic inflation factor (lambda)
p_lambda <- gwas_analyzed %>%
  dplyr::filter(!is.na(p_value), p_value > 0) %>%
  group_by(variable, model) %>%
  summarize(
    lambda = QCEWAS::P_lambda(p = p_value),
    .groups = "drop"
  ) %>%
  group_by(variable) %>%
  mutate(
    best_model = model[which.min(abs(1 - lambda))],
    is_best = model == best_model
  ) %>%
  ungroup()

cat("  Genomic inflation factors calculated\n")
cat("  Best models (closest to lambda = 1):\n")
print(p_lambda %>% dplyr::filter(is_best) %>% dplyr::select(variable, model, lambda))


# Create QQ plots -------------------------------------------------------------

# Create QQ plots -------------------------------------------------------------

cat("\nGenerating QQ plots...\n")

# Create a multi-panel QQ plot using qqman
# We'll create individual QQ plots for each variable and model combination

# Set up the plotting layout
n_vars <- length(env_vars)
n_cols <- min(3, n_vars)
n_rows <- ceiling(n_vars / n_cols)

# Create PDF with QQ plots
pdf(file.path(figures_dir, "egwas_qqplots.pdf"), width = 14, height = 4 * n_rows)

par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

# Create QQ plots for each variable
for (var in env_vars) {
  cat("  Creating QQ plot for", var, "...\n")
  
  # Get data for this variable and the best model
  best_model_var <- p_lambda %>%
    dplyr::filter(variable == var, is_best) %>%
    pull(model)
  
  if (length(best_model_var) == 0) {
    best_model_var <- "model2"  # Default to model2 if no best model found
  }
  
  var_data <- gwas_analyzed %>%
    dplyr::filter(variable == var, model == best_model_var) %>%
    dplyr::filter(!is.na(p_value), p_value > 0)
  
  if (nrow(var_data) > 0) {
    # Get lambda value
    lambda_val <- p_lambda %>%
      dplyr::filter(variable == var, model == best_model_var) %>%
      pull(lambda)
    
    if (length(lambda_val) == 0) lambda_val <- NA
    
    # Create QQ plot using qqman
    qq(var_data$p_value, 
       main = paste0(var, " (", best_model_var, ", λ=", round(lambda_val, 3), ")"),
       col = "blue4",
       cex = 0.5)
  } else {
    plot.new()
    text(0.5, 0.5, paste("No data for", var), cex = 1.5)
  }
}

dev.off()

cat("  QQ plots saved to:", file.path(figures_dir, "egwas_qqplots.pdf"), "\n")

# Also create a comparison plot showing all models for a few key variables
if (length(env_vars) > 0) {
  
  cat("  Creating model comparison QQ plots...\n")
  
  # Select up to 4 variables to show model comparisons
  vars_to_plot <- head(env_vars, 4)
  
  pdf(file.path(figures_dir, "egwas_qqplots_model_comparison.pdf"), 
      width = 16, height = 4 * length(vars_to_plot))
  
  par(mfrow = c(length(vars_to_plot), 6), mar = c(4, 4, 3, 1))
  
  for (var in vars_to_plot) {
    for (mod in unique(egwas_models$model)) {
      
      var_data <- gwas_analyzed %>%
        dplyr::filter(variable == var, model == mod) %>%
        dplyr::filter(!is.na(p_value), p_value > 0)
      
      if (nrow(var_data) > 0) {
        lambda_val <- p_lambda %>%
          dplyr::filter(variable == var, model == mod) %>%
          pull(lambda)
        
        is_best_model <- p_lambda %>%
          dplyr::filter(variable == var, model == mod) %>%
          pull(is_best)
        
        if (length(lambda_val) == 0) lambda_val <- NA
        if (length(is_best_model) == 0) is_best_model <- FALSE
        
        col_use <- if (is_best_model) "red3" else "blue4"
        
        qq(var_data$p_value,
           main = paste0(var, " - ", mod, "\nλ=", round(lambda_val, 3)),
           col = col_use,
           cex = 0.4)
      } else {
        plot.new()
        text(0.5, 0.5, "No data", cex = 1)
      }
    }
  }
  
  dev.off()
  
  cat("  Model comparison QQ plots saved to:", 
      file.path(figures_dir, "egwas_qqplots_model_comparison.pdf"), "\n")
}


# Create lambda boxplots ---------------------------------------------------

cat("\nGenerating lambda value boxplots...\n")

# Create boxplot showing lambda values for each variable across models
lambda_boxplot <- p_lambda %>%
  ggplot(aes(x = variable, y = lambda, fill = model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = model), width = 0.2, size = 2, alpha = 0.8) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = c(0.95, 1.05), linetype = "dotted", color = "orange", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Genomic Inflation Factor (λ) Across Models",
    subtitle = "Red dashed line = ideal (λ=1.0); Orange dotted lines = acceptable range (0.95-1.05)",
    x = "Environmental Variable",
    y = "Lambda (λ)",
    fill = "Model",
    color = "Model"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(figures_dir, "lambda_boxplot.png"),
  plot = lambda_boxplot,
  width = 12,
  height = 6,
  dpi = 300
)

cat("  Lambda boxplot saved to:", file.path(figures_dir, "lambda_boxplot.png"), "\n")

# Also create a heatmap showing lambda values
lambda_heatmap <- p_lambda %>%
  ggplot(aes(x = model, y = variable, fill = lambda)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(lambda, 2)), size = 3) +
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 1.0,
    limits = c(0.8, max(1.2, max(p_lambda$lambda, na.rm = TRUE))),
    name = "Lambda (λ)"
  ) +
  labs(
    title = "Genomic Inflation Factor Heatmap",
    subtitle = "White = ideal (λ=1.0); Blue = under-inflation; Red = over-inflation",
    x = "Model",
    y = "Environmental Variable"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(figures_dir, "lambda_heatmap.png"),
  plot = lambda_heatmap,
  width = 10,
  height = max(6, length(env_vars) * 0.4),
  dpi = 300
)

cat("  Lambda heatmap saved to:", file.path(figures_dir, "lambda_heatmap.png"), "\n")

# Create a summary plot showing which model is best for each variable
best_model_summary <- p_lambda %>%
  dplyr::filter(is_best) %>%
  dplyr::count(model) %>%
  ggplot(aes(x = model, y = n, fill = model)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Best Model Selection Summary",
    subtitle = paste0("Best model = closest to λ=1.0 (n=", length(env_vars), " variables)"),
    x = "Model",
    y = "Number of Variables",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = file.path(figures_dir, "best_model_summary.png"),
  plot = best_model_summary,
  width = 10,
  height = 6,
  dpi = 300
)

cat("  Best model summary saved to:", file.path(figures_dir, "best_model_summary.png"), "\n")


# Identify significant markers ------------------------------------------------

cat("\nIdentifying significant markers...\n")

# Use model2 (K + 3 PCs) as the default best model
# You can modify this based on your p_lambda results
sig_markers <- gwas_analyzed %>%
  dplyr::filter(model == "model2") %>%
  group_by(variable) %>%
  group_split() %>%
  map_dfr(function(df) {
    var_name <- df$variable[1]
    p_vals <- df$p_value[!is.na(df$p_value)]
    
    if (length(p_vals) == 0) return(tibble())
    
    # Calculate FDR
    # Try Storey's q-value first (more powerful for eGWAS); fall back to
    # Benjamini-Hochberg if pi0 estimation fails (too few markers, weird
    # p-value distribution, etc.)
    qobj <- tryCatch(
      qvalue::qvalue(p_vals),
      error = function(e) NULL,
      warning = function(w) NULL
    )
    
    df_with_q <- df %>%
      dplyr::filter(!is.na(p_value))
    
    if (!is.null(qobj)) {
      df_with_q$q_value <- qobj$qvalues
      method_used <- "Storey q-value"
    } else {
      df_with_q$q_value <- p.adjust(df_with_q$p_value, method = "BH")
      method_used <- "BH (qvalue fallback)"
    }
    
    sig <- df_with_q %>%
      dplyr::filter(q_value <= fdr_threshold) %>%
      arrange(p_value)
    
    if (nrow(sig) > 0) {
      cat("  ", var_name, ":", nrow(sig), "significant markers [", method_used, "]\n")
    }
    
    return(sig)
  })

cat("\nTotal significant markers:", nrow(sig_markers), "\n")


# Create Manhattan plots ------------------------------------------------------

if (nrow(sig_markers) > 0) {
  cat("\nGenerating Manhattan plots...\n")
  
  # Create Manhattan plots using qqman for better visualization
  pdf(file.path(figures_dir, "egwas_manhattan.pdf"), width = 14, height = 4 * n_rows)
  
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))
  
  for (var in env_vars) {
    cat("  Creating Manhattan plot for", var, "...\n")
    
    # Get data for this variable and the best model
    best_model_var <- p_lambda %>%
      dplyr::filter(variable == var, is_best) %>%
      pull(model)
    
    if (length(best_model_var) == 0) {
      best_model_var <- "model2"
    }
    
    var_data <- gwas_analyzed %>%
      dplyr::filter(variable == var, model == best_model_var) %>%
      dplyr::filter(!is.na(p_value), p_value > 0)
    
    if (nrow(var_data) > 0) {
      # Prepare data for manhattan plot
      # qqman expects: SNP, CHR, BP, P
      manhattan_data <- var_data %>%
        dplyr::select(marker, chrom, pos, p_value) %>%
        dplyr::rename(SNP = marker, CHR = chrom, BP = pos, P = p_value)
      
      # Convert chromosome to numeric if needed
      if (!is.numeric(manhattan_data$CHR)) {
        manhattan_data$CHR <- as.numeric(factor(manhattan_data$CHR))
      }
      
      # Create Manhattan plot
      tryCatch({
        manhattan(manhattan_data,
                 chr = "CHR",
                 bp = "BP",
                 p = "P",
                 snp = "SNP",
                 main = paste0(var, " (", best_model_var, ")"),
                 col = c("blue4", "orange3"),
                 suggestiveline = -log10(0.05),
                 genomewideline = -log10(fdr_threshold / nrow(manhattan_data)),
                 cex = 0.5)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error plotting", var), cex = 1)
      })
    } else {
      plot.new()
      text(0.5, 0.5, paste("No data for", var), cex = 1.5)
    }
  }
  
  dev.off()
  
  cat("  Manhattan plots saved to:", file.path(figures_dir, "egwas_manhattan.pdf"), "\n")
} else {
  cat("\nNo significant markers found - skipping Manhattan plots\n")
}


# Calculate variance explained by significant markers ------------------------

cat("\nCalculating variance explained...\n")

variance_explained <- list()

if (nrow(sig_markers) > 0) {
  
  # Merge phenotype with PCs
  pheno_with_pcs <- pheno %>%
    dplyr::left_join(K_eigens, by = c("individual" = "sample_id"))
  
  for (var_name in unique(sig_markers$variable)) {
    
    cat("  Processing", var_name, "...\n")
    
    # Get significant markers for this variable
    sig_mar_var <- sig_markers %>%
      dplyr::filter(variable == var_name) %>%
      pull(marker)
    
    if (length(sig_mar_var) == 0) next
    
    # Prepare model components
    y <- pheno_with_pcs[[var_name]]
    X <- model.matrix(~ PC1 + PC2 + PC3, data = pheno_with_pcs)
    Z <- model.matrix(~ -1 + individual, data = pheno_with_pcs)
    
    # Fit mixed model to get residuals
    mm_fit <- mixed.solve(y = y, Z = Z, K = K, X = X, method = "REML")
    
    # Calculate adjusted values (residuals after accounting for structure)
    y_adj <- y - ((X %*% mm_fit$beta) + (Z %*% mm_fit$u))
    
    # Get genotypes for significant markers
    geno_sig <- geno_mat_final[, sig_mar_var, drop = FALSE]
    
    # Fit model with significant markers
    model_data <- data.frame(y_adj = as.numeric(y_adj), geno_sig)
    
    if (ncol(geno_sig) > 1) {
      marker_model <- lm(y_adj ~ ., data = model_data)
      # Stepwise selection
      marker_model_step <- step(marker_model, trace = 0)
    } else {
      marker_model_step <- lm(y_adj ~ ., data = model_data)
    }
    
    # Calculate variance components
    pop_structure_r2 <- var(X %*% mm_fit$beta) / var(y)
    polygenic_r2 <- mm_fit$Vu / var(y)
    residual_r2 <- 1 - (pop_structure_r2 + polygenic_r2)
    marker_r2 <- summary(marker_model_step)$adj.r.squared
    
    variance_explained[[var_name]] <- tibble(
      variable = var_name,
      n_sig_markers = length(sig_mar_var),
      n_retained_markers = length(coef(marker_model_step)) - 1,
      pop_structure_var_exp = pop_structure_r2,
      polygenic_var_exp = polygenic_r2,
      markers_var_exp = marker_r2 * residual_r2,
      total_var_exp = pop_structure_r2 + polygenic_r2 + (marker_r2 * residual_r2)
    )
  }
  
  variance_df <- bind_rows(variance_explained)
  
  cat("\nVariance explained summary:\n")
  print(variance_df)
  
} else {
  cat("  No significant markers found - skipping variance calculation\n")
  variance_df <- tibble()
}


# Save results ----------------------------------------------------------------

cat("\nSaving results...\n")

# Save R objects (including genotype and environmental data for visualizations)
save(
  gwas_combined,
  gwas_analyzed,
  p_lambda,
  sig_markers,
  variance_df,
  egwas_models,
  K,
  snp_info_final,
  geno_mat_final,  # Needed for allele-environment boxplots
  env_data,        # Needed for allele-environment boxplots and maps
  file = file.path(output_dir, "egwas_results.RData")
)

# Save significant markers to CSV
if (nrow(sig_markers) > 0) {
  # Remove any list columns before saving to CSV
  # Explicitly select the columns we want to save
  sig_markers_clean <- sig_markers %>%
    dplyr::select(any_of(c("marker", "chrom", "pos", "model", "variable", 
                           "score", "p_value", "unif_score"))) %>%
    as.data.frame()
  
  write.csv(
    sig_markers_clean,
    file = file.path(output_dir, "significant_markers.csv"),
    row.names = FALSE
  )
  cat("  Significant markers saved to:", 
      file.path(output_dir, "significant_markers.csv"), "\n")
}

# Save variance explained to CSV
if (nrow(variance_df) > 0) {
  write.csv(
    variance_df,
    file = file.path(output_dir, "variance_explained.csv"),
    row.names = FALSE
  )
  cat("  Variance components saved to:", 
      file.path(output_dir, "variance_explained.csv"), "\n")
}

# Save lambda values
p_lambda_clean <- p_lambda %>%
  dplyr::select(any_of(c("variable", "model", "lambda", "best_model", "is_best"))) %>%
  as.data.frame()

write.csv(
  p_lambda_clean,
  file = file.path(output_dir, "genomic_inflation_factors.csv"),
  row.names = FALSE
)

cat("\nAnalysis complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("Figures saved to:", figures_dir, "\n")
