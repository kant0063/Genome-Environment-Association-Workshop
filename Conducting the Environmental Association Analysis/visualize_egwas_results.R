################################################################################
# Environmental GWAS Results Visualization
# 
# This script creates comprehensive visualizations of environmental GWAS results
# using the significant markers identified from the environmental_gwas_from_vcf.R
# script.
#
# Required input files:
# - output/egwas_results.RData (from environmental_gwas_from_vcf.R)
# - output/significant_markers.csv (from environmental_gwas_from_vcf.R)
#
################################################################################

# Load required packages ------------------------------------------------------

library(tidyverse)
library(patchwork)
library(ggrepel)
library(scales)
library(viridis)

# Note: MASS (loaded by tidyverse dependencies) can mask dplyr::select()
# All select() calls use explicit dplyr:: prefix to avoid conflicts

# Set plot theme
theme_set(theme_bw())

# Configuration ---------------------------------------------------------------

# Paths to input files
results_file <- "output/egwas_results.RData"
sig_markers_file <- "output/significant_markers.csv"

# Output directory for figures
fig_dir <- "figures/visualizations"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
fdr_threshold <- 0.20  # Should match the threshold used in main script


# Load data -------------------------------------------------------------------

cat("Loading environmental GWAS results...\n")

# Load the R objects
load(results_file)

# Load significant markers
sig_markers <- read.csv(sig_markers_file, stringsAsFactors = FALSE)

cat("  Loaded", nrow(sig_markers), "significant marker-variable associations\n")
cat("  Unique markers:", n_distinct(sig_markers$marker), "\n")
cat("  Unique variables:", n_distinct(sig_markers$variable), "\n")


# Summary statistics ----------------------------------------------------------

cat("\nCalculating summary statistics...\n")

# Number of significant associations per variable
sig_per_variable <- sig_markers %>%
  group_by(variable) %>%
  summarize(
    n_markers = n_distinct(marker),
    n_associations = n(),
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE)
  ) %>%
  arrange(desc(n_markers))

cat("  Variables with most associations:\n")
print(head(sig_per_variable, 5))

# Number of associations per marker
sig_per_marker <- sig_markers %>%
  group_by(marker, chrom, pos) %>%
  summarize(
    n_variables = n_distinct(variable),
    variables = paste(variable, collapse = ", "),
    mean_score = mean(score, na.rm = TRUE),
    max_score = max(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_variables))

cat("\n  Markers associated with most variables:\n")
print(head(sig_per_marker, 5))


# Visualization 1: Summary of significant associations -----------------------

cat("\nCreating summary visualizations...\n")

# Bar plot: Number of significant markers per variable
p1 <- ggplot(sig_per_variable, aes(x = reorder(variable, n_markers), y = n_markers)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = n_markers), hjust = -0.3, size = 3) +
  coord_flip() +
  labs(
    title = "Number of Significant Markers per Environmental Variable",
    subtitle = paste0("FDR threshold = ", fdr_threshold),
    x = "Environmental Variable",
    y = "Number of Significant Markers"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "01_markers_per_variable.png"),
  plot = p1,
  width = 10,
  height = max(6, n_distinct(sig_markers$variable) * 0.3),
  dpi = 300
)

# Bar plot: Number of variables per marker (top 20)
p2 <- sig_per_marker %>%
  head(20) %>%
  ggplot(aes(x = reorder(marker, n_variables), y = n_variables)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = n_variables), hjust = -0.3, size = 3) +
  coord_flip() +
  labs(
    title = "Top 20 Markers Associated with Multiple Variables",
    subtitle = "Markers showing pleiotropy across environmental traits",
    x = "Marker",
    y = "Number of Associated Variables"
  ) +
  theme(
    axis.text.y = element_text(size = 8, family = "mono"),
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "02_pleiotropic_markers.png"),
  plot = p2,
  width = 10,
  height = 8,
  dpi = 300
)

cat("  Summary plots saved\n")


# Visualization 2: Distribution of association scores ------------------------

cat("Creating score distribution plots...\n")

# Histogram of -log10(p) scores
p3 <- ggplot(sig_markers, aes(x = score)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = median(sig_markers$score, na.rm = TRUE), 
             linetype = "dashed", color = "red", size = 1) +
  labs(
    title = "Distribution of Association Scores for Significant Markers",
    subtitle = "Red line = median score",
    x = "-log10(p-value)",
    y = "Count"
  )

ggsave(
  filename = file.path(fig_dir, "03_score_distribution.png"),
  plot = p3,
  width = 8,
  height = 6,
  dpi = 300
)

# Boxplot of scores by variable
p4 <- sig_markers %>%
  left_join(sig_per_variable %>% dplyr::select(variable, n_markers), by = "variable") %>%
  ggplot(aes(x = reorder(variable, -n_markers), y = score)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  coord_flip() +
  labs(
    title = "Distribution of Association Scores by Variable",
    x = "Environmental Variable",
    y = "-log10(p-value)"
  ) +
  theme(axis.text.y = element_text(size = 9))

ggsave(
  filename = file.path(fig_dir, "04_scores_by_variable.png"),
  plot = p4,
  width = 10,
  height = max(6, n_distinct(sig_markers$variable) * 0.3),
  dpi = 300
)

cat("  Score distribution plots saved\n")


# Visualization 3: Chromosome distribution -----------------------------------

cat("Creating chromosome distribution plots...\n")

# Count markers per chromosome
markers_per_chrom <- sig_markers %>%
  distinct(marker, chrom) %>%
  count(chrom, name = "n_markers") %>%
  arrange(chrom)

p5 <- ggplot(markers_per_chrom, aes(x = factor(chrom), y = n_markers)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_text(aes(label = n_markers), vjust = -0.5, size = 3) +
  labs(
    title = "Distribution of Significant Markers Across Chromosomes",
    x = "Chromosome",
    y = "Number of Unique Markers"
  ) +
  theme(panel.grid.major.x = element_blank())

ggsave(
  filename = file.path(fig_dir, "05_markers_per_chromosome.png"),
  plot = p5,
  width = 10,
  height = 6,
  dpi = 300
)

# Chromosome positions of significant markers
# First, categorize variables by type
sig_markers_categorized <- sig_markers %>%
  distinct(marker, chrom, pos, variable) %>%
  mutate(
    env_type = case_when(
      # Temperature variables: bio1 through bio11
      variable %in% c("bio_01", "bio_02", "bio_03", "bio_04", "bio_05", "bio_06", 
                      "bio_07", "bio_08", "bio_09", "bio_10", "bio_11") ~ "Temperature",
      # Precipitation variables: bio12 through bio19
      variable %in% c("bio_12", "bio_13", "bio_14", "bio_15", "bio_16", 
                      "bio_17", "bio_18", "bio_19") ~ "Precipitation",
      # Everything else is soil/other
      TRUE ~ "Soil/Other"
    )
  )

# Print summary to verify categorization
cat("\nVariable categorization for p6:\n")
print(table(sig_markers_categorized$env_type))
cat("\nBreakdown by variable:\n")
print(sig_markers_categorized %>% 
        count(variable, env_type) %>% 
        arrange(env_type, variable))

p6 <- ggplot(sig_markers_categorized, 
             aes(x = pos / 1e6, y = factor(chrom), 
                 color = env_type, shape = variable)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(
    name = "Environmental Type",
    values = c(
      "Temperature" = "#d62728",      # Red
      "Precipitation" = "#1f77b4",    # Blue
      "Soil/Other" = "#2ca02c"        # Green
    )
  ) +
  scale_shape_manual(
    name = "Variable",
    values = rep(c(16, 17, 15, 18, 8, 9, 10, 11, 12, 13, 14, 7, 6, 5, 4, 3, 2, 1, 0), 
                 length.out = n_distinct(sig_markers_categorized$variable))
  ) +
  labs(
    title = "Genomic Positions of Significant Markers",
    subtitle = "Colored by environmental type, shaped by variable",
    x = "Position (Mb)",
    y = "Chromosome"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "gray90"),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 4)),
    shape = guide_legend(order = 2, ncol = 1)
  )

ggsave(
  filename = file.path(fig_dir, "06_marker_positions.png"),
  plot = p6,
  width = 14,
  height = 8,
  dpi = 300
)

cat("  Marker positions plot saved (with color coding by environmental type)\n")

cat("  Chromosome distribution plots saved\n")


# Visualization 4: Heatmap of associations -----------------------------------

cat("Creating association heatmap...\n")

# Create a presence/absence matrix for marker-variable associations
assoc_matrix <- sig_markers %>%
  distinct(marker, variable) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = variable, values_from = present, values_fill = 0)

# For visualization, select top markers by number of associations
top_markers <- sig_per_marker %>%
  head(30) %>%
  pull(marker)

# Create heatmap data
heatmap_data <- sig_markers %>%
  dplyr::filter(marker %in% top_markers) %>%
  dplyr::select(marker, variable, score)

p7 <- ggplot(heatmap_data, aes(x = variable, y = marker, fill = score)) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_viridis(name = "-log10(p)", option = "plasma") +
  labs(
    title = "Association Strength Heatmap",
    subtitle = "Top 30 markers by number of associations",
    x = "Environmental Variable",
    y = "Marker"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7, family = "mono"),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(fig_dir, "07_association_heatmap.png"),
  plot = p7,
  width = 12,
  height = 10,
  dpi = 300
)

cat("  Association heatmap saved\n")


# Visualization 5: Network-style plot of marker-variable relationships -------

cat("Creating network visualization...\n")

# For visualization, use top markers and variables
top_variables <- sig_per_variable %>%
  head(10) %>%
  pull(variable)

network_data <- sig_markers %>%
  dplyr::filter(marker %in% top_markers[1:15], variable %in% top_variables) %>%
  dplyr::select(marker, variable, score)

# Create a simplified circular plot
p8 <- ggplot(network_data, aes(x = variable, y = score, group = marker, color = marker)) +
  geom_line(alpha = 0.4, size = 1) +
  geom_point(size = 3, alpha = 0.8) +
  coord_polar() +
  labs(
    title = "Marker-Variable Association Network",
    subtitle = "Top 15 pleiotropic markers × Top 10 variables",
    y = "-log10(p-value)",
    color = "Marker"
  ) +
  theme(
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    legend.text = element_text(size = 7, family = "mono")
  )

ggsave(
  filename = file.path(fig_dir, "08_network_plot.png"),
  plot = p8,
  width = 12,
  height = 10,
  dpi = 300
)

cat("  Network plot saved\n")


# Visualization 9: Allele-environment boxplots ----------------------------------

cat("Creating allele-environment association boxplots...\n")

# Wrap in tryCatch to capture any errors
tryCatch({
  
  # Check what objects are available after loading the RData file
  cat("  Checking available data objects...\n")
  available_objects <- ls()
  cat("    Total objects loaded:", length(available_objects), "\n")
  cat("    Objects: ", paste(head(available_objects, 10), collapse = ", "), 
      if(length(available_objects) > 10) "..." else "", "\n", sep = "")
  
  # The main GWAS script saves these objects
  # We need: genotype matrix and environmental data with coordinates
  # Check for different possible names
  
  has_geno <- FALSE
  has_env <- FALSE
  
  # Check for genotype matrix (could have different names)
  if (exists("geno_mat_final")) {
    geno_matrix <- geno_mat_final
    has_geno <- TRUE
    cat("  ✓ Found genotype matrix: geno_mat_final\n")
    cat("    Dimensions:", nrow(geno_mat_final), "samples ×", ncol(geno_mat_final), "markers\n")
  } else if (exists("geno_mat")) {
    geno_matrix <- geno_mat
    has_geno <- TRUE
    cat("  ✓ Found genotype matrix: geno_mat\n")
    cat("    Dimensions:", nrow(geno_mat), "samples ×", ncol(geno_mat), "markers\n")
  } else {
    cat("  ✗ No genotype matrix found\n")
    cat("    Looked for: geno_mat_final, geno_mat\n")
    cat("    These are needed for figures 09 and 11\n")
  }
  
  # Check for environmental data (could have different names)
  if (exists("env_data")) {
    env_dataset <- env_data
    has_env <- TRUE
    cat("  ✓ Found environmental data: env_data\n")
    cat("    Dimensions:", nrow(env_data), "samples ×", ncol(env_data), "variables\n")
  } else if (exists("pheno")) {
    env_dataset <- pheno
    has_env <- TRUE
    cat("  ✓ Found environmental data: pheno\n")
    cat("    Dimensions:", nrow(pheno), "samples ×", ncol(pheno), "variables\n")
  } else {
    cat("  ✗ No environmental data found\n")
    cat("    Looked for: env_data, pheno\n")
    cat("    These are needed for figures 09 and 11\n")
  }
  
  # Summary message
  if (!has_geno || !has_env) {
    cat("\n")
    cat("  ═══════════════════════════════════════════════════════════\n")
    cat("  NOTE: Advanced visualizations (figures 09, 11) not available\n")
    cat("  ═══════════════════════════════════════════════════════════\n")
    cat("  To enable these visualizations:\n")
    cat("  1. Re-run environmental_gwas_from_vcf.R (updated version saves these)\n")
    cat("  2. Or manually add to the save() command:\n")
    cat("     save(..., geno_mat_final, env_data, file='egwas_results.RData')\n")
    cat("  ═══════════════════════════════════════════════════════════\n\n")
  }
  
  # Proceed if we have both datasets
  if (has_geno && has_env) {
    
    cat("  Both genotype and environmental data found - creating visualizations\n")
    
    # Select top 15 pleiotropic markers and top 10 variables
    top_15_markers <- sig_per_marker %>%
      head(15) %>%
      pull(marker)
    
    top_10_variables <- sig_per_variable %>%
      head(10) %>%
      pull(variable)
    
    cat("  Top 15 markers:", paste(head(top_15_markers, 3), collapse = ", "), "...\n")
    cat("  Top 10 variables:", paste(head(top_10_variables, 3), collapse = ", "), "...\n")
    
    # Check if markers are in the genotype matrix
    markers_available <- top_15_markers[top_15_markers %in% colnames(geno_matrix)]
    
    if (length(markers_available) > 0) {
      
      cat("  Found", length(markers_available), "markers in genotype matrix\n")
      
      # Create a data frame with genotypes
      geno_data <- geno_matrix[, markers_available, drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("sample_id")
      
      # Check column names in environmental data
      cat("  Environmental data columns:", paste(head(colnames(env_dataset), 5), collapse = ", "), "...\n")
      
      # Identify the sample ID column
      id_col <- "sample_id"
      if (!"sample_id" %in% colnames(env_dataset)) {
        if ("individual" %in% colnames(env_dataset)) {
          id_col <- "individual"
          env_dataset <- env_dataset %>% dplyr::rename(sample_id = individual)
        } else {
          cat("  Warning: Cannot identify sample ID column in environmental data\n")
          cat("  Available columns:", paste(colnames(env_dataset), collapse = ", "), "\n")
        }
      }
      
      # Get environmental data for top variables that exist in the dataset
      vars_available <- top_10_variables[top_10_variables %in% colnames(env_dataset)]
      
      if (length(vars_available) == 0) {
        cat("  Warning: None of the top variables found in environmental data\n")
        cat("  Looking for:", paste(top_10_variables, collapse = ", "), "\n")
        cat("  Available:", paste(setdiff(colnames(env_dataset), "sample_id"), collapse = ", "), "\n")
      } else {
        cat("  Found", length(vars_available), "variables in environmental data\n")
      }
      
      # Get environmental data for available variables
      env_data_subset <- env_dataset %>%
        dplyr::select(sample_id, any_of(vars_available))
      
      # Merge genotype and environmental data
      merged_data <- geno_data %>%
        left_join(env_data_subset, by = "sample_id")
      
      cat("  Merged data dimensions:", nrow(merged_data), "samples ×", ncol(merged_data), "columns\n")
      
      # Check if we have any data to plot
      if (length(vars_available) > 0 && length(markers_available) > 0) {
        
        cat("  Creating plots with", length(markers_available), "markers ×", length(vars_available), "variables\n")
        
        cat("  Creating plots with", length(markers_available), "markers ×", length(vars_available), "variables\n")
        
        # Reshape for plotting
        plot_data <- merged_data %>%
          pivot_longer(
            cols = all_of(markers_available),
            names_to = "marker",
            values_to = "genotype"
          ) %>%
          pivot_longer(
            cols = all_of(vars_available),
            names_to = "variable",
            values_to = "env_value"
          ) %>%
          # Convert genotype to factor with meaningful labels
          mutate(
            genotype_class = case_when(
              genotype <= -0.5 ~ "Ref/Ref",
              genotype >= 0.5 ~ "Alt/Alt",
              TRUE ~ "Het"
            ),
            genotype_class = factor(genotype_class, 
                                    levels = c("Ref/Ref", "Het", "Alt/Alt"))
          ) %>%
          # Only keep marker-variable combinations that are significant
          dplyr::filter(!is.na(env_value), !is.na(genotype))
        
        # Filter to only significant combinations if we have sig_markers data
        if (nrow(plot_data) > 0) {
          plot_data <- plot_data %>%
            inner_join(
              sig_markers %>% 
                dplyr::filter(marker %in% markers_available, variable %in% vars_available) %>%
                dplyr::select(marker, variable, score),
              by = c("marker", "variable")
            )
          
          cat("  Plot data prepared:", nrow(plot_data), "observations\n")
        } else {
          cat("  Warning: No valid data for plotting\n")
        }
        
        if (nrow(plot_data) > 0) {
          
          cat("  Creating boxplot visualization...\n")
          
          # Create the boxplot
          p9 <- ggplot(plot_data, 
                       aes(x = genotype_class, y = env_value, fill = genotype_class)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.3, size = 0.5) +
            facet_grid(variable ~ marker, scales = "free_y", 
                       labeller = label_wrap_gen(width = 10)) +
            scale_fill_manual(
              name = "Genotype",
              values = c(
                "Ref/Ref" = "#3498db",   # Blue
                "Het" = "#95a5a6",        # Gray
                "Alt/Alt" = "#e74c3c"     # Red
              )
            ) +
            labs(
              title = "Environmental Values by Allelic State",
              subtitle = "Top 15 pleiotropic markers × Top 10 associated variables",
              x = "Genotype",
              y = "Environmental Value"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              axis.text.y = element_text(size = 7),
              strip.text.x = element_text(size = 7, angle = 90),
              strip.text.y = element_text(size = 7),
              legend.position = "bottom",
              panel.grid.minor = element_blank()
            )
          
          cat("  Calculating R² values...\n")
          
          # Calculate R² for each marker-variable combination
          r2_stats <- plot_data %>%
            group_by(marker, variable) %>%
            summarize(
              r2 = summary(lm(env_value ~ genotype))$r.squared,
              .groups = "drop"
            ) %>%
            mutate(
              label = paste0("R² = ", sprintf("%.3f", r2)),
              # Position for label (top right of each facet)
              x_pos = 2.5,
              y_pos = Inf
            )
          
          cat("  Adding R² labels to plot...\n")
          
          # Add R² labels to the plot
          p9 <- p9 + 
            geom_text(
              data = r2_stats,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE,
              hjust = 1, vjust = 1.5,
              size = 2.5,
              fontface = "bold",
              color = "black"
            )
          
          cat("  Saving full boxplot figure...\n")
          
          ggsave(
            filename = file.path(fig_dir, "09_allele_environment_boxplots.png"),
            plot = p9,
            width = 18,
            height = 14,
            dpi = 300
          )
          
          cat("  ✓ Allele-environment boxplots saved (with R² values)\n")
          
          cat("  Creating simplified version (top 5×5)...\n")
          
          # Also create a simplified version with just top 5 markers x top 5 variables
          plot_data_simple <- plot_data %>%
            dplyr::filter(
              marker %in% top_15_markers[1:5],
              variable %in% top_10_variables[1:5]
            )
          
          # Calculate R² for simplified version
          r2_stats_simple <- r2_stats %>%
            dplyr::filter(
              marker %in% top_15_markers[1:5],
              variable %in% top_10_variables[1:5]
            )
          
          p9_simple <- ggplot(plot_data_simple, 
                              aes(x = genotype_class, y = env_value, fill = genotype_class)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
            facet_grid(variable ~ marker, scales = "free_y") +
            scale_fill_manual(
              name = "Genotype",
              values = c(
                "Ref/Ref" = "#3498db",   # Blue
                "Het" = "#95a5a6",        # Gray
                "Alt/Alt" = "#e74c3c"     # Red
              )
            ) +
            labs(
              title = "Environmental Values by Allelic State (Simplified)",
              subtitle = "Top 5 pleiotropic markers × Top 5 associated variables",
              x = "Genotype",
              y = "Environmental Value"
            ) +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
              strip.text = element_text(size = 9),
              legend.position = "bottom",
              panel.grid.minor = element_blank()
            )
          
          # Add R² labels to simplified plot
          p9_simple <- p9_simple + 
            geom_text(
              data = r2_stats_simple,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE,
              hjust = 1, vjust = 1.5,
              size = 3.5,
              fontface = "bold",
              color = "black"
            )
          
          ggsave(
            filename = file.path(fig_dir, "09_allele_environment_boxplots_simple.png"),
            plot = p9_simple,
            width = 12,
            height = 10,
            dpi = 300
          )
          
          cat("  Simplified allele-environment boxplots saved (with R² values)\n")
          
          # Create summary statistics including R²
          allele_stats <- plot_data %>%
            group_by(marker, variable, genotype_class) %>%
            summarize(
              n = n(),
              mean_env = mean(env_value, na.rm = TRUE),
              sd_env = sd(env_value, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            left_join(r2_stats %>% dplyr::select(marker, variable, r2), 
                      by = c("marker", "variable")) %>%
            arrange(marker, variable, genotype_class)
          
          write.csv(
            allele_stats,
            file = file.path(fig_dir, "table_allele_environment_stats.csv"),
            row.names = FALSE
          )
          
          cat("  Allele-environment statistics table saved\n")
          
          
          # Create geographic map of allelic states for top 5 pleiotropic markers ----
          
          cat("  Creating geographic distribution maps...\n")
          cat("  Checking for coordinate columns in environmental data...\n")
          cat("    Available columns:", paste(colnames(env_dataset), collapse = ", "), "\n")
          
          # Check if we have latitude/longitude data (fully case-insensitive)
          coord_cols <- colnames(env_dataset)
          coord_cols_lower <- tolower(coord_cols)
          
          # Look for latitude variations (lat, latitude, etc.)
          lat_patterns <- c("^lat$", "^latitude$", "^lat_", "^latitude_")
          lat_idx <- which(sapply(lat_patterns, function(p) any(grepl(p, coord_cols_lower))))
          if (length(lat_idx) > 0) {
            lat_matches <- unique(unlist(lapply(lat_patterns[lat_idx], function(p) {
              coord_cols[grepl(p, coord_cols_lower)]
            })))
            lat_col <- lat_matches[1]
            has_lat <- TRUE
          } else {
            has_lat <- FALSE
          }
          
          # Look for longitude variations (lon, long, longitude, etc.)
          lon_patterns <- c("^lon$", "^long$", "^longitude$", "^lon_", "^long_", "^longitude_")
          lon_idx <- which(sapply(lon_patterns, function(p) any(grepl(p, coord_cols_lower))))
          if (length(lon_idx) > 0) {
            lon_matches <- unique(unlist(lapply(lon_patterns[lon_idx], function(p) {
              coord_cols[grepl(p, coord_cols_lower)]
            })))
            lon_col <- lon_matches[1]
            has_lon <- TRUE
          } else {
            has_lon <- FALSE
          }
          
          has_coords <- has_lat && has_lon
          
          if (has_coords) {
            cat("    ✓ Found latitude column:", lat_col, "\n")
            cat("    ✓ Found longitude column:", lon_col, "\n")
            
            # Rename to standard names for easier processing
            if (tolower(lat_col) != "latitude") {
              env_dataset <- env_dataset %>% dplyr::rename(latitude = !!sym(lat_col))
              cat("    Renamed '", lat_col, "' to 'latitude'\n", sep = "")
            } else if (lat_col != "latitude") {
              # Column is some variation of "latitude" but not exactly "latitude"
              env_dataset <- env_dataset %>% dplyr::rename(latitude = !!sym(lat_col))
              cat("    Renamed '", lat_col, "' to 'latitude'\n", sep = "")
            }
            
            if (tolower(lon_col) != "longitude") {
              env_dataset <- env_dataset %>% dplyr::rename(longitude = !!sym(lon_col))
              cat("    Renamed '", lon_col, "' to 'longitude'\n", sep = "")
            } else if (lon_col != "longitude") {
              # Column is some variation of "longitude" but not exactly "longitude"
              env_dataset <- env_dataset %>% dplyr::rename(longitude = !!sym(lon_col))
              cat("    Renamed '", lon_col, "' to 'longitude'\n", sep = "")
            }
            
            # Check if coordinates are valid
            coord_summary <- env_dataset %>%
              summarize(
                n_samples = n(),
                n_with_coords = sum(!is.na(latitude) & !is.na(longitude)),
                lat_range = paste(round(range(latitude, na.rm = TRUE), 2), collapse = " to "),
                lon_range = paste(round(range(longitude, na.rm = TRUE), 2), collapse = " to ")
              )
            
            cat("    Samples with coordinates:", coord_summary$n_with_coords, "of", coord_summary$n_samples, "\n")
            cat("    Latitude range:", coord_summary$lat_range, "\n")
            cat("    Longitude range:", coord_summary$lon_range, "\n")
            
            if (coord_summary$n_with_coords == 0) {
              cat("    ✗ No samples have valid coordinates\n")
              has_coords <- FALSE
            }
          } else {
            cat("    ✗ Latitude/longitude columns not found\n")
            if (!has_lat) {
              cat("      Missing: latitude column\n")
              cat("      Looked for: lat, Lat, LAT, latitude, Latitude, LATITUDE, etc.\n")
            }
            if (!has_lon) {
              cat("      Missing: longitude column\n")
              cat("      Looked for: lon, Lon, LON, long, Long, LONG, longitude, Longitude, LONGITUDE, etc.\n")
            }
          }
          
          if (has_coords) {
            
            # Get top 5 pleiotropic markers that are available
            top_5_markers <- markers_available[1:min(5, length(markers_available))]
            
            cat("  Using", length(top_5_markers), "markers for geographic maps\n")
            
            # Prepare data for mapping
            map_data <- geno_data %>%
              dplyr::select(sample_id, all_of(top_5_markers)) %>%
              pivot_longer(
                cols = all_of(top_5_markers),
                names_to = "marker",
                values_to = "genotype"
              ) %>%
              mutate(
                genotype_class = case_when(
                  genotype <= -0.5 ~ "Ref/Ref",
                  genotype >= 0.5 ~ "Alt/Alt",
                  TRUE ~ "Het"
                ),
                genotype_class = factor(genotype_class, 
                                        levels = c("Ref/Ref", "Het", "Alt/Alt"))
              ) %>%
              left_join(
                env_dataset %>% dplyr::select(sample_id, any_of(c("latitude", "longitude"))),
                by = "sample_id"
              ) %>%
              dplyr::filter(!is.na(latitude), !is.na(longitude))
            
            if (nrow(map_data) > 0) {
              cat("  Map data prepared:", nrow(map_data), "observations\n")
              
              cat("  Creating geographic map plot...\n")
              
              # Create the map
              p_map <- ggplot(map_data, aes(x = longitude, y = latitude, 
                                            color = genotype_class, shape = genotype_class)) +
                geom_point(size = 3, alpha = 0.8) +
                facet_wrap(~ marker, ncol = 3) +
                scale_color_manual(
                  name = "Genotype",
                  values = c(
                    "Ref/Ref" = "#3498db",   # Blue
                    "Het" = "#95a5a6",        # Gray
                    "Alt/Alt" = "#e74c3c"     # Red
                  )
                ) +
                scale_shape_manual(
                  name = "Genotype",
                  values = c("Ref/Ref" = 16, "Het" = 17, "Alt/Alt" = 15)
                ) +
                labs(
                  title = "Geographic Distribution of Allelic States",
                  subtitle = "Top 5 pleiotropic markers",
                  x = "Longitude",
                  y = "Latitude"
                ) +
                theme_bw() +
                theme(
                  legend.position = "bottom",
                  strip.text = element_text(size = 9, face = "bold"),
                  panel.grid.minor = element_blank()
                ) +
                coord_fixed(ratio = 1.3)  # Adjust aspect ratio for better map appearance
              
              cat("  Saving geographic map...\n")
              
              ggsave(
                filename = file.path(fig_dir, "11_geographic_allele_distribution.png"),
                plot = p_map,
                width = 14,
                height = 10,
                dpi = 300
              )
              
              cat("  ✓ Geographic allele distribution map saved\n")
              
              # Create an enhanced version with map background if mapdata is available
              if (requireNamespace("maps", quietly = TRUE)) {
                
                library(maps)
                
                # Get world map data
                world_map <- map_data("world")
                
                # Determine map bounds with some padding
                lon_range <- range(map_data$longitude, na.rm = TRUE)
                lat_range <- range(map_data$latitude, na.rm = TRUE)
                lon_padding <- diff(lon_range) * 0.1
                lat_padding <- diff(lat_range) * 0.1
                
                p_map_bg <- ggplot() +
                  geom_polygon(
                    data = world_map,
                    aes(x = long, y = lat, group = group),
                    fill = "gray95", color = "gray80", linewidth = 0.2
                  ) +
                  geom_point(
                    data = map_data,
                    aes(x = longitude, y = latitude, color = genotype_class, shape = genotype_class),
                    size = 3.5, alpha = 0.8
                  ) +
                  facet_wrap(~ marker, ncol = 3) +
                  scale_color_manual(
                    name = "Genotype",
                    values = c(
                      "Ref/Ref" = "#3498db",   # Blue
                      "Het" = "#95a5a6",        # Gray
                      "Alt/Alt" = "#e74c3c"     # Red
                    )
                  ) +
                  scale_shape_manual(
                    name = "Genotype",
                    values = c("Ref/Ref" = 16, "Het" = 17, "Alt/Alt" = 15)
                  ) +
                  coord_fixed(
                    ratio = 1.3,
                    xlim = c(lon_range[1] - lon_padding, lon_range[2] + lon_padding),
                    ylim = c(lat_range[1] - lat_padding, lat_range[2] + lat_padding)
                  ) +
                  labs(
                    title = "Geographic Distribution of Allelic States",
                    subtitle = "Top 5 pleiotropic markers with geographic context",
                    x = "Longitude",
                    y = "Latitude"
                  ) +
                  theme_bw() +
                  theme(
                    legend.position = "bottom",
                    strip.text = element_text(size = 9, face = "bold"),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "aliceblue")
                  )
                
                ggsave(
                  filename = file.path(fig_dir, "11_geographic_allele_distribution_with_map.png"),
                  plot = p_map_bg,
                  width = 14,
                  height = 10,
                  dpi = 300
                )
                
                cat("  Geographic allele distribution map with background saved\n")
              }
              
              # Calculate geographic patterns statistics
              geo_stats <- map_data %>%
                group_by(marker, genotype_class) %>%
                summarize(
                  n = n(),
                  mean_lat = mean(latitude, na.rm = TRUE),
                  sd_lat = sd(latitude, na.rm = TRUE),
                  mean_lon = mean(longitude, na.rm = TRUE),
                  sd_lon = sd(longitude, na.rm = TRUE),
                  .groups = "drop"
                )
              
              write.csv(
                geo_stats,
                file = file.path(fig_dir, "table_geographic_allele_stats.csv"),
                row.names = FALSE
              )
              
              cat("  Geographic allele statistics table saved\n")
              
            }  # Close if (nrow(map_data) > 0)
            else {
              cat("  Warning: No valid data for geographic maps (no samples with coordinates)\n")
            }
            
          } else {  # Close if (has_coords)
            cat("  Note: Latitude/longitude not found in environmental data\n")
            cat("  Skipping geographic allele distribution maps\n")
          }
          
        } else {  # Close if (length(vars_available) > 0 && length(markers_available) > 0)
          cat("  Warning: No valid marker-variable combinations for plotting\n")
        }
        
      } else {  # Close if (length(markers_available) > 0)
        cat("  Warning: Markers not found in genotype matrix\n")
        cat("  Available markers:", paste(head(colnames(geno_matrix), 5), collapse = ", "), "...\n")
      }
      
    } else {  # Close if (has_geno && has_env) and start else
      cat("  Note: Required data not found\n")
      if (!has_geno) cat("    - Missing: genotype matrix\n")
      if (!has_env) cat("    - Missing: environmental data\n")
      cat("  Skipping allele-environment boxplots and geographic maps\n")
      cat("  To generate these plots:\n")
      cat("    1. Ensure geno_mat_final and env_data are saved in egwas_results.RData\n")
      cat("    2. Or manually load them before running this script\n")
    }  # Close else block
    
  }  # Close if (has_geno && has_env) statement
  
}, error = function(e) {  # Close tryCatch main block, start error handler
  cat("\n")
  cat("  ✗✗✗ ERROR in Visualization 9 section ✗✗✗\n")
  cat("  Error message:", conditionMessage(e), "\n")
  cat("  Error location:", paste(deparse(conditionCall(e)), collapse = " "), "\n")
  cat("\n")
  cat("  This section creates figures 09 and 11.\n")
  cat("  The script will continue with remaining visualizations.\n")
  cat("\n")
  traceback()
}, warning = function(w) {
  cat("  ⚠ Warning in Visualization 9:", conditionMessage(w), "\n")
})

cat("  Visualization 9 section complete\n\n")


# Visualization 10: Effect size comparison (if alleles available) -------------

# This section would require allele frequency data
# Placeholder for future enhancement


# Visualization 11: Combined summary figure -----------------------------------

cat("Creating combined summary figure...\n")

# Create a multi-panel summary figure
summary_plot <- (p1 | p5) / (p3 | p4)

ggsave(
  filename = file.path(fig_dir, "10_combined_summary.png"),
  plot = summary_plot,
  width = 16,
  height = 12,
  dpi = 300
)

cat("  Combined summary figure saved\n")


# Export summary tables -------------------------------------------------------

cat("\nExporting summary tables...\n")

# Table 1: Summary by variable
write.csv(
  sig_per_variable,
  file = file.path(fig_dir, "table_summary_by_variable.csv"),
  row.names = FALSE
)

# Table 2: Summary by marker
write.csv(
  sig_per_marker,
  file = file.path(fig_dir, "table_summary_by_marker.csv"),
  row.names = FALSE
)

# Table 3: Chromosome distribution
write.csv(
  markers_per_chrom,
  file = file.path(fig_dir, "table_chromosome_distribution.csv"),
  row.names = FALSE
)

cat("  Summary tables exported\n")


# Generate HTML report --------------------------------------------------------

cat("\nGenerating HTML summary report...\n")

html_report <- sprintf('
<!DOCTYPE html>
<html>
<head>
    <title>Environmental GWAS Results Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .stat-box { 
            background-color: white; 
            padding: 20px; 
            margin: 20px 0; 
            border-radius: 8px; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .stat-row { display: flex; justify-content: space-around; margin: 20px 0; }
        .stat-item { 
            text-align: center; 
            padding: 20px;
            background-color: #ecf0f1;
            border-radius: 5px;
            flex: 1;
            margin: 0 10px;
        }
        .stat-value { font-size: 36px; font-weight: bold; color: #3498db; }
        .stat-label { font-size: 14px; color: #7f8c8d; margin-top: 5px; }
        img { max-width: 100%%; height: auto; margin: 20px 0; border-radius: 5px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        table { border-collapse: collapse; width: 100%%; margin: 20px 0; background-color: white; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .footer { margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; color: #7f8c8d; font-size: 12px; }
    </style>
</head>
<body>
    <h1>Environmental GWAS Results Summary</h1>
    
    <div class="stat-box">
        <h2>Overview Statistics</h2>
        <div class="stat-row">
            <div class="stat-item">
                <div class="stat-value">%d</div>
                <div class="stat-label">Total Associations</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">%d</div>
                <div class="stat-label">Unique Markers</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">%d</div>
                <div class="stat-label">Environmental Variables</div>
            </div>
            <div class="stat-item">
                <div class="stat-value">%d</div>
                <div class="stat-label">Chromosomes</div>
            </div>
        </div>
    </div>
    
    <div class="stat-box">
        <h2>Top Environmental Variables</h2>
        <table>
            <tr>
                <th>Variable</th>
                <th>Markers</th>
                <th>Associations</th>
                <th>Max Score</th>
            </tr>
            %s
        </table>
    </div>
    
    <div class="stat-box">
        <h2>Key Visualizations</h2>
        
        <h3>1. Markers per Variable</h3>
        <img src="01_markers_per_variable.png" alt="Markers per Variable">
        
        <h3>2. Pleiotropic Markers</h3>
        <img src="02_pleiotropic_markers.png" alt="Pleiotropic Markers">
        
        <h3>3. Score Distribution</h3>
        <img src="03_score_distribution.png" alt="Score Distribution">
        
        <h3>4. Chromosome Distribution</h3>
        <img src="05_markers_per_chromosome.png" alt="Chromosome Distribution">
        
        <h3>5. Association Heatmap</h3>
        <img src="07_association_heatmap.png" alt="Association Heatmap">
        
        <h3>6. Allele-Environment Associations (Simplified)</h3>
        <img src="09_allele_environment_boxplots_simple.png" alt="Allele-Environment Boxplots">
        <p style="font-size: 12px; color: #7f8c8d; margin-top: 10px;">
            Boxplots showing environmental values for different genotypes at top associated markers.
            Blue = Reference homozygote, Gray = Heterozygote, Red = Alternate homozygote.
            R² values show variance explained by each marker.
        </p>
        
        <h3>7. Geographic Allele Distribution</h3>
        <img src="11_geographic_allele_distribution.png" alt="Geographic Allele Distribution" onerror="this.style.display=\'none\'">
        <p style="font-size: 12px; color: #7f8c8d; margin-top: 10px;">
            Geographic distribution of genotypes for top 5 pleiotropic markers.
            Reveals spatial patterns of local adaptation.
        </p>
        
        <h3>8. Combined Summary</h3>
        <img src="10_combined_summary.png" alt="Combined Summary">
    </div>
    
    <div class="footer">
        Generated: %s<br>
        FDR Threshold: %.2f<br>
        Input file: %s
    </div>
</body>
</html>
',
                       nrow(sig_markers),
                       n_distinct(sig_markers$marker),
                       n_distinct(sig_markers$variable),
                       n_distinct(sig_markers$chrom),
                       paste(
                         apply(head(sig_per_variable, 5), 1, function(row) {
                           sprintf(
                             "<tr><td>%s</td><td>%s</td><td>%s</td><td>%.2f</td></tr>",
                             row["variable"], row["n_markers"], row["n_associations"], as.numeric(row["max_score"])
                           )
                         }),
                         collapse = "\n            "
                       ),
                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                       fdr_threshold,
                       basename(sig_markers_file)
)


writeLines(html_report, file.path(fig_dir, "summary_report.html"))

cat("  HTML report generated: ", file.path(fig_dir, "summary_report.html"), "\n")


# Print summary ---------------------------------------------------------------

cat("\n" , rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nSummary:\n")
cat("  Total associations:", nrow(sig_markers), "\n")
cat("  Unique markers:", n_distinct(sig_markers$marker), "\n")
cat("  Environmental variables:", n_distinct(sig_markers$variable), "\n")
cat("  Chromosomes represented:", n_distinct(sig_markers$chrom), "\n")
cat("\nOutput files saved to:", fig_dir, "\n")
cat("  - Up to 12 visualization figures (PNG)\n")
cat("  - Up to 5 summary tables (CSV)\n")
cat("  - 1 HTML report\n")
cat("\nKey visualizations:\n")
cat("  - Basic summaries (01-05)\n")
cat("  - Marker positions (06)\n")
cat("  - Association patterns (07-08)\n")
cat("  - Allele-environment effects (09)\n")
cat("  - Combined summary (10)\n")
cat("  - Geographic allele distribution (11, if lat/long available)\n")
cat("\nTo view results:\n")
cat("  - Open:", file.path(fig_dir, "summary_report.html"), "\n")
cat(rep("=", 70), "\n\n", sep = "")