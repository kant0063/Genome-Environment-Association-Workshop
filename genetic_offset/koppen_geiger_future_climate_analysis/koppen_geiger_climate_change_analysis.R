# ============================================================================
# Calculate Köppen-Geiger Climate Class Changes Between Time Periods
# ============================================================================
# This script extracts Köppen-Geiger climate classifications for sample 
# locations from a CSV file and calculates changes between 1991-2020 and 
# 2041-2070 under SSP585
#
# Data source: Beck et al. (2023) - High-resolution (1 km) Köppen-Geiger maps
# https://doi.org/10.1038/s41597-023-02549-6
# ============================================================================

# Load required packages
if (!require("terra", quietly = TRUE)) {
  stop("Package 'terra' is required but not installed. Install with: install.packages('terra')")
}
if (!require("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required but not installed. Install with: install.packages('dplyr')")
}
if (!require("tidyr", quietly = TRUE)) {
  stop("Package 'tidyr' is required but not installed. Install with: install.packages('tidyr')")
}
if (!require("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required but not installed. Install with: install.packages('ggplot2')")
}
if (!require("readr", quietly = TRUE)) {
  stop("Package 'readr' is required but not installed. Install with: install.packages('readr')")
}
if (!require("viridis", quietly = TRUE)) {
  stop("Package 'viridis' is required but not installed. Install with: install.packages('viridis')")
}

# Explicitly load terra to ensure methods are available
library(terra)

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# 1. Path to your input CSV file
# Your CSV should contain at minimum: longitude, latitude, and population columns
# Example columns: sample_id, longitude, latitude, population
input_csv <- "your_samples.csv"

# 2. Column names in your CSV (adjust these to match your actual column names)
lon_column <- "longitude"       # Name of longitude column
lat_column <- "latitude"        # Name of latitude column
population_column <- "population"  # Name of population/grouping column
id_column <- "sample_id"        # Optional: Name of sample ID column (set to NULL if none)

# 3. Paths to the Köppen-Geiger rasters
# Download from: https://figshare.com/articles/dataset/21789074
# File: koppen_geiger_tif.zip (125 MB)
# Note: The baseline period closest to 1991-2000 in the dataset is 1991-2020
kg_historical_path <- "koppen_geiger/1991_2020/koppen_geiger_0p00833333.tif"
kg_future_path <- "koppen_geiger/2041_2070/ssp585/koppen_geiger_0p00833333.tif"

# 4. Output file prefix (optional)
output_prefix <- "koppen_geiger_analysis"

# ============================================================================
# FUNCTION: Validate CSV columns
# ============================================================================
validate_csv <- function(data, lon_col, lat_col, pop_col, id_col = NULL) {
  required_cols <- c(lon_col, lat_col, pop_col)
  if (!is.null(id_col)) {
    required_cols <- c(id_col, required_cols)
  }
  
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("ERROR: Missing required columns in CSV:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Check for NA values in coordinates
  n_missing_lon <- sum(is.na(data[[lon_col]]))
  n_missing_lat <- sum(is.na(data[[lat_col]]))
  
  if (n_missing_lon > 0 || n_missing_lat > 0) {
    warning(paste("WARNING:", n_missing_lon, "rows with missing longitude,", 
                  n_missing_lat, "rows with missing latitude"))
  }
  
  # Check coordinate ranges
  lon_range <- range(data[[lon_col]], na.rm = TRUE)
  lat_range <- range(data[[lat_col]], na.rm = TRUE)
  
  if (lon_range[1] < -180 || lon_range[2] > 180) {
    warning("WARNING: Longitude values outside valid range (-180 to 180)")
  }
  
  if (lat_range[1] < -90 || lat_range[2] > 90) {
    warning("WARNING: Latitude values outside valid range (-90 to 90)")
  }
  
  cat("CSV validation passed!\n")
  cat("  Rows:", nrow(data), "\n")
  cat("  Longitude range:", round(lon_range[1], 3), "to", round(lon_range[2], 3), "\n")
  cat("  Latitude range:", round(lat_range[1], 3), "to", round(lat_range[2], 3), "\n")
  cat("  Unique populations:", length(unique(data[[pop_col]])), "\n\n")
  
  return(TRUE)
}

# ============================================================================
# KÖPPEN-GEIGER CLASSIFICATION LEGEND
# ============================================================================
# Create a lookup table for Köppen-Geiger classes
# Based on legend.txt from Beck et al. (2023)
kg_legend <- data.frame(
  code = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
           17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
  class = c("Ocean", "Af", "Am", "Aw", "BWh", "BWk", "BSh", "BSk", "Csa", 
            "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc", 
            "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", 
            "Dfa", "Dfb", "Dfc", "Dfd", "ET", "EF"),
  description = c(
    "Ocean",
    "Tropical, rainforest",
    "Tropical, monsoon",
    "Tropical, savannah",
    "Arid, desert, hot",
    "Arid, desert, cold",
    "Arid, steppe, hot",
    "Arid, steppe, cold",
    "Temperate, dry summer, hot summer",
    "Temperate, dry summer, warm summer",
    "Temperate, dry summer, cold summer",
    "Temperate, dry winter, hot summer",
    "Temperate, dry winter, warm summer",
    "Temperate, dry winter, cold summer",
    "Temperate, no dry season, hot summer",
    "Temperate, no dry season, warm summer",
    "Temperate, no dry season, cold summer",
    "Cold, dry summer, hot summer",
    "Cold, dry summer, warm summer",
    "Cold, dry summer, cold summer",
    "Cold, dry summer, very cold winter",
    "Cold, dry winter, hot summer",
    "Cold, dry winter, warm summer",
    "Cold, dry winter, cold summer",
    "Cold, dry winter, very cold winter",
    "Cold, no dry season, hot summer",
    "Cold, no dry season, warm summer",
    "Cold, no dry season, cold summer",
    "Cold, no dry season, very cold winter",
    "Polar, tundra",
    "Polar, frost"
  ),
  major_class = c(
    "Ocean",
    "A", "A", "A",  # Tropical
    "B", "B", "B", "B",  # Arid
    "C", "C", "C", "C", "C", "C", "C", "C", "C",  # Temperate
    "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D",  # Cold
    "E", "E"  # Polar
  )
)

# ============================================================================
# MAIN ANALYSIS
# ============================================================================
cat(strrep("=", 78), "\n")
cat("Köppen-Geiger Climate Change Analysis\n")
cat(strrep("=", 78), "\n\n")

# ----------------------------------------------------------------------------
# STEP 1: Load and validate input CSV
# ----------------------------------------------------------------------------
cat("STEP 1: Loading input CSV file...\n")
cat("File:", input_csv, "\n\n")

if (!file.exists(input_csv)) {
  stop("ERROR: Input CSV file not found: ", input_csv)
}

# Read the CSV file
sample_data <- read_csv(input_csv, show_col_types = FALSE)

# Validate the CSV
validate_csv(sample_data, lon_column, lat_column, population_column, id_column)

# Standardize column names for easier processing
sample_data <- sample_data %>%
  rename(
    longitude = !!sym(lon_column),
    latitude = !!sym(lat_column),
    population = !!sym(population_column)
  )

# Add sample_id if not provided
if (is.null(id_column) || !id_column %in% names(sample_data)) {
  sample_data$sample_id <- 1:nrow(sample_data)
  cat("Note: No sample ID column found. Created sequential IDs.\n\n")
} else {
  sample_data <- sample_data %>%
    rename(sample_id = !!sym(id_column))
}

# Remove rows with missing coordinates
n_before <- nrow(sample_data)
sample_data <- sample_data %>%
  filter(!is.na(longitude) & !is.na(latitude))
n_after <- nrow(sample_data)

if (n_before != n_after) {
  cat("Removed", n_before - n_after, "rows with missing coordinates\n\n")
}

# ----------------------------------------------------------------------------
# STEP 2: Load Köppen-Geiger rasters
# ----------------------------------------------------------------------------
cat("STEP 2: Loading Köppen-Geiger rasters...\n")

if (!file.exists(kg_historical_path)) {
  stop("ERROR: Historical Köppen-Geiger raster not found: ", kg_historical_path)
}

if (!file.exists(kg_future_path)) {
  stop("ERROR: Future Köppen-Geiger raster not found: ", kg_future_path)
}

# Load historical (1991-2020) and future (2041-2070, SSP585) rasters
kg_historical <- terra::rast(kg_historical_path)
kg_future <- terra::rast(kg_future_path)

cat("  Historical raster (1991-2020):", nrow(kg_historical), "x", ncol(kg_historical), "pixels\n")
cat("  Future raster (2041-2070, SSP585):", nrow(kg_future), "x", ncol(kg_future), "pixels\n\n")

# ----------------------------------------------------------------------------
# STEP 3: Extract Köppen-Geiger classes for sample locations
# ----------------------------------------------------------------------------
cat("STEP 3: Extracting Köppen-Geiger classes for", nrow(sample_data), "locations...\n")

# Create spatial points from sample data
# Use try-catch to handle potential CRS issues
sample_points <- tryCatch({
  terra::vect(sample_data, 
              geom = c("longitude", "latitude"), 
              crs = "EPSG:4326")
}, warning = function(w) {
  cat("  Warning during spatial point creation:", conditionMessage(w), "\n")
  cat("  Attempting to create points without explicit CRS...\n")
  # Try without explicit CRS
  terra::vect(sample_data, 
              geom = c("longitude", "latitude"))
}, error = function(e) {
  stop("Error creating spatial points: ", conditionMessage(e))
})

# Extract historical Köppen-Geiger classes
cat("  Extracting historical climate classes...\n")
kg_hist_values <- terra::extract(kg_historical, sample_points)
sample_data$kg_historical_code <- kg_hist_values[, 2]

# Extract future Köppen-Geiger classes
cat("  Extracting future climate classes...\n")
kg_future_values <- terra::extract(kg_future, sample_points)
sample_data$kg_future_code <- kg_future_values[, 2]

cat("  Extraction complete!\n\n")

# ----------------------------------------------------------------------------
# STEP 4: Add climate class labels
# ----------------------------------------------------------------------------
cat("STEP 4: Adding climate class labels...\n")

# Merge with legend to get class names and descriptions
sample_data <- sample_data %>%
  left_join(
    kg_legend %>% 
      select(code, class, description, major_class) %>% 
      rename(
        kg_historical_code = code,
        kg_historical_class = class,
        kg_historical_description = description,
        kg_historical_major = major_class
      ),
    by = "kg_historical_code"
  ) %>%
  left_join(
    kg_legend %>% 
      select(code, class, description, major_class) %>% 
      rename(
        kg_future_code = code,
        kg_future_class = class,
        kg_future_description = description,
        kg_future_major = major_class
      ),
    by = "kg_future_code"
  )

cat("  Labels added!\n\n")

# ----------------------------------------------------------------------------
# STEP 5: Calculate climate change metrics
# ----------------------------------------------------------------------------
cat("STEP 5: Calculating climate change metrics...\n")

# Determine if climate class changed
sample_data <- sample_data %>%
  mutate(
    climate_changed = kg_historical_class != kg_future_class,
    major_class_changed = kg_historical_major != kg_future_major,
    change_type = case_when(
      kg_historical_class == kg_future_class ~ "No change",
      TRUE ~ paste0(kg_historical_class, " → ", kg_future_class)
    ),
    major_change_type = case_when(
      kg_historical_major == kg_future_major ~ "No change",
      TRUE ~ paste0(kg_historical_major, " → ", kg_future_major)
    )
  )

cat("  Metrics calculated!\n\n")

# ----------------------------------------------------------------------------
# STEP 6: Create spatial maps
# ----------------------------------------------------------------------------
cat("STEP 6: Creating spatial maps...\n")

# Create a bounding box around sample points with buffer
lon_range <- range(sample_data$longitude, na.rm = TRUE)
lat_range <- range(sample_data$latitude, na.rm = TRUE)

# Add 10% buffer to bounding box
lon_buffer <- diff(lon_range) * 0.1
lat_buffer <- diff(lat_range) * 0.1

bbox_ext <- terra::ext(
  lon_range[1] - lon_buffer,
  lon_range[2] + lon_buffer,
  lat_range[1] - lat_buffer,
  lat_range[2] + lat_buffer
)

# Crop rasters to study area
cat("  Cropping rasters to study region...\n")
kg_hist_crop <- terra::crop(kg_historical, bbox_ext)
kg_fut_crop <- terra::crop(kg_future, bbox_ext)

# Convert rasters to data frames for ggplot
cat("  Converting rasters to data frames...\n")
kg_hist_df <- as.data.frame(kg_hist_crop, xy = TRUE)
kg_fut_df <- as.data.frame(kg_fut_crop, xy = TRUE)

# Rename the value columns
names(kg_hist_df)[3] <- "kg_code"
names(kg_fut_df)[3] <- "kg_code"

# Add climate class labels
kg_hist_df <- kg_hist_df %>%
  left_join(kg_legend %>% select(code, class, description) %>% rename(kg_code = code),
            by = "kg_code")

kg_fut_df <- kg_fut_df %>%
  left_join(kg_legend %>% select(code, class, description) %>% rename(kg_code = code),
            by = "kg_code")

# Create color palette for Köppen-Geiger classes
# Using standard Köppen-Geiger colors
kg_colors <- c(
  "Ocean" = "#96C8FF",
  "Af" = "#960000", "Am" = "#FF0000", "Aw" = "#FFC8C8",  # Tropical (red)
  "BWh" = "#FFFF64", "BWk" = "#FFDC64", "BSh" = "#F0B432", "BSk" = "#F0DC64",  # Arid (yellow/tan)
  "Csa" = "#FFFF00", "Csb" = "#C8C800", "Csc" = "#969600",  # Temperate dry summer (yellow-green)
  "Cwa" = "#96FF96", "Cwb" = "#64C864", "Cwc" = "#329632",  # Temperate dry winter (green)
  "Cfa" = "#C8FF50", "Cfb" = "#64FF50", "Cfc" = "#32C800",  # Temperate no dry season (bright green)
  "Dsa" = "#FF6EFF", "Dsb" = "#C85AC8", "Dsc" = "#963C96", "Dsd" = "#966496",  # Cold dry summer (magenta)
  "Dwa" = "#ABD6FF", "Dwb" = "#6EA6FF", "Dwc" = "#4678C8", "Dwd" = "#326496",  # Cold dry winter (blue)
  "Dfa" = "#00FFFF", "Dfb" = "#50C8F0", "Dfc" = "#0096C8", "Dfd" = "#006496",  # Cold no dry season (cyan)
  "ET" = "#B4B4B4", "EF" = "#646464"  # Polar (gray)
)

# Map 1: Historical climate (1991-2020)
cat("  Creating historical climate map...\n")
p_hist_map <- ggplot() +
  geom_raster(data = kg_hist_df %>% filter(!is.na(class)), 
              aes(x = x, y = y, fill = class)) +
  geom_point(data = sample_data, 
             aes(x = longitude, y = latitude, shape = population),
             size = 2.5, color = "black", stroke = 1.2) +
  geom_point(data = sample_data, 
             aes(x = longitude, y = latitude, color = population),
             size = 1.5) +
  scale_fill_manual(values = kg_colors, 
                    name = "Köppen-Geiger\nClimate Class",
                    na.value = "white") +
  scale_color_viridis_d(name = "Population") +
  scale_shape_manual(values = 15:(15 + length(unique(sample_data$population)) - 1),
                     name = "Population") +
  coord_fixed() +
  labs(title = "A) Historical Climate Classification (1991-2020)",
       subtitle = paste0("Study region with ", nrow(sample_data), " sample locations"),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NA))

ggsave(paste0(output_prefix, "_map_A_historical.png"), p_hist_map, 
       width = 14, height = 10, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_map_A_historical.png\n"))

# Map 2: Future climate (2041-2070, SSP585)
cat("  Creating future climate map...\n")
p_fut_map <- ggplot() +
  geom_raster(data = kg_fut_df %>% filter(!is.na(class)), 
              aes(x = x, y = y, fill = class)) +
  geom_point(data = sample_data, 
             aes(x = longitude, y = latitude, shape = population),
             size = 2.5, color = "black", stroke = 1.2) +
  geom_point(data = sample_data, 
             aes(x = longitude, y = latitude, color = population),
             size = 1.5) +
  scale_fill_manual(values = kg_colors, 
                    name = "Köppen-Geiger\nClimate Class",
                    na.value = "white") +
  scale_color_viridis_d(name = "Population") +
  scale_shape_manual(values = 15:(15 + length(unique(sample_data$population)) - 1),
                     name = "Population") +
  coord_fixed() +
  labs(title = "B) Future Climate Classification (2041-2070, SSP585)",
       subtitle = paste0("Projected climate change scenario with ", nrow(sample_data), " sample locations"),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NA))

ggsave(paste0(output_prefix, "_map_B_future.png"), p_fut_map, 
       width = 14, height = 10, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_map_B_future.png\n"))

# Map 3: Change map showing areas that changed climate class
cat("  Creating climate change map...\n")

# Calculate change for the raster
kg_change_rast <- kg_hist_crop != kg_fut_crop
kg_change_df <- as.data.frame(kg_change_rast, xy = TRUE)
names(kg_change_df)[3] <- "changed"

# Add sample data with change information
sample_data_map <- sample_data %>%
  mutate(change_label = ifelse(climate_changed, 
                               paste0(kg_historical_class, " → ", kg_future_class),
                               "No change"))

p_change_map <- ggplot() +
  geom_raster(data = kg_change_df %>% filter(!is.na(changed)), 
              aes(x = x, y = y, fill = changed)) +
  geom_point(data = sample_data_map, 
             aes(x = longitude, y = latitude, shape = climate_changed),
             size = 3, color = "black", stroke = 1.5) +
  geom_point(data = sample_data_map, 
             aes(x = longitude, y = latitude, color = climate_changed),
             size = 2) +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#95a5a6"),
                    name = "Climate Class\nChanged",
                    labels = c("No change", "Changed")) +
  scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#7f8c8d"),
                     name = "Sample\nChanged",
                     labels = c("No", "Yes")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                     name = "Sample\nChanged",
                     labels = c("No", "Yes")) +
  coord_fixed() +
  labs(title = "Climate Change Map: 1991-2020 to 2041-2070 (SSP585)",
       subtitle = paste0("Red areas = climate class changed; Gray areas = no change | ",
                        "Triangles = samples that changed; Circles = no change"),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NA))

ggsave(paste0(output_prefix, "_map_change.png"), p_change_map, 
       width = 14, height = 10, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_map_change.png\n"))

# Map 4: Side-by-side comparison (combined A/B figure)
cat("  Creating side-by-side comparison map...\n")

# Prepare data for faceting
kg_hist_df_labeled <- kg_hist_df %>%
  filter(!is.na(class)) %>%
  mutate(period = "A) Historical (1991-2020)")

kg_fut_df_labeled <- kg_fut_df %>%
  filter(!is.na(class)) %>%
  mutate(period = "B) Future (2041-2070, SSP585)")

kg_combined_df <- bind_rows(kg_hist_df_labeled, kg_fut_df_labeled)

sample_data_facet <- bind_rows(
  sample_data %>% mutate(period = "A) Historical (1991-2020)"),
  sample_data %>% mutate(period = "B) Future (2041-2070, SSP585)")
)

p_combined_map <- ggplot() +
  geom_raster(data = kg_combined_df, 
              aes(x = x, y = y, fill = class)) +
  geom_point(data = sample_data_facet, 
             aes(x = longitude, y = latitude),
             size = 2, color = "black", shape = 21, fill = "yellow", stroke = 1) +
  scale_fill_manual(values = kg_colors, 
                    name = "Climate Class",
                    na.value = "white") +
  coord_fixed() +
  facet_wrap(~period, ncol = 2) +
  labs(title = "Köppen-Geiger Climate Classification: Historical vs Future",
       subtitle = paste0("Study region with ", nrow(sample_data), " sample locations (yellow points)"),
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        strip.text = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "lightblue", color = NA))

ggsave(paste0(output_prefix, "_map_AB_combined.png"), p_combined_map, 
       width = 16, height = 8, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_map_AB_combined.png\n\n"))

# ----------------------------------------------------------------------------
# STEP 7: Summarize results
# ----------------------------------------------------------------------------
cat("STEP 7: Summarizing results...\n\n")

# Overall summary
overall_summary <- sample_data %>%
  summarize(
    total_samples = n(),
    samples_changed = sum(climate_changed, na.rm = TRUE),
    pct_changed = round(100 * samples_changed / total_samples, 1),
    major_class_changed = sum(major_class_changed, na.rm = TRUE),
    pct_major_changed = round(100 * major_class_changed / total_samples, 1)
  )

cat(strrep("=", 60), "\n")
cat("OVERALL SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("Total samples:", overall_summary$total_samples, "\n")
cat("Samples with climate change:", overall_summary$samples_changed, 
    "(", overall_summary$pct_changed, "%)\n")
cat("Major class changes:", overall_summary$major_class_changed,
    "(", overall_summary$pct_major_changed, "%)\n")
cat(strrep("=", 60), "\n\n")

# Summary by population
population_summary <- sample_data %>%
  group_by(population) %>%
  summarize(
    n_samples = n(),
    n_changed = sum(climate_changed, na.rm = TRUE),
    pct_changed = round(100 * n_changed / n_samples, 1),
    n_major_changed = sum(major_class_changed, na.rm = TRUE),
    pct_major_changed = round(100 * n_major_changed / n_samples, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_changed))

cat("SUMMARY BY POPULATION\n")
print(population_summary, n = Inf)
cat("\n")

# Detailed change summary by population
change_detail <- sample_data %>%
  filter(climate_changed) %>%
  group_by(population, change_type) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(population, desc(n))

cat("DETAILED CLIMATE TRANSITIONS BY POPULATION\n")
print(change_detail, n = 30)
cat("\n")

# Most common climate transitions overall
top_transitions <- sample_data %>%
  filter(climate_changed) %>%
  count(change_type, sort = TRUE) %>%
  head(15)

cat("TOP 15 MOST COMMON CLIMATE TRANSITIONS\n")
print(top_transitions, n = Inf)
cat("\n")

# ----------------------------------------------------------------------------
# STEP 8: Create visualizations
# ----------------------------------------------------------------------------
cat("STEP 8: Creating visualizations...\n")

# 1. Bar plot of climate change by population
p1 <- ggplot(population_summary, aes(x = reorder(population, pct_changed), y = pct_changed)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(pct_changed, "%")), 
            hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "Percentage of Samples with Climate Class Change by Population",
       subtitle = "1991-2020 to 2041-2070 (SSP585)",
       x = "Population",
       y = "% Samples with Climate Change") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(paste0(output_prefix, "_by_population.png"), p1, width = 10, height = 6, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_by_population.png\n"))

# 2. Stacked bar chart showing climate transitions
change_counts <- sample_data %>%
  count(population, climate_changed) %>%
  mutate(climate_changed = ifelse(climate_changed, "Changed", "Unchanged"))

p2 <- ggplot(change_counts, aes(x = reorder(population, -n), y = n, fill = climate_changed)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Changed" = "#e74c3c", "Unchanged" = "#95a5a6")) +
  labs(title = "Climate Class Changes by Population",
       subtitle = "1991-2020 to 2041-2070 (SSP585)",
       x = "Population",
       y = "Number of Samples",
       fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, face = "bold"))

ggsave(paste0(output_prefix, "_stacked_bar.png"), p2, width = 10, height = 6, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_stacked_bar.png\n"))

# 3. Most common climate transitions
if(nrow(top_transitions) > 0) {
  p3 <- ggplot(top_transitions, aes(x = reorder(change_type, n), y = n)) +
    geom_bar(stat = "identity", fill = "#3498db") +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    coord_flip() +
    labs(title = "Most Common Climate Transitions",
         subtitle = "1991-2020 to 2041-2070 (SSP585)",
         x = "Climate Transition",
         y = "Number of Samples") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"))
  
  ggsave(paste0(output_prefix, "_top_transitions.png"), p3, width = 10, height = 8, dpi = 300)
  cat("  Saved:", paste0(output_prefix, "_top_transitions.png\n"))
}

# 4. Major climate class changes
major_class_summary <- sample_data %>%
  group_by(population) %>%
  summarize(
    n_samples = n(),
    n_major_changed = sum(major_class_changed, na.rm = TRUE),
    pct_major_changed = round(100 * n_major_changed / n_samples, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_major_changed))

p4 <- ggplot(major_class_summary, aes(x = reorder(population, pct_major_changed), 
                                       y = pct_major_changed)) +
  geom_bar(stat = "identity", fill = "#e67e22") +
  geom_text(aes(label = paste0(pct_major_changed, "%")), 
            hjust = -0.2, size = 3.5) +
  coord_flip() +
  labs(title = "Major Climate Class Changes by Population",
       subtitle = "1991-2020 to 2041-2070 (SSP585) | A=Tropical, B=Arid, C=Temperate, D=Cold, E=Polar",
       x = "Population",
       y = "% Samples with Major Class Change") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave(paste0(output_prefix, "_major_classes.png"), p4, width = 10, height = 6, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_major_classes.png\n"))

# 5. Detailed minor climate class transitions by population
# Create a heatmap or stacked bar showing specific transitions
transition_matrix <- sample_data %>%
  filter(climate_changed) %>%
  count(kg_historical_class, kg_future_class) %>%
  arrange(desc(n))

# Top 20 transitions for visualization
top_20_transitions <- transition_matrix %>%
  head(20)

if(nrow(top_20_transitions) > 0) {
  p5 <- ggplot(top_20_transitions, 
               aes(x = reorder(paste0(kg_historical_class, " → ", kg_future_class), n), 
                   y = n)) +
    geom_bar(stat = "identity", fill = "#9b59b6") +
    geom_text(aes(label = n), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(title = "Top 20 Minor Climate Class Transitions",
         subtitle = "1991-2020 to 2041-2070 (SSP585) | Detailed Köppen-Geiger Classes",
         x = "Climate Transition",
         y = "Number of Samples") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 9))
  
  ggsave(paste0(output_prefix, "_minor_transitions_top20.png"), p5, 
         width = 10, height = 10, dpi = 300)
  cat("  Saved:", paste0(output_prefix, "_minor_transitions_top20.png\n"))
}

# 6. Minor class transitions by population - faceted plot
pop_transitions <- sample_data %>%
  filter(climate_changed) %>%
  group_by(population, kg_historical_class, kg_future_class) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(population) %>%
  arrange(desc(n)) %>%
  slice_head(n = 5) %>%  # Top 5 transitions per population
  ungroup() %>%
  mutate(transition = paste0(kg_historical_class, " → ", kg_future_class))

if(nrow(pop_transitions) > 0) {
  p6 <- ggplot(pop_transitions, 
               aes(x = reorder(transition, n), y = n, fill = population)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~population, scales = "free_y", ncol = 2) +
    labs(title = "Top 5 Minor Climate Transitions per Population",
         subtitle = "1991-2020 to 2041-2070 (SSP585)",
         x = "Climate Transition",
         y = "Number of Samples") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          strip.text = element_text(face = "bold"))
  
  ggsave(paste0(output_prefix, "_minor_transitions_by_pop.png"), p6, 
         width = 12, height = 8, dpi = 300)
  cat("  Saved:", paste0(output_prefix, "_minor_transitions_by_pop.png\n"))
}

# 7. Transition matrix heatmap (if not too many classes)
# Count unique classes in the data
unique_hist_classes <- unique(sample_data$kg_historical_class[sample_data$climate_changed])
unique_fut_classes <- unique(sample_data$kg_future_class[sample_data$climate_changed])

# Only create heatmap if reasonable number of classes (< 15 each)
if(length(unique_hist_classes) <= 15 && length(unique_fut_classes) <= 15) {
  transition_heatmap_data <- sample_data %>%
    filter(climate_changed) %>%
    count(kg_historical_class, kg_future_class)
  
  p7 <- ggplot(transition_heatmap_data, 
               aes(x = kg_future_class, y = kg_historical_class, fill = n)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = n), color = "white", size = 3) +
    scale_fill_gradient(low = "#3498db", high = "#e74c3c", name = "Count") +
    labs(title = "Climate Class Transition Matrix",
         subtitle = "1991-2020 to 2041-2070 (SSP585)",
         x = "Future Climate Class (2041-2070)",
         y = "Historical Climate Class (1991-2020)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  
  ggsave(paste0(output_prefix, "_transition_heatmap.png"), p7, 
         width = 10, height = 8, dpi = 300)
  cat("  Saved:", paste0(output_prefix, "_transition_heatmap.png\n"))
}

# 8. Distribution of climate classes - before and after
class_distribution <- sample_data %>%
  select(kg_historical_class, kg_future_class) %>%
  pivot_longer(cols = everything(), 
               names_to = "period", 
               values_to = "class") %>%
  mutate(period = ifelse(period == "kg_historical_class", 
                        "Historical (1991-2020)", 
                        "Future (2041-2070)")) %>%
  count(period, class) %>%
  group_by(period) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

# Only show classes present in either period
p8 <- ggplot(class_distribution, 
             aes(x = reorder(class, -n), y = n, fill = period)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(percentage, 1)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = c("Historical (1991-2020)" = "#3498db", 
                                "Future (2041-2070)" = "#e67e22")) +
  labs(title = "Climate Class Distribution: Historical vs Future",
       subtitle = "1991-2020 to 2041-2070 (SSP585) | Numbers show percentage of total",
       x = "Climate Class",
       y = "Number of Samples",
       fill = "Period") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

ggsave(paste0(output_prefix, "_class_distribution.png"), p8, 
       width = 12, height = 7, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_class_distribution.png\n"))

# 9. Alluvial/Sankey-style plot showing transitions (simplified version)
# Show top N historical classes and their transitions
top_n_classes <- 8  # Show top 8 most common historical classes

top_classes_data <- sample_data %>%
  count(kg_historical_class, sort = TRUE) %>%
  head(top_n_classes) %>%
  pull(kg_historical_class)

transition_flow <- sample_data %>%
  filter(kg_historical_class %in% top_classes_data) %>%
  count(kg_historical_class, kg_future_class, climate_changed) %>%
  mutate(transition_label = ifelse(climate_changed, 
                                   paste0(kg_historical_class, " → ", kg_future_class),
                                   "No change"))

p9 <- ggplot(transition_flow, 
             aes(x = kg_historical_class, y = n, fill = kg_future_class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = paste0("Climate Class Transitions (Top ", top_n_classes, " Historical Classes)"),
       subtitle = "1991-2020 to 2041-2070 (SSP585) | Stacked bars show destination classes",
       x = "Historical Climate Class (1991-2020)",
       y = "Number of Samples",
       fill = "Future Class\n(2041-2070)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  scale_fill_viridis_d(option = "turbo")

ggsave(paste0(output_prefix, "_transition_flow.png"), p9, 
       width = 12, height = 7, dpi = 300)
cat("  Saved:", paste0(output_prefix, "_transition_flow.png\n\n"))

# ----------------------------------------------------------------------------
# STEP 9: Export results
# ----------------------------------------------------------------------------
cat("STEP 9: Exporting results...\n")

# Export full results with all columns
write_csv(sample_data, paste0(output_prefix, "_full_results.csv"))
cat("  Saved:", paste0(output_prefix, "_full_results.csv\n"))

# Export summary by population
write_csv(population_summary, paste0(output_prefix, "_population_summary.csv"))
cat("  Saved:", paste0(output_prefix, "_population_summary.csv\n"))

# Export detailed change types
write_csv(change_detail, paste0(output_prefix, "_change_details.csv"))
cat("  Saved:", paste0(output_prefix, "_change_details.csv\n"))

# Export major class summary
write_csv(major_class_summary, paste0(output_prefix, "_major_class_summary.csv"))
cat("  Saved:", paste0(output_prefix, "_major_class_summary.csv\n"))

# Export overall summary
write_csv(overall_summary, paste0(output_prefix, "_overall_summary.csv"))
cat("  Saved:", paste0(output_prefix, "_overall_summary.csv\n"))

# Export top transitions
write_csv(top_transitions, paste0(output_prefix, "_top_transitions.csv"))
cat("  Saved:", paste0(output_prefix, "_top_transitions.csv\n\n"))

# ----------------------------------------------------------------------------
# ANALYSIS COMPLETE
# ----------------------------------------------------------------------------
cat(strrep("=", 78), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 78), "\n\n")

cat("Output files created:\n")
cat("  CSV files:\n")
cat("    -", paste0(output_prefix, "_full_results.csv\n"))
cat("    -", paste0(output_prefix, "_population_summary.csv\n"))
cat("    -", paste0(output_prefix, "_change_details.csv\n"))
cat("    -", paste0(output_prefix, "_major_class_summary.csv\n"))
cat("    -", paste0(output_prefix, "_overall_summary.csv\n"))
cat("    -", paste0(output_prefix, "_top_transitions.csv\n"))
cat("\n")
cat("  Spatial Maps:\n")
cat("    -", paste0(output_prefix, "_map_A_historical.png\n"))
cat("    -", paste0(output_prefix, "_map_B_future.png\n"))
cat("    -", paste0(output_prefix, "_map_change.png\n"))
cat("    -", paste0(output_prefix, "_map_AB_combined.png\n"))
cat("\n")
cat("  Statistical Plots:\n")
cat("    -", paste0(output_prefix, "_by_population.png\n"))
cat("    -", paste0(output_prefix, "_stacked_bar.png\n"))
cat("    -", paste0(output_prefix, "_top_transitions.png\n"))
cat("    -", paste0(output_prefix, "_major_classes.png\n"))
cat("    -", paste0(output_prefix, "_minor_transitions_top20.png\n"))
cat("    -", paste0(output_prefix, "_minor_transitions_by_pop.png\n"))
if(length(unique_hist_classes) <= 15 && length(unique_fut_classes) <= 15) {
  cat("    -", paste0(output_prefix, "_transition_heatmap.png\n"))
}
cat("    -", paste0(output_prefix, "_class_distribution.png\n"))
cat("    -", paste0(output_prefix, "_transition_flow.png\n"))
cat("\n")

cat("Summary:\n")
cat("  ", overall_summary$samples_changed, "out of", overall_summary$total_samples, 
    "samples (", overall_summary$pct_changed, "%) changed climate class\n")
cat("  ", overall_summary$major_class_changed, "samples (", 
    overall_summary$pct_major_changed, "%) changed major climate class\n\n")
