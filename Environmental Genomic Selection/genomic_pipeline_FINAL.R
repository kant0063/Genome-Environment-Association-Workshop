################################################################################
# GENOMIC SELECTION PIPELINE - COMPLETE STANDALONE VERSION
# 
# This single script includes:
# - Alternative VCF conversion (transpose first, then convert)
# - Automatic sample matching
# - Parallel processing
# - All optimizations
#
# Version: 1.3 (Final - December 2024)
################################################################################

# Required Libraries
suppressPackageStartupMessages({
  library(rrBLUP)
  library(SNPRelate)
  library(dplyr)
  library(ggplot2)
  library(raster)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(reshape2)
  
  # For heatmap with multiple fill scales
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    cat("Note: ggnewscale not installed. Installing now...\n")
    install.packages("ggnewscale", repos = "https://cloud.r-project.org")
  }
  library(ggnewscale)
  
  # Optional: BGLR for BayesCpi
  if (!requireNamespace("BGLR", quietly = TRUE)) {
    cat("Note: BGLR not installed. BayesCpi method will not be available.\n")
    cat("Install with: install.packages('BGLR')\n")
  } else {
    library(BGLR)
  }
})

################################################################################
# VCF CONVERSION - ALTERNATIVE METHOD (Transpose First)
################################################################################

vcf_to_rrblup_fast <- function(vcf_file, output_file = NULL, 
                               sample_subset = NULL, n_threads = 1) {
  
  cat("Converting VCF to rrBLUP format (alternative method)...\n")
  
  library(SNPRelate)
  
  # Convert to GDS
  gds_file <- tempfile(fileext = ".gds")
  cat("  Step 1: VCF to GDS...\n")
  snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only", verbose = FALSE)
  
  # Open and get sample list
  genofile <- snpgdsOpen(gds_file)
  all_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
  
  # Handle sample subsetting
  if (!is.null(sample_subset)) {
    use_samples <- intersect(sample_subset, all_samples)
    missing <- setdiff(sample_subset, all_samples)
    if (length(missing) > 0) {
      cat("    Note:", length(missing), "requested samples not in VCF (skipping)\n")
    }
    if (length(use_samples) == 0) {
      snpgdsClose(genofile)
      unlink(gds_file)
      stop("ERROR: None of the requested samples found in VCF!")
    }
  } else {
    use_samples <- all_samples
  }
  
  cat("  Step 2: Extracting", length(use_samples), "samples...\n")
  
  # Extract genotypes
  geno <- snpgdsGetGeno(genofile, sample.id = use_samples, 
                        verbose = FALSE, with.id = TRUE)
  
  # Get SNP info
  snp_id <- read.gdsn(index.gdsn(genofile, "snp.id"))
  snp_chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
  snp_pos <- read.gdsn(index.gdsn(genofile, "snp.position"))
  snp_allele <- read.gdsn(index.gdsn(genofile, "snp.allele"))
  
  snpgdsClose(genofile)
  
  cat("  Step 3: Processing genotype matrix...\n")
  
  # geno$genotype is SNPs x Samples from snpgdsGetGeno
  n_snps <- length(geno$snp.id)
  n_samples <- length(geno$sample.id)
  
  cat("    Input:", n_snps, "SNPs x", n_samples, "samples\n")
  
  # Step 3a: Transpose first to get Samples x SNPs
  cat("    Transposing...\n")
  gt_t <- t(geno$genotype)
  
  cat("    After transpose:", nrow(gt_t), "x", ncol(gt_t), "\n")
  
  # Step 3b: Create a completely NEW matrix with explicit dimensions
  cat("    Creating new matrix for conversion...\n")
  result <- matrix(as.numeric(gt_t), nrow = n_samples, ncol = n_snps)
  
  cat("    New matrix dimensions:", nrow(result), "x", ncol(result), "\n")
  
  # Step 3c: Convert encoding on the new matrix
  cat("    Converting encoding...\n")
  result[result == 0] <- -1  # 0/0 -> -1
  result[result == 2] <- 1   # 1/1 -> 1
  # 1 -> 0 (already correct for het)
  # NA -> NA (already correct)
  
  cat("  Step 4: Adding names...\n")
  cat("    Result dimensions:", nrow(result), "rows x", ncol(result), "cols\n")
  cat("    Sample IDs to assign:", length(geno$sample.id), "\n")
  cat("    SNP IDs to assign:", length(geno$snp.id), "\n")
  
  # Verify dimensions
  if (nrow(result) != n_samples || ncol(result) != n_snps) {
    cat("    ERROR: Dimension mismatch!\n")
    cat("      Expected:", n_samples, "x", n_snps, "\n")
    cat("      Got:", nrow(result), "x", ncol(result), "\n")
    stop("Matrix dimensions don't match expected size")
  }
  
  # Assign names
  rownames(result) <- geno$sample.id
  colnames(result) <- geno$snp.id
  
  cat("    ✓ Names assigned successfully\n")
  
  # Write hapmap file if requested
  if (!is.null(output_file)) {
    cat("  Step 5: Writing hapmap file...\n")
    
    chrom_unique <- unique(snp_chr)
    chrom_num <- match(snp_chr, chrom_unique)
    
    # For hapmap format, we need SNPs as rows, samples as columns
    hapmap <- data.frame(
      `rs#` = snp_id,
      allele = snp_allele,
      chrom = chrom_num,
      pos = snp_pos,
      t(result),  # Transpose back: Samples x SNPs -> SNPs x Samples
      check.names = FALSE
    )
    
    write.table(hapmap, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("    ✓ Hapmap written to:", output_file, "\n")
  }
  
  # Clean up
  unlink(gds_file)
  
  cat("  ✓ Conversion complete:", nrow(result), "samples x", ncol(result), "SNPs\n\n")
  
  return(result)
}

################################################################################
# CORE SELECTION - OPTIMIZED
################################################################################

greedy_core_selection_fast <- function(geno_matrix, n_core, method = "IBS",
                                       max_snps = 10000, n_threads = 1, seed = 123) {
  
  cat("=== Core Selection ===\n")
  
  # Check for missing data
  n_missing <- sum(is.na(geno_matrix))
  if (n_missing > 0) {
    cat("  Warning: Genotype matrix contains", n_missing, "missing values\n")
    cat("  Note: Using mean imputation for distance calculation\n")
    # Simple mean imputation for distance calculation only
    geno_clean <- geno_matrix
    for (j in 1:ncol(geno_clean)) {
      if (any(is.na(geno_clean[, j]))) {
        geno_clean[is.na(geno_clean[, j]), j] <- mean(geno_clean[, j], na.rm = TRUE)
      }
    }
  } else {
    geno_clean <- geno_matrix
  }
  
  cat("  Using", ncol(geno_clean), "SNPs for core selection\n")
  
  # Sample SNPs if too many
  if (!is.null(max_snps) && ncol(geno_clean) > max_snps) {
    cat("  Sampling", max_snps, "SNPs for distance calculation...\n")
    set.seed(seed)
    geno_dist <- geno_clean[, sample(ncol(geno_clean), max_snps)]
  } else {
    geno_dist <- geno_clean
  }
  
  # Calculate distance matrix
  cat("  Computing genetic distances (", method, ")...\n", sep = "")
  
  if (method == "IBS") {
    # Use SNPRelate for fast IBS calculation
    gds_file <- tempfile(fileext = ".gds")
    snpgdsCreateGeno(gds_file, 
                     genmat = t(geno_dist),
                     sample.id = rownames(geno_dist),
                     snp.id = colnames(geno_dist),
                     snpfirstdim = TRUE)
    genofile <- snpgdsOpen(gds_file)
    ibs <- snpgdsIBS(genofile, num.thread = n_threads, verbose = FALSE)
    snpgdsClose(genofile)
    unlink(gds_file)
    dist_matrix <- 1 - ibs$ibs
  } else {
    dist_matrix <- as.matrix(dist(geno_dist, method = method))
  }
  
  # Greedy selection
  cat("  Performing greedy selection...\n")
  set.seed(seed)
  selected_idx <- sample(nrow(dist_matrix), 1)
  
  pb <- txtProgressBar(min = 0, max = n_core, style = 3)
  while (length(selected_idx) < n_core) {
    remaining_idx <- setdiff(1:nrow(dist_matrix), selected_idx)
    min_distances <- apply(dist_matrix[remaining_idx, selected_idx, drop = FALSE], 1, min)
    next_idx <- remaining_idx[which.max(min_distances)]
    selected_idx <- c(selected_idx, next_idx)
    setTxtProgressBar(pb, length(selected_idx))
  }
  close(pb)
  
  core_ids <- rownames(geno_clean)[selected_idx]
  cat("\n  ✓ Selected", n_core, "core samples\n\n")
  
  return(core_ids)
}

################################################################################
# ENVIRONMENTAL DATA EXTRACTION
################################################################################

extract_worldclim_data <- function(coords, worldclim_path, variables = "bio", 
                                   use_all_tifs = FALSE, output_file = NULL) {
  
  cat("=== Extracting Environmental Data ===\n")
  cat("  Looking in:", worldclim_path, "\n")
  
  # Check if path exists
  if (!dir.exists(worldclim_path)) {
    # Try making it absolute
    worldclim_path <- normalizePath(worldclim_path, mustWork = FALSE)
    cat("  Trying absolute path:", worldclim_path, "\n")
    
    if (!dir.exists(worldclim_path)) {
      stop("Directory not found: ", worldclim_path)
    }
  }
  
  # Find ALL GeoTIFF files in the directory
  cat("  Scanning for GeoTIFF files...\n")
  all_tif_files <- list.files(worldclim_path, pattern = "\\.tif$|\\.tiff$", 
                               full.names = TRUE, ignore.case = TRUE)
  
  if (length(all_tif_files) == 0) {
    stop("No .tif or .tiff files found in: ", worldclim_path)
  }
  
  cat("  Found", length(all_tif_files), "GeoTIFF files\n")
  
  # Decide which files to use
  if (use_all_tifs) {
    cat("  Using ALL", length(all_tif_files), "GeoTIFF files\n")
    files_to_use <- all_tif_files
    
    cat("  Files:\n")
    for (i in 1:min(10, length(files_to_use))) {
      cat("    ", basename(files_to_use[i]), "\n")
    }
    if (length(files_to_use) > 10) {
      cat("    ... and", length(files_to_use) - 10, "more\n")
    }
    
  } else {
    # Try to identify bioclimatic variables (bio_1 through bio_19)
    bio_pattern <- "bio[_-]?(\\d+)\\.tif"
    bio_files <- all_tif_files[grepl(bio_pattern, all_tif_files, ignore.case = TRUE)]
    
    if (length(bio_files) > 0) {
      cat("  Identified", length(bio_files), "bioclimatic variable files\n")
      
      # Sort them by number to ensure correct order
      bio_numbers <- as.numeric(gsub(".*bio[_-]?(\\d+)\\.tif.*", "\\1", 
                                     bio_files, ignore.case = TRUE))
      bio_files <- bio_files[order(bio_numbers)]
      
      cat("  Using bioclimatic files:\n")
      for (i in 1:min(5, length(bio_files))) {
        cat("    ", basename(bio_files[i]), "\n")
      }
      if (length(bio_files) > 5) {
        cat("    ... and", length(bio_files) - 5, "more\n")
      }
      
      # Show what's NOT being used
      other_files <- setdiff(all_tif_files, bio_files)
      if (length(other_files) > 0) {
        cat("  Note:", length(other_files), "other .tif files found but not used\n")
        cat("        Set use_all_tifs=TRUE to include all files\n")
      }
      
      files_to_use <- bio_files
      
    } else {
      # No bioclimatic pattern found, use ALL tif files
      cat("  No bioclimatic pattern detected, using ALL", length(all_tif_files), "files\n")
      
      cat("  Files:\n")
      for (i in 1:min(10, length(all_tif_files))) {
        cat("    ", basename(all_tif_files[i]), "\n")
      }
      if (length(all_tif_files) > 10) {
        cat("    ... and", length(all_tif_files) - 10, "more\n")
      }
      
      files_to_use <- all_tif_files
    }
  }
  
  # Load the rasters
  cat("  Loading raster stack...\n")
  tryCatch({
    climate_layers <- stack(files_to_use)
    cat("  ✓ Loaded", nlayers(climate_layers), "layers\n")
  }, error = function(e) {
    stop("Error loading rasters: ", e$message)
  })
  
  # Extract values
  cat("  Extracting data for", nrow(coords), "locations...\n")
  coords_sp <- coords[, c("Longitude", "Latitude")]
  
  tryCatch({
    climate_values <- extract(climate_layers, coords_sp)
    cat("  ✓ Extracted", ncol(climate_values), "variables\n")
    
    # Check for missing values
    missing_per_sample <- rowSums(is.na(climate_values))
    if (any(missing_per_sample > 0)) {
      cat("  Warning:", sum(missing_per_sample > 0), "samples have missing values\n")
      cat("          This may be due to coordinates outside raster extent\n")
    }
    
  }, error = function(e) {
    stop("Error extracting data: ", e$message, 
         "\nCheck that coordinates are in the same CRS as rasters (usually WGS84)")
  })
  
  # Combine with coordinates
  result <- cbind(coords, climate_values)
  
  # Save if requested
  if (!is.null(output_file)) {
    write.csv(result, output_file, row.names = FALSE)
    cat("  ✓ Saved to:", output_file, "\n")
  }
  
  cat("  ✓ Extraction complete\n\n")
  return(result)
}

################################################################################
# PARALLEL CROSS-VALIDATION
################################################################################

kfold_cross_validation_parallel <- function(geno, pheno, train_set, 
                                            k_fold = 10, n_reps = 25,
                                            method = "rrBLUP", kernel = NULL,
                                            n_cores = NULL) {
  
  cat("=== Cross-Validation (Parallel) ===\n")
  
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  cat("  Using", n_cores, "CPU cores\n")
  cat("  Method:", method, "\n")
  
  # Pre-compute kernels if needed
  kernel_gauss <- NULL
  kernel_exp <- NULL
  
  if (method == "GBLUP-Gaussian") {
    cat("  Computing Gaussian kernel...\n")
    geno_scaled <- scale(geno, center = TRUE, scale = TRUE)
    geno_scaled[is.na(geno_scaled)] <- 0
    D_all <- as.matrix(dist(geno_scaled))
    h <- median(D_all[D_all > 0])
    if (is.na(h) || h == 0 || is.infinite(h)) h <- 1
    kernel_gauss <- exp(-D_all^2 / (2 * h^2))
    diag(kernel_gauss) <- diag(kernel_gauss) + 1e-6
  } else if (method == "GBLUP-Exponential") {
    cat("  Computing Exponential kernel...\n")
    geno_scaled <- scale(geno, center = TRUE, scale = TRUE)
    geno_scaled[is.na(geno_scaled)] <- 0
    D_all <- as.matrix(dist(geno_scaled))
    h <- median(D_all[D_all > 0])
    if (is.na(h) || h == 0 || is.infinite(h)) h <- 1
    kernel_exp <- exp(-D_all / h)
    diag(kernel_exp) <- diag(kernel_exp) + 1e-6
  } else if (method == "GBLUP") {
    cat("  Computing relationship matrix...\n")
    kernel <- A.mat(geno)
  } else if (method == "BayesCpi") {
    cat("  Note: BayesCpi will be fitted separately for each fold\n")
  }
  
  n_traits <- ncol(pheno)
  results <- data.frame(
    trait = character(n_traits),
    method = character(n_traits),
    mean_accuracy = numeric(n_traits),
    sd_accuracy = numeric(n_traits),
    stringsAsFactors = FALSE
  )
  
  # Setup parallel
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # Export variables based on method
  if (method == "GBLUP-Gaussian") {
    clusterExport(cl, c("geno", "pheno", "train_set", "method", "kernel_gauss", "k_fold"), 
                  envir = environment())
  } else if (method == "GBLUP-Exponential") {
    clusterExport(cl, c("geno", "pheno", "train_set", "method", "kernel_exp", "k_fold"), 
                  envir = environment())
  } else if (method == "GBLUP") {
    clusterExport(cl, c("geno", "pheno", "train_set", "method", "kernel", "k_fold"), 
                  envir = environment())
  } else {
    clusterExport(cl, c("geno", "pheno", "train_set", "method", "k_fold"), 
                  envir = environment())
  }
  
  # Load packages on cluster
  if (method == "BayesCpi") {
    clusterEvalQ(cl, {library(rrBLUP); library(BGLR)})
  } else {
    clusterEvalQ(cl, library(rrBLUP))
  }
  
  for (t in 1:n_traits) {
    trait_name <- colnames(pheno)[t]
    cat("  Trait:", trait_name, "\n")
    
    # Get training samples with complete data for this trait
    train_complete <- train_set[!is.na(train_set[, t]), , drop = FALSE]
    
    cat("    Training samples with data:", nrow(train_complete), "\n")
    
    if (nrow(train_complete) < k_fold) {
      cat("    Skipping - not enough samples (need at least", k_fold, ")\n")
      next
    }
    
    # Make sure we have sample IDs
    train_ids_complete <- rownames(train_complete)
    
    if (is.null(train_ids_complete) || any(is.na(train_ids_complete))) {
      cat("    ERROR: Training set missing rownames\n")
      next
    }
    
    # Verify all training IDs are in geno and pheno
    missing_in_geno <- setdiff(train_ids_complete, rownames(geno))
    missing_in_pheno <- setdiff(train_ids_complete, rownames(pheno))
    
    if (length(missing_in_geno) > 0) {
      cat("    ERROR:", length(missing_in_geno), "training samples not in genotype matrix\n")
      next
    }
    if (length(missing_in_pheno) > 0) {
      cat("    ERROR:", length(missing_in_pheno), "training samples not in phenotype matrix\n")
      next
    }
    
    # Parallel CV with error handling
    tryCatch({
      # Collect errors for diagnosis
      all_errors <- c()
      
      accuracies <- foreach(rep = 1:n_reps, .combine = c, .packages = "rrBLUP",
                           .errorhandling = "pass") %dopar% {
        
        tryCatch({
          # Create folds
          n_samples <- length(train_ids_complete)
          fold_assignments <- cut(sample(n_samples), breaks = k_fold, labels = FALSE)
          
          fold_predictions <- list()
          
          for (fold_num in 1:k_fold) {
            # Get validation and training IDs for this fold
            val_indices <- which(fold_assignments == fold_num)
            
            if (length(val_indices) == 0) next
            
            val_ids <- train_ids_complete[val_indices]
            train_ids_fold <- train_ids_complete[-val_indices]
            
            if (length(train_ids_fold) < 2) next  # Need at least 2 training samples
            
            # Get training phenotypes
            y_train <- pheno[train_ids_fold, trait_name]
            
            # Train model and predict
            if (method == "rrBLUP") {
              solve_out <- mixed.solve(y = y_train, Z = geno[train_ids_fold, ], 
                                       SE = FALSE, return.Hinv = FALSE)
              pred_matrix <- geno[val_ids, , drop = FALSE] %*% solve_out$u
              predictions <- as.vector(pred_matrix)
              names(predictions) <- val_ids
              
            } else if (method == "GBLUP-Gaussian") {
              K_subset <- kernel_gauss[c(train_ids_fold, val_ids), c(train_ids_fold, val_ids)]
              y_all <- rep(NA, length(c(train_ids_fold, val_ids)))
              names(y_all) <- c(train_ids_fold, val_ids)
              y_all[train_ids_fold] <- y_train
              pheno_df <- data.frame(id = names(y_all), trait = y_all)
              solve_out <- kin.blup(data = pheno_df, geno = "id", pheno = "trait", K = K_subset)
              predictions <- solve_out$g[val_ids]
              
            } else if (method == "GBLUP-Exponential") {
              K_subset <- kernel_exp[c(train_ids_fold, val_ids), c(train_ids_fold, val_ids)]
              y_all <- rep(NA, length(c(train_ids_fold, val_ids)))
              names(y_all) <- c(train_ids_fold, val_ids)
              y_all[train_ids_fold] <- y_train
              pheno_df <- data.frame(id = names(y_all), trait = y_all)
              solve_out <- kin.blup(data = pheno_df, geno = "id", pheno = "trait", K = K_subset)
              predictions <- solve_out$g[val_ids]
              
            } else if (method == "GBLUP") {
              K_subset <- kernel[c(train_ids_fold, val_ids), c(train_ids_fold, val_ids)]
              y_all <- rep(NA, length(c(train_ids_fold, val_ids)))
              names(y_all) <- c(train_ids_fold, val_ids)
              y_all[train_ids_fold] <- y_train
              pheno_df <- data.frame(id = names(y_all), trait = y_all)
              solve_out <- kin.blup(data = pheno_df, geno = "id", pheno = "trait", K = K_subset)
              predictions <- solve_out$g[val_ids]
              
            } else if (method == "BayesCpi") {
              # BayesCπ using BGLR
              fit <- BGLR(y = y_train, ETA = list(list(X = geno[train_ids_fold, ], model = "BayesC")),
                         nIter = 5000, burnIn = 1000, verbose = FALSE)
              pred_matrix <- geno[val_ids, , drop = FALSE] %*% fit$ETA[[1]]$b
              predictions <- as.vector(pred_matrix)
              names(predictions) <- val_ids
            }
            
            fold_predictions[[fold_num]] <- predictions
          }
          
          # Combine predictions
          all_gebv <- unlist(fold_predictions)
          
          # Calculate correlation
          if (length(all_gebv) > 5 && !any(is.na(all_gebv))) {
            observed <- pheno[names(all_gebv), trait_name]
            if (!any(is.na(observed))) {
              cor(all_gebv, observed, use = "complete.obs")
            } else {
              NA
            }
          } else {
            NA
          }
        }, error = function(e) {
          # Return error message instead of NA
          return(list(error = e$message))
        })
      }
      
      # Check for errors
      error_results <- sapply(accuracies, function(x) is.list(x) && !is.null(x$error))
      if (any(error_results)) {
        # Show first few unique errors
        errors <- unique(sapply(accuracies[error_results], function(x) x$error))
        cat("    CV ERRORS detected:\n")
        for (i in 1:min(3, length(errors))) {
          cat("      -", errors[i], "\n")
        }
      }
      
      # Keep only numeric results
      accuracies <- unlist(accuracies[!error_results])
      accuracies <- accuracies[!is.na(accuracies)]
      
      if (length(accuracies) > 0) {
        results[t, ] <- list(trait_name, method, mean(accuracies), sd(accuracies))
        cat("    Accuracy:", round(mean(accuracies), 3), 
            "±", round(sd(accuracies), 3), 
            "(", length(accuracies), "/", n_reps, "successful reps )\n")
      } else {
        cat("    ERROR: All", n_reps, "CV reps failed\n")
      }
      
    }, error = function(e) {
      cat("    OUTER ERROR in CV:", e$message, "\n")
    })
  }
  
  stopCluster(cl)
  
  cat("  ✓ Cross-validation complete\n\n")
  return(list(results = na.omit(results), k_fold = k_fold, 
              n_reps = n_reps, method = method))
}

################################################################################
# GENOMIC SELECTION
################################################################################

genomic_selection <- function(geno_train, geno_pred, pheno_train, method = "rrBLUP") {
  
  cat("=== Genomic Selection ===\n")
  
  # Align SNPs
  common_snps <- intersect(colnames(geno_train), colnames(geno_pred))
  geno_train <- geno_train[, common_snps]
  geno_pred <- geno_pred[, common_snps]
  
  cat("  Training samples:", nrow(geno_train), "\n")
  cat("  Prediction samples:", nrow(geno_pred), "\n")
  cat("  Using", length(common_snps), "common SNPs\n")
  cat("  Method:", method, "\n")
  
  n_traits <- ncol(pheno_train)
  
  # Calculate GEBVs for all prediction samples
  gebv_pred <- matrix(NA, nrow = nrow(geno_pred), ncol = n_traits)
  rownames(gebv_pred) <- rownames(geno_pred)
  colnames(gebv_pred) <- colnames(pheno_train)
  
  for (t in 1:n_traits) {
    if (method == "rrBLUP") {
      solve_out <- mixed.solve(y = pheno_train[, t], Z = geno_train, 
                              SE = FALSE, return.Hinv = FALSE)
      gebv_pred[, t] <- geno_pred %*% solve_out$u
      
    } else if (method == "GBLUP-Gaussian") {
      # Gaussian kernel: K = exp(-D^2 / (2*h^2))
      # Scale genotypes for better distance calculation
      geno_train_scaled <- scale(geno_train, center = TRUE, scale = TRUE)
      geno_pred_scaled <- scale(geno_pred, center = TRUE, scale = TRUE)
      
      # Replace any NAs from scaling with 0
      geno_train_scaled[is.na(geno_train_scaled)] <- 0
      geno_pred_scaled[is.na(geno_pred_scaled)] <- 0
      
      # Compute Euclidean distance
      D_train <- as.matrix(dist(geno_train_scaled))
      h <- median(D_train[D_train > 0])  # Exclude zeros
      
      if (is.na(h) || h == 0 || is.infinite(h)) {
        cat("    Warning: Invalid bandwidth, using h=1\n")
        h <- 1
      }
      
      K_train <- exp(-D_train^2 / (2 * h^2))
      # Add small constant to diagonal for numerical stability
      diag(K_train) <- diag(K_train) + 1e-6
      
      # Compute distances between prediction and training
      all_scaled <- rbind(geno_pred_scaled, geno_train_scaled)
      D_all <- as.matrix(dist(all_scaled))
      D_pred_train <- D_all[1:nrow(geno_pred), 
                            (nrow(geno_pred)+1):nrow(all_scaled)]
      K_pred_train <- exp(-D_pred_train^2 / (2 * h^2))
      
      # Fit kinship BLUP on training
      y_train <- pheno_train[, t]
      names(y_train) <- rownames(geno_train)
      pheno_df <- data.frame(id = names(y_train), trait = y_train)
      
      tryCatch({
        solve_out <- kin.blup(data = pheno_df, geno = "id", pheno = "trait", 
                             K = K_train, GAUSS = FALSE)
        # Predict using kernel with regularized inverse
        K_train_reg <- K_train + diag(1e-4, nrow(K_train))
        gebv_pred[, t] <- K_pred_train %*% solve(K_train_reg) %*% solve_out$g
      }, error = function(e) {
        cat("    Warning: Gaussian kernel failed for trait", colnames(pheno_train)[t], "\n")
        gebv_pred[, t] <- mean(y_train, na.rm = TRUE)
      })
      
    } else if (method == "GBLUP-Exponential") {
      # Exponential kernel: K = exp(-D / h)
      geno_train_scaled <- scale(geno_train, center = TRUE, scale = TRUE)
      geno_pred_scaled <- scale(geno_pred, center = TRUE, scale = TRUE)
      
      # Replace any NAs from scaling with 0
      geno_train_scaled[is.na(geno_train_scaled)] <- 0
      geno_pred_scaled[is.na(geno_pred_scaled)] <- 0
      
      D_train <- as.matrix(dist(geno_train_scaled))
      h <- median(D_train[D_train > 0])
      
      if (is.na(h) || h == 0 || is.infinite(h)) {
        cat("    Warning: Invalid bandwidth, using h=1\n")
        h <- 1
      }
      
      K_train <- exp(-D_train / h)
      diag(K_train) <- diag(K_train) + 1e-6
      
      all_scaled <- rbind(geno_pred_scaled, geno_train_scaled)
      D_all <- as.matrix(dist(all_scaled))
      D_pred_train <- D_all[1:nrow(geno_pred), 
                            (nrow(geno_pred)+1):nrow(all_scaled)]
      K_pred_train <- exp(-D_pred_train / h)
      
      y_train <- pheno_train[, t]
      names(y_train) <- rownames(geno_train)
      pheno_df <- data.frame(id = names(y_train), trait = y_train)
      
      tryCatch({
        solve_out <- kin.blup(data = pheno_df, geno = "id", pheno = "trait", 
                             K = K_train, GAUSS = FALSE)
        K_train_reg <- K_train + diag(1e-4, nrow(K_train))
        gebv_pred[, t] <- K_pred_train %*% solve(K_train_reg) %*% solve_out$g
      }, error = function(e) {
        cat("    Warning: Exponential kernel failed for trait", colnames(pheno_train)[t], "\n")
        gebv_pred[, t] <- mean(y_train, na.rm = TRUE)
      })
      
    } else if (method == "BayesCpi") {
      # BayesCπ using BGLR
      y <- pheno_train[, t]
      
      # Run BGLR
      fit <- BGLR(y = y, ETA = list(list(X = geno_train, model = "BayesC")),
                  nIter = 12000, burnIn = 2000, verbose = FALSE)
      
      # Predict
      gebv_pred[, t] <- geno_pred %*% fit$ETA[[1]]$b
      
    } else {
      stop("Unknown method: ", method)
    }
  }
  
  # Separate training samples if they're in prediction set
  train_ids <- rownames(geno_train)
  is_training <- rownames(gebv_pred) %in% train_ids
  
  gebv_train <- gebv_pred[is_training, , drop = FALSE]
  
  cat("  ✓ GEBVs calculated for", n_traits, "traits\n")
  cat("  ✓ Training samples in results:", sum(is_training), "\n")
  cat("  ✓ Novel samples in results:", sum(!is_training), "\n\n")
  
  return(list(
    gebv_all = gebv_pred,  # All samples
    gebv_train = gebv_train,  # Just training samples
    gebv_pred = gebv_pred  # All samples (same as gebv_all)
  ))
}

################################################################################
# VISUALIZATION
################################################################################

plot_cv_results <- function(cv_results) {
  all_results <- do.call(rbind, lapply(cv_results, function(x) x$results))
  all_results <- na.omit(all_results)
  
  if (nrow(all_results) == 0) {
    cat("  Warning: No CV results to plot (all NA)\n")
    # Create empty plot
    p <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, 
               label = "No cross-validation results\n(all accuracies were NA)", 
               size = 6) +
      theme_void()
    return(p)
  }
  
  # Define colors for methods
  method_colors <- c(
    "rrBLUP" = "#E41A1C",
    "GBLUP" = "#377EB8",
    "GBLUP-Gaussian" = "#4DAF4A",
    "GBLUP-Exponential" = "#984EA3",
    "BayesCpi" = "#FF7F00"
  )
  
  ggplot(all_results, aes(x = trait, y = mean_accuracy, color = method, group = method)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = mean_accuracy - sd_accuracy,
                     ymax = mean_accuracy + sd_accuracy), 
                  width = 0.3, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = method_colors, name = "Method") +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "top") +
    labs(title = "Cross-Validation Prediction Accuracy", 
         x = "Environmental Trait", 
         y = "Prediction Accuracy (r)") +
    ylim(-0.2, 1)
}

################################################################################
# GEBV HEATMAP WITH COUNTRY CLUSTERING
################################################################################

assign_country <- function(lat, lon) {
  # Simple country assignment based on lat/lon
  # You can improve this with rnaturalearth or other packages
  
  if (is.na(lat) || is.na(lon)) return("Unknown")
  
  # Europe/Middle East
  if (lon >= -10 && lon <= 50 && lat >= 35 && lat <= 70) {
    if (lon >= 25 && lon <= 45 && lat >= 35 && lat <= 42) return("Turkey")
    if (lon >= 6 && lon <= 19 && lat >= 36 && lat <= 48) return("Italy")
    if (lon >= -10 && lon <= 5 && lat >= 36 && lat <= 44) return("Spain")
    if (lon >= 20 && lon <= 30 && lat >= 35 && lat <= 42) return("Greece")
    if (lon >= -5 && lon <= 10 && lat >= 42 && lat <= 52) return("France")
    return("Europe-Other")
  }
  
  # North Africa/Middle East
  if (lon >= -10 && lon <= 60 && lat >= 15 && lat <= 35) {
    if (lon >= 25 && lon <= 37 && lat >= 22 && lat <= 32) return("Egypt")
    if (lon >= 35 && lon <= 48 && lat >= 29 && lat <= 38) return("Syria-Iraq")
    if (lon >= 44 && lon <= 64 && lat >= 24 && lat <= 32) return("Iran")
    return("North-Africa-Middle-East")
  }
  
  # South Asia
  if (lon >= 60 && lon <= 95 && lat >= 5 && lat <= 40) {
    if (lon >= 68 && lon <= 80 && lat >= 23 && lat <= 37) return("Pakistan")
    if (lon >= 68 && lon <= 98 && lat >= 6 && lat <= 36) return("India")
    if (lon >= 80 && lon <= 90 && lat >= 20 && lat <= 28) return("Bangladesh")
    return("South-Asia")
  }
  
  # East Asia
  if (lon >= 95 && lon <= 145 && lat >= 20 && lat <= 50) {
    if (lon >= 100 && lon <= 125 && lat >= 18 && lat <= 54) return("China")
    return("East-Asia")
  }
  
  # Americas
  if (lon >= -125 && lon <= -30) {
    if (lat >= 25 && lat <= 50) return("USA-Canada")
    if (lat >= 15 && lat <= 35) return("Mexico-Central-America")
    if (lat >= -55 && lat <= 15) return("South-America")
    return("Americas")
  }
  
  # Australia/Oceania
  if (lon >= 110 && lon <= 180 && lat >= -45 && lat <= -10) return("Australia")
  
  # Africa
  if (lon >= -20 && lon <= 55 && lat >= -35 && lat <= 15) return("Sub-Saharan-Africa")
  
  return("Other")
}

assign_country_specific <- function(lat, lon) {
  # Comprehensive country assignment with global coverage
  # Check most specific boundaries first
  
  if (is.na(lat) || is.na(lon)) return("Unknown")
  
  # ============================================================================
  # AFRICA - Comprehensive coverage
  # ============================================================================
  
  # NORTH AFRICA
  if (lon >= -18 && lon <= -1 && lat >= 21 && lat <= 37) return("Morocco")
  if (lon >= -9 && lon <= -5 && lat >= 27 && lat <= 29) return("Western-Sahara")
  if (lon >= -2 && lon <= 12 && lat >= 18 && lat <= 38) return("Algeria")
  if (lon >= 7 && lon <= 12 && lat >= 30 && lat <= 38) return("Tunisia")
  if (lon >= 9 && lon <= 26 && lat >= 19 && lat <= 34) return("Libya")
  if (lon >= 25 && lon <= 37 && lat >= 22 && lat <= 32) return("Egypt")
  if (lon >= 36 && lon <= 44 && lat >= 3 && lat <= 19) return("Sudan")
  if (lon >= 32 && lon <= 49 && lat >= 3 && lat <= 18) return("South-Sudan-Ethiopia-Eritrea")
  
  # WEST AFRICA (Atlantic coast)
  if (lon >= -18 && lon <= -11 && lat >= 12 && lat <= 17) return("Mauritania")
  if (lon >= -17 && lon <= -11 && lat >= 12 && lat <= 17) return("Senegal-Gambia")
  if (lon >= -16 && lon <= -10 && lat >= 10 && lat <= 13) return("Guinea-Bissau-Guinea")
  if (lon >= -14 && lon <= -7 && lat >= 7 && lat <= 11) return("Sierra-Leone-Liberia")
  if (lon >= -9 && lon <= 3 && lat >= 4 && lat <= 12) return("Ivory-Coast-Ghana")
  if (lon >= 0 && lon <= 4 && lat >= 6 && lat <= 12) return("Togo-Benin")
  if (lon >= 2 && lon <= 15 && lat >= 4 && lat <= 14) return("Nigeria-Cameroon")
  
  # CENTRAL AFRICA
  if (lon >= 8 && lon <= 19 && lat >= -3 && lat <= 4) return("Gabon-Congo-Eq-Guinea")
  if (lon >= 11 && lon <= 32 && lat >= -13 && lat <= 6) return("DR-Congo-CAR")
  if (lon >= 28 && lon <= 35 && lat >= -4 && lat <= 2) return("Rwanda-Burundi-Uganda")
  
  # EAST AFRICA
  if (lon >= 33 && lon <= 42 && lat >= -12 && lat <= 0) return("Kenya-Tanzania")
  if (lon >= 40 && lon <= 52 && lat >= -1 && lat <= 12) return("Somalia-Djibouti")
  if (lon >= 29 && lon <= 41 && lat >= -26 && lat <= -10) return("Mozambique-Malawi-Zimbabwe")
  if (lon >= 25 && lon <= 34 && lat >= -23 && lat <= -15) return("Zambia")
  
  # SOUTHERN AFRICA
  if (lon >= 21 && lon <= 30 && lat >= -29 && lat <= -17) return("Botswana")
  if (lon >= 11 && lon <= 26 && lat >= -29 && lat <= -17) return("Namibia-Angola")
  if (lon >= 16 && lon <= 33 && lat >= -35 && lat <= -22) return("South-Africa-Lesotho-Swaziland")
  if (lon >= 43 && lon <= 51 && lat >= -26 && lat <= -12) return("Madagascar")
  
  # SAHEL (interior)
  if (lon >= -5 && lon <= 5 && lat >= 10 && lat <= 18) return("Mali-Burkina-Faso-Niger")
  if (lon >= 13 && lon <= 24 && lat >= 7 && lat <= 24) return("Chad")
  
  # ============================================================================
  # EUROPE
  # ============================================================================
  
  # Mediterranean
  if (lon >= 26 && lon <= 45 && lat >= 36 && lat <= 43) return("Turkey")
  if (lon >= 6 && lon <= 19 && lat >= 36 && lat <= 48) return("Italy")
  if (lon >= -10 && lon <= -6 && lat >= 37 && lat <= 43) return("Portugal")
  if (lon >= -10 && lon <= 5 && lat >= 36 && lat <= 44) return("Spain")
  if (lon >= 20 && lon <= 30 && lat >= 35 && lat <= 42) return("Greece")
  if (lon >= -5 && lon <= 10 && lat >= 42 && lat <= 52) return("France")
  if (lon >= 7 && lon <= 12 && lat >= 45 && lat <= 48) return("Switzerland")
  
  # Western/Central Europe
  if (lon >= 5 && lon <= 16 && lat >= 47 && lat <= 56) return("Germany")
  if (lon >= -11 && lon <= 2 && lat >= 50 && lat <= 61) return("UK-Ireland")
  if (lon >= 3 && lon <= 8 && lat >= 49 && lat <= 54) return("Belgium-Netherlands-Luxembourg")
  if (lon >= 9 && lon <= 13 && lat >= 54 && lat <= 59) return("Denmark")
  if (lon >= 10 && lon <= 25 && lat >= 55 && lat <= 72) return("Sweden-Norway")
  if (lon >= 20 && lon <= 32 && lat >= 59 && lat <= 71) return("Finland")
  if (lon >= 19 && lon <= 25 && lat >= 49 && lat <= 55) return("Poland")
  if (lon >= 9 && lon <= 18 && lat >= 46 && lat <= 49) return("Austria")
  if (lon >= 13 && lon <= 20 && lat >= 48 && lat <= 51) return("Czech-Slovakia")
  if (lon >= 16 && lon <= 23 && lat >= 45 && lat <= 49) return("Hungary")
  
  # Balkans
  if (lon >= 13 && lon <= 17 && lat >= 42 && lat <= 47) return("Slovenia-Croatia-Bosnia")
  if (lon >= 18 && lon <= 23 && lat >= 39 && lat <= 44) return("Albania-Macedonia-Montenegro")
  if (lon >= 19 && lon <= 24 && lat >= 41 && lat <= 47) return("Serbia-Kosovo")
  if (lon >= 22 && lon <= 30 && lat >= 41 && lat <= 45) return("Bulgaria")
  if (lon >= 20 && lon <= 30 && lat >= 44 && lat <= 49) return("Romania")
  
  # Eastern Europe
  if (lon >= 23 && lon <= 41 && lat >= 44 && lat <= 53) return("Ukraine")
  if (lon >= 23 && lon <= 33 && lat >= 51 && lat <= 57) return("Belarus")
  if (lon >= 21 && lon <= 29 && lat >= 53 && lat <= 60) return("Baltic-States")
  if (lon >= 19 && lon <= 25 && lat >= 53 && lat <= 56) return("Lithuania")
  if (lon >= 27 && lon <= 180 && lat >= 41 && lat <= 82) return("Russia")
  
  # ============================================================================
  # MIDDLE EAST
  # ============================================================================
  
  if (lon >= 35 && lon <= 43 && lat >= 33 && lat <= 38) return("Syria")
  if (lon >= 35 && lon <= 37 && lat >= 33 && lat <= 35) return("Lebanon")
  if (lon >= 35 && lon <= 40 && lat >= 29 && lat <= 34) return("Jordan")
  if (lon >= 34 && lon <= 36 && lat >= 31 && lat <= 33) return("Israel-Palestine")
  if (lon >= 43 && lon <= 49 && lat >= 29 && lat <= 38) return("Iraq")
  if (lon >= 44 && lon <= 64 && lat >= 25 && lat <= 41) return("Iran")
  if (lon >= 46 && lon <= 57 && lat >= 13 && lat <= 32) return("Saudi-Arabia")
  if (lon >= 51 && lon <= 57 && lat >= 22 && lat <= 27) return("UAE-Oman")
  if (lon >= 46 && lon <= 52 && lat >= 28 && lat <= 31) return("Kuwait")
  if (lon >= 50 && lon <= 52 && lat >= 25 && lat <= 27) return("Bahrain-Qatar")
  if (lon >= 42 && lon <= 55 && lat >= 12 && lat <= 19) return("Yemen")
  
  # ============================================================================
  # CENTRAL ASIA
  # ============================================================================
  
  if (lon >= 46 && lon <= 56 && lat >= 36 && lat <= 48) return("Turkmenistan")
  if (lon >= 52 && lon <= 70 && lat >= 36 && lat <= 46) return("Uzbekistan")
  if (lon >= 66 && lon <= 81 && lat >= 36 && lat <= 44) return("Tajikistan-Kyrgyzstan")
  if (lon >= 46 && lon <= 88 && lat >= 40 && lat <= 56) return("Kazakhstan")
  
  # ============================================================================
  # SOUTH ASIA
  # ============================================================================
  
  if (lon >= 60 && lon <= 72 && lat >= 29 && lat <= 39) return("Afghanistan")
  if (lon >= 60 && lon <= 78 && lat >= 23 && lat <= 38) return("Pakistan")
  if (lon >= 68 && lon <= 98 && lat >= 6 && lat <= 37) return("India")
  if (lon >= 88 && lon <= 93 && lat >= 20 && lat <= 27) return("Bangladesh")
  if (lon >= 79 && lon <= 82 && lat >= 5 && lat <= 10) return("Sri-Lanka")
  if (lon >= 80 && lon <= 89 && lat >= 26 && lat <= 31) return("Nepal")
  if (lon >= 88 && lon <= 93 && lat >= 26 && lat <= 29) return("Bhutan")
  if (lon >= 72 && lon <= 74 && lat >= 4 && lat <= 8) return("Maldives")
  
  # ============================================================================
  # EAST ASIA
  # ============================================================================
  
  if (lon >= 73 && lon <= 135 && lat >= 18 && lat <= 54) return("China")
  if (lon >= 105 && lon <= 112 && lat >= 21 && lat <= 24) return("Hong-Kong-Macau")
  if (lon >= 120 && lon <= 122 && lat >= 22 && lat <= 26) return("Taiwan")
  if (lon >= 124 && lon <= 132 && lat >= 33 && lat <= 44) return("Korea")
  if (lon >= 129 && lon <= 146 && lat >= 30 && lat <= 46) return("Japan")
  if (lon >= 87 && lon <= 120 && lat >= 41 && lat <= 53) return("Mongolia")
  
  # ============================================================================
  # SOUTHEAST ASIA
  # ============================================================================
  
  if (lon >= 92 && lon <= 102 && lat >= 10 && lat <= 29) return("Myanmar-Thailand")
  if (lon >= 102 && lon <= 110 && lat >= 10 && lat <= 24) return("Laos-Vietnam")
  if (lon >= 102 && lon <= 108 && lat >= 10 && lat <= 15) return("Cambodia")
  if (lon >= 95 && lon <= 106 && lat >= -6 && lat <= 6) return("Sumatra-W-Malaysia")
  if (lon >= 106 && lon <= 120 && lat >= -8 && lat <= 8) return("Java-Borneo")
  if (lon >= 117 && lon <= 141 && lat >= -11 && lat <= 2) return("Sulawesi-Papua-Indonesia")
  if (lon >= 99 && lon <= 105 && lat >= 1 && lat <= 8) return("Peninsular-Malaysia")
  if (lon >= 103 && lon <= 105 && lat >= 1 && lat <= 2) return("Singapore")
  if (lon >= 114 && lon <= 116 && lat >= 4 && lat <= 6) return("Brunei")
  if (lon >= 117 && lon <= 127 && lat >= 5 && lat <= 22) return("Philippines")
  if (lon >= 96 && lon <= 106 && lat >= -9 && lat <= 0) return("Timor-Leste")
  
  # ============================================================================
  # AMERICAS - NORTH
  # ============================================================================
  
  if (lon >= -170 && lon <= -130 && lat >= 55 && lat <= 72) return("Alaska")
  if (lon >= -141 && lon <= -95 && lat >= 60 && lat <= 84) return("Northern-Canada")
  if (lon >= -130 && lon <= -52 && lat >= 41 && lat <= 60) return("Canada")
  if (lon >= -125 && lon <= -65 && lat >= 25 && lat <= 50) return("USA")
  if (lon >= -170 && lon <= -154 && lat >= 18 && lat <= 23) return("Hawaii")
  
  # CENTRAL AMERICA & CARIBBEAN
  if (lon >= -119 && lon <= -86 && lat >= 14 && lat <= 33) return("Mexico")
  if (lon >= -92 && lon <= -87 && lat >= 13 && lat <= 18) return("Guatemala-Belize")
  if (lon >= -91 && lon <= -83 && lat >= 12 && lat <= 16) return("Honduras-El-Salvador-Nicaragua")
  if (lon >= -86 && lon <= -82 && lat >= 8 && lat <= 12) return("Costa-Rica-Panama")
  if (lon >= -85 && lon <= -60 && lat >= 10 && lat <= 24) return("Caribbean-Islands")
  
  # SOUTH AMERICA
  if (lon >= -74 && lon <= -66 && lat >= -5 && lat <= 13) return("Colombia")
  if (lon >= -73 && lon <= -60 && lat >= 0 && lat <= 13) return("Venezuela-Guyana-Suriname")
  if (lon >= -82 && lon <= -75 && lat >= -5 && lat <= 2) return("Ecuador")
  if (lon >= -82 && lon <= -68 && lat >= -19 && lat <= 0) return("Peru")
  if (lon >= -75 && lon <= -47 && lat >= -34 && lat <= 6) return("Brazil")
  if (lon >= -70 && lon <= -53 && lat >= -24 && lat <= -10) return("Bolivia-Paraguay")
  if (lon >=-76 && lon <= -66 && lat >= -56 && lat <= -17) return("Chile")
  if (lon >= -74 && lon <= -53 && lat >= -56 && lat <= -21) return("Argentina-Uruguay")
  
  # ============================================================================
  # OCEANIA
  # ============================================================================
  
  if (lon >= 110 && lon <= 155 && lat >= -45 && lat <= -10) return("Australia")
  if (lon >= 165 && lon <= 180 && lat >= -48 && lat <= -34) return("New-Zealand")
  if (lon >= 140 && lon <= 180 && lat >= -25 && lat <= 10) return("Pacific-Islands")
  if (lon >= 130 && lon <= 155 && lat >= -12 && lat <= -1) return("Papua-New-Guinea")
  
  return("Other-Unassigned")
}

plot_gebv_heatmap <- function(gebv_matrix, coords, method_name = "", 
                              trait_subset = NULL, use_specific_countries = FALSE) {
  
  # Subset traits if specified
  if (!is.null(trait_subset)) {
    gebv_matrix <- gebv_matrix[, trait_subset, drop = FALSE]
  }
  
  # Assign countries/regions to ALL samples
  sample_ids <- rownames(gebv_matrix)
  
  if (use_specific_countries) {
    countries <- sapply(sample_ids, function(id) {
      if (id %in% rownames(coords)) {
        assign_country_specific(coords[id, "Latitude"], coords[id, "Longitude"])
      } else {
        "Unknown"
      }
    })
    title_prefix <- "Specific Countries"
  } else {
    countries <- sapply(sample_ids, function(id) {
      if (id %in% rownames(coords)) {
        assign_country(coords[id, "Latitude"], coords[id, "Longitude"])
      } else {
        "Unknown"
      }
    })
    title_prefix <- "Geographic Regions"
  }
  
  # Create country color palette
  unique_countries <- sort(unique(countries))
  n_countries <- length(unique_countries)
  
  # Use distinct colors - expand palette for more countries
  country_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                     "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
                     "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                     "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                     "#E6AB02", "#A6761D", "#666666", "#B3DE69", "#FCCDE5")
  
  if (n_countries > length(country_colors)) {
    country_colors <- rainbow(n_countries)
  } else {
    country_colors <- country_colors[1:n_countries]
  }
  names(country_colors) <- unique_countries
  
  # Scale GEBVs by trait (Z-scores)
  gebv_scaled <- scale(gebv_matrix)
  
  # Order by country, then cluster within country
  country_order <- order(countries)
  gebv_ordered <- gebv_scaled[country_order, ]
  countries_ordered <- countries[country_order]
  
  # Within each country, cluster by similarity
  country_blocks <- split(1:length(countries_ordered), countries_ordered)
  final_order <- unlist(lapply(names(country_blocks), function(ctry) {
    idx <- country_blocks[[ctry]]
    if (length(idx) > 1) {
      sub_gebv <- gebv_ordered[idx, , drop = FALSE]
      tryCatch({
        hc <- hclust(dist(sub_gebv))
        idx[hc$order]
      }, error = function(e) {
        idx
      })
    } else {
      idx
    }
  }))
  
  gebv_final <- gebv_ordered[final_order, ]
  countries_final <- countries_ordered[final_order]
  sample_names_final <- rownames(gebv_final)
  
  # Create data frame for main heatmap
  gebv_long <- reshape2::melt(gebv_final)
  colnames(gebv_long) <- c("Sample", "Trait", "Zscore")
  gebv_long$Sample_num <- as.numeric(factor(gebv_long$Sample, 
                                           levels = sample_names_final))
  
  # Create annotation data for country bar
  anno_df <- data.frame(
    Sample_num = 1:length(sample_names_final),
    Country = countries_final,
    x_start = -2,
    x_end = -0.5
  )
  
  # Main heatmap
  p <- ggplot() +
    # Country annotation bar (left side)
    geom_rect(data = anno_df, 
              aes(xmin = x_start, xmax = x_end, 
                  ymin = Sample_num - 0.5, ymax = Sample_num + 0.5,
                  fill = Country),
              color = NA) +
    scale_fill_manual(values = country_colors, 
                     name = ifelse(use_specific_countries, "Country", "Geographic Region")) +
    # Main heatmap with new fill scale
    ggnewscale::new_scale_fill() +
    geom_tile(data = gebv_long, 
              aes(x = as.numeric(factor(Trait)), y = Sample_num, fill = Zscore)) +
    scale_fill_gradient2(low = "#053061", mid = "#F7F7F7", high = "#67001F",
                        midpoint = 0, 
                        limits = c(-3, 3),
                        oob = scales::squish,
                        name = "Z-score") +
    scale_x_continuous(breaks = 1:ncol(gebv_final),
                      labels = colnames(gebv_final),
                      expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-2, ncol(gebv_final) + 0.5)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(size = 10),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    ) +
    labs(
      title = paste0("GEBVs Heatmap Grouped by ", title_prefix),
      subtitle = paste0(nrow(gebv_final), " accessions .. ", ncol(gebv_final), " traits"),
      x = "Trait",
      y = "Accession"
    )
  
  return(p)
}
################################################################################
# MAIN PIPELINE
################################################################################

run_genomic_pipeline_optimized <- function(vcf_file, coords_file, worldclim_path,
                                          n_core = 25, output_dir = "results",
                                          k_fold = 10, n_reps = 25,
                                          methods = c("rrBLUP"),
                                          n_cores = NULL, max_snps_core = 10000,
                                          use_all_tifs = FALSE,
                                          max_missing = 0,
                                          impute_method = "mean") {
  
  if (is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("================================================================================\n")
  cat("                    GENOMIC SELECTION PIPELINE v1.3                            \n")
  cat("================================================================================\n\n")
  
  # STEP 0: Match samples between VCF and coordinates
  cat("STEP 0: Matching samples...\n")
  coords <- read.csv(coords_file)
  coords <- coords[!duplicated(coords$ID), ]
  
  gds_temp <- tempfile(fileext = ".gds")
  snpgdsVCF2GDS(vcf_file, gds_temp, method = "biallelic.only", verbose = FALSE)
  genofile <- snpgdsOpen(gds_temp)
  vcf_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
  snpgdsClose(genofile)
  unlink(gds_temp)
  
  cat("  Total VCF samples:", length(vcf_samples), "\n")
  cat("  Samples with coordinates:", length(coords$ID), "\n")
  
  samples_with_coords <- intersect(vcf_samples, coords$ID)
  cat("  Samples with both genotypes AND coordinates:", length(samples_with_coords), "\n\n")
  
  if (length(samples_with_coords) == 0) {
    stop("ERROR: No samples have both genotypes and coordinates!")
  }
  
  # STEP 1: Convert VCF for ALL samples
  cat("STEP 1: Converting VCF (ALL samples)...\n")
  geno_all <- vcf_to_rrblup_fast(vcf_file, 
                                 file.path(output_dir, "genotypes_all.txt"),
                                 sample_subset = vcf_samples,  # ALL VCF samples
                                 n_threads = n_cores)
  
  # STEP 2: Extract environmental data for samples with coordinates
  cat("STEP 2: Extracting environmental data...\n")
  coords_filtered <- coords[coords$ID %in% samples_with_coords, ]
  env_data <- extract_worldclim_data(coords_filtered, worldclim_path, 
                                     use_all_tifs = use_all_tifs,
                                     output_file = file.path(output_dir, "env_data.csv"))
  
  # STEP 3: Identify training samples (complete environmental data)
  cat("STEP 3: Identifying training samples (with complete environmental data)...\n")
  rownames(env_data) <- env_data$ID
  
  # Extract ALL environmental/climate columns (all numeric except ID, Longitude, Latitude)
  exclude_cols <- c("ID", "Longitude", "Latitude")
  pheno_cols <- which(!colnames(env_data) %in% exclude_cols & 
                      sapply(env_data, is.numeric))
  
  if (length(pheno_cols) == 0) {
    stop("ERROR: No environmental data columns found!")
  }
  
  pheno_all <- as.matrix(env_data[, pheno_cols])
  
  cat("  Total environmental variables extracted:", ncol(pheno_all), "\n")
  cat("  Variable names:", paste(head(colnames(pheno_all), 10), collapse=", "))
  if (ncol(pheno_all) > 10) cat(", ...")
  cat("\n")
  cat("  Samples with coordinates:", nrow(pheno_all), "\n")
  
  # Find samples with COMPLETE environmental data (traits should be complete)
  missing_per_sample <- rowSums(is.na(pheno_all))
  samples_complete_env <- rownames(pheno_all)[missing_per_sample == 0]
  
  cat("  Samples with COMPLETE environmental data (all traits):", length(samples_complete_env), "\n")
  
  # Find samples that have BOTH complete env data AND genotypes
  samples_both <- intersect(samples_complete_env, rownames(geno_all))
  
  cat("  Samples with BOTH complete env data AND genotypes:", length(samples_both), "\n")
  
  if (length(samples_both) < n_core) {
    stop("ERROR: Only ", length(samples_both), 
         " samples have both complete env data and genotypes, but you requested ", 
         n_core, " core samples!")
  }
  
  # Filter to complete samples for training
  pheno_complete <- pheno_all[samples_both, ]
  geno_complete <- geno_all[samples_both, ]
  
  # Filter SNPs based on missing data tolerance
  total_snps <- ncol(geno_complete)
  missing_per_snp <- colSums(is.na(geno_complete))
  
  if (max_missing == 0) {
    # Strict: no missing genotypes allowed
    snps_keep <- missing_per_snp == 0
    cat("  SNP filtering: Strict (0% missing allowed)\n")
  } else if (max_missing < 1) {
    # Proportion: e.g., max_missing=0.1 means allow SNPs with up to 10% missing samples
    max_missing_samples <- ceiling(nrow(geno_complete) * max_missing)
    snps_keep <- missing_per_snp <= max_missing_samples
    cat("  SNP filtering: Allow ≤", round(max_missing*100, 1), "% missing (≤", 
        max_missing_samples, "samples)\n")
  } else {
    # Absolute count: e.g., max_missing=5 means allow SNPs missing in up to 5 samples
    snps_keep <- missing_per_snp <= max_missing
    cat("  SNP filtering: Allow SNPs missing in ≤", max_missing, "samples\n")
  }
  
  geno_complete <- geno_complete[, snps_keep]
  
  # Show SNP filtering results
  cat("  SNPs before filtering:", total_snps, "\n")
  cat("  SNPs after filtering:", ncol(geno_complete), 
      "(", round(100*ncol(geno_complete)/total_snps, 1), "% retained)\n")
  
  if (ncol(geno_complete) < 1000) {
    warning("Very few SNPs remaining (", ncol(geno_complete), "). Consider increasing max_missing.")
  }
  
  cat("  Training pool: ", nrow(geno_complete), " samples, ", ncol(geno_complete), 
      " SNPs, ", ncol(pheno_complete), " traits\n\n", sep="")
  
  # STEP 4: Core selection from samples with complete data
  cat("STEP 4: Selecting core from samples with complete environmental data...\n")
  if (n_core > nrow(geno_complete)) {
    n_core <- floor(nrow(geno_complete) * 0.8)
    cat("  Adjusted n_core to", n_core, "\n")
  }
  
  core_ids <- greedy_core_selection_fast(geno_complete, n_core, 
                                         max_snps = max_snps_core, 
                                         n_threads = n_cores)
  
  # Save core info
  core_info <- env_data[core_ids, c("ID", "Longitude", "Latitude")]
  write.csv(core_info, file.path(output_dir, "core_ids.csv"), row.names = FALSE)
  cat("  ✓ Selected", length(core_ids), "core samples (all with complete env data)\n")
  cat("  ✓ Core IDs saved with coordinates\n\n")
  
  # STEP 5: Prepare training and prediction sets
  cat("STEP 5: Preparing training and prediction sets...\n")
  
  # Training: core samples with complete data
  geno_train <- geno_complete[core_ids, ]
  pheno_train <- pheno_complete[core_ids, ]
  
  # Prediction: ALL samples from VCF (align SNPs with training)
  common_snps <- intersect(colnames(geno_train), colnames(geno_all))
  geno_train <- geno_train[, common_snps]
  geno_pred <- geno_all[, common_snps]
  
  cat("  Training set:", nrow(geno_train), "samples (core with complete env data)\n")
  cat("  Prediction set:", nrow(geno_pred), "samples (ALL from VCF)\n")
  cat("  Common SNPs:", ncol(geno_train), "\n")
  cat("  Traits:", ncol(pheno_train), "\n\n")
  
  # Verify all core samples have complete data
  if (any(is.na(pheno_train))) {
    stop("ERROR: Core training samples have missing environmental data!")
  }
  
  # STEP 5b: Impute missing genotypes in prediction set
  cat("STEP 5b: Imputing missing genotypes in prediction set...\n")
  
  # Check for missing data
  n_missing_train <- sum(is.na(geno_train))
  n_missing_pred <- sum(is.na(geno_pred))
  
  cat("  Missing genotypes in training set:", n_missing_train, 
      "(", round(100*n_missing_train/(nrow(geno_train)*ncol(geno_train)), 2), "%)\n")
  cat("  Missing genotypes in prediction set:", n_missing_pred, 
      "(", round(100*n_missing_pred/(nrow(geno_pred)*ncol(geno_pred)), 2), "%)\n")
  
  if (n_missing_train > 0 || n_missing_pred > 0) {
    cat("  Imputation method:", impute_method, "\n")
    
    if (impute_method == "mean") {
      # Mean imputation
      # Impute training set
      if (n_missing_train > 0) {
        for (j in 1:ncol(geno_train)) {
          if (any(is.na(geno_train[, j]))) {
            mean_val <- mean(geno_train[, j], na.rm = TRUE)
            geno_train[is.na(geno_train[, j]), j] <- mean_val
          }
        }
      }
      
      # Impute prediction set using training set means
      if (n_missing_pred > 0) {
        for (j in 1:ncol(geno_pred)) {
          if (any(is.na(geno_pred[, j]))) {
            # Use mean from training set for this SNP
            if (colnames(geno_pred)[j] %in% colnames(geno_train)) {
              mean_val <- mean(geno_train[, colnames(geno_pred)[j]], na.rm = TRUE)
            } else {
              mean_val <- mean(geno_pred[, j], na.rm = TRUE)
            }
            geno_pred[is.na(geno_pred[, j]), j] <- mean_val
          }
        }
      }
      
    } else if (impute_method == "median") {
      # Median imputation
      if (n_missing_train > 0) {
        for (j in 1:ncol(geno_train)) {
          if (any(is.na(geno_train[, j]))) {
            median_val <- median(geno_train[, j], na.rm = TRUE)
            geno_train[is.na(geno_train[, j]), j] <- median_val
          }
        }
      }
      
      if (n_missing_pred > 0) {
        for (j in 1:ncol(geno_pred)) {
          if (any(is.na(geno_pred[, j]))) {
            if (colnames(geno_pred)[j] %in% colnames(geno_train)) {
              median_val <- median(geno_train[, colnames(geno_pred)[j]], na.rm = TRUE)
            } else {
              median_val <- median(geno_pred[, j], na.rm = TRUE)
            }
            geno_pred[is.na(geno_pred[, j]), j] <- median_val
          }
        }
      }
      
    } else if (impute_method == "mode") {
      # Mode imputation (most common value)
      get_mode <- function(x) {
        x <- x[!is.na(x)]
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }
      
      if (n_missing_train > 0) {
        for (j in 1:ncol(geno_train)) {
          if (any(is.na(geno_train[, j]))) {
            mode_val <- get_mode(geno_train[, j])
            geno_train[is.na(geno_train[, j]), j] <- mode_val
          }
        }
      }
      
      if (n_missing_pred > 0) {
        for (j in 1:ncol(geno_pred)) {
          if (any(is.na(geno_pred[, j]))) {
            if (colnames(geno_pred)[j] %in% colnames(geno_train)) {
              mode_val <- get_mode(geno_train[, colnames(geno_pred)[j]])
            } else {
              mode_val <- get_mode(geno_pred[, j])
            }
            geno_pred[is.na(geno_pred[, j]), j] <- mode_val
          }
        }
      }
      
    } else {
      stop("Unknown imputation method: ", impute_method, 
           ". Choose from: 'mean', 'median', 'mode'")
    }
    
    cat("  ✓ Imputation complete\n")
    cat("  Training set: 0 missing (", nrow(geno_train), " samples × ", 
        ncol(geno_train), " SNPs)\n", sep="")
    cat("  Prediction set: 0 missing (", nrow(geno_pred), " samples × ", 
        ncol(geno_pred), " SNPs)\n\n", sep="")
  } else {
    cat("  No missing data - imputation not needed\n\n")
  }
  
  # STEP 6: Cross-validation (on training set only)
  cat("STEP 6: Cross-validation on training set...\n")
  cv_results <- list()
  for (method in methods) {
    cat("\n  Method:", method, "\n")
    cv_results[[method]] <- kfold_cross_validation_parallel(
      geno_train, pheno_train, pheno_train, k_fold, n_reps, method, n_cores = n_cores)
    saveRDS(cv_results[[method]], file.path(output_dir, paste0("cv_", method, ".rds")))
  }
  
  # STEP 7: Genomic selection (train on core, predict on ALL VCF samples)
  cat("\nSTEP 7: Genomic selection (train on core, predict on ALL VCF samples)...\n")
  cat("  Training samples:", nrow(geno_train), "(all with complete environmental data)\n")
  cat("  Prediction samples:", nrow(geno_pred), "(all samples from original VCF)\n\n")
  
  gs_results_all <- list()
  for (method in methods) {
    cat("  Running", method, "...\n")
    gs_results_all[[method]] <- genomic_selection(
      geno_train = geno_train, 
      geno_pred = geno_pred,
      pheno_train = pheno_train,
      method = method
    )
    
    # Save method-specific results
    write.csv(gs_results_all[[method]]$gebv_all, 
              file.path(output_dir, paste0("GEBVs_", method, "_all.csv")))
    write.csv(gs_results_all[[method]]$gebv_train, 
              file.path(output_dir, paste0("GEBVs_", method, "_training.csv")))
  }
  
  cat("\n  ✓ GEBVs for all", nrow(geno_pred), "VCF samples saved\n")
  cat("  ✓ Results saved for", length(methods), "methods\n\n")
  
  # Plot CV results
  cat("\nGenerating plots...\n")
  p_cv <- plot_cv_results(cv_results)
  ggsave(file.path(output_dir, "cv_plot.png"), p_cv, width = 14, height = 10, dpi = 300)
  cat("  ✓ CV plot saved\n")
  
  # Plot GEBV heatmaps for each method
  cat("  Generating GEBV heatmaps (all samples, all traits)...\n")
  
  # Prepare coordinates matrix with sample IDs as rownames
  coords_matrix <- coords_filtered
  rownames(coords_matrix) <- coords_matrix$ID
  
  for (method in methods) {
    cat("    -", method, "heatmaps...\n")
    gebv_data <- gs_results_all[[method]]$gebv_all
    
    cat("      (showing all", ncol(gebv_data), "traits)\n")
    
    # Regional heatmap - ALL TRAITS
    tryCatch({
      cat("      - Regional grouping...\n")
      p_regional <- plot_gebv_heatmap(gebv_data, coords_matrix, 
                                      method_name = method,
                                      trait_subset = NULL,  # Use ALL traits
                                      use_specific_countries = FALSE)
      ggsave(file.path(output_dir, paste0("gebv_heatmap_regional_", method, ".png")), 
             p_regional, width = 18, height = 12, dpi = 300)
    }, error = function(e) {
      cat("        Warning: Could not create regional heatmap:", e$message, "\n")
    })
    
    # Country-specific heatmap - ALL TRAITS
    tryCatch({
      cat("      - Country-specific grouping...\n")
      p_country <- plot_gebv_heatmap(gebv_data, coords_matrix, 
                                     method_name = method,
                                     trait_subset = NULL,  # Use ALL traits
                                     use_specific_countries = TRUE)
      ggsave(file.path(output_dir, paste0("gebv_heatmap_country_", method, ".png")), 
             p_country, width = 18, height = 12, dpi = 300)
    }, error = function(e) {
      cat("        Warning: Could not create country heatmap:", e$message, "\n")
    })
  }
  cat("  ✓ GEBV heatmaps saved (both regional and country-specific, all traits)\n")
  
  cat("\n================================================================================\n")
  cat("                           PIPELINE COMPLETE!                                   \n")
  cat("================================================================================\n")
  cat("Results saved to:", output_dir, "\n")
  cat("Methods compared:", paste(methods, collapse = ", "), "\n")
  cat("\nOutput files:\n")
  cat("  - cv_plot.png (cross-validation accuracy)\n")
  for (method in methods) {
    cat("  - gebv_heatmap_regional_", method, ".png (GEBVs by region)\n", sep="")
    cat("  - gebv_heatmap_country_", method, ".png (GEBVs by country)\n", sep="")
  }
  cat("\n")
  
  return(list(geno_train = geno_train, geno_pred = geno_pred, 
              pheno_train = pheno_train, core_ids = core_ids,
              cv_results = cv_results, 
              gebvs = gs_results_all))
}

################################################################################
# SCRIPT LOADED
################################################################################

cat("\n================================================================================\n")
cat("  Genomic Selection Pipeline v1.8 - LOADED\n")
cat("================================================================================\n")
cat("\nFeatures:\n")
cat("  ✓ Fixed VCF conversion\n")
cat("  ✓ Automatic sample matching\n")
cat("  ✓ Parallel processing\n")
cat("  ✓ Flexible SNP missing data tolerance\n")
cat("  ✓ Genotype imputation (mean/median/mode)\n")
cat("  ✓ Multiple prediction methods:\n")
cat("      - rrBLUP (ridge regression)\n")
cat("      - GBLUP (genomic BLUP)\n")
cat("      - GBLUP-Gaussian (Gaussian kernel)\n")
cat("      - GBLUP-Exponential (Exponential kernel)\n")
cat("      - BayesCpi (Bayesian variable selection)\n")
cat("  ✓ Comparative CV visualization\n")
cat("  ✓ Dual GEBV heatmaps (regional + country-specific)\n")
cat("  ✓ All traits included in heatmaps\n\n")
cat("Usage:\n")
cat("  results <- run_genomic_pipeline_optimized(\n")
cat("    vcf_file = 'genotypes.vcf.gz',\n")
cat("    coords_file = 'coordinates.csv',\n")
cat("    worldclim_path = 'worldclim/',\n")
cat("    n_core = 50,\n")
cat("    output_dir = 'results',\n")
cat("    methods = c('rrBLUP', 'GBLUP-Gaussian', 'BayesCpi'),\n")
cat("    max_missing = 0.05,  # Allow SNPs with up to 5% missing data\n")
cat("    impute_method = 'mean',  # 'mean', 'median', or 'mode'\n")
cat("    n_cores = 10,\n")
cat("    use_all_tifs = TRUE\n")
cat("  )\n\n")
cat("Note: Heatmaps will include ALL environmental traits\n")
cat("      Install BGLR for BayesCpi: install.packages('BGLR')\n\n")
