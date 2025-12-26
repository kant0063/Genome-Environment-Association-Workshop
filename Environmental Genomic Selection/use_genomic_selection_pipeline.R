# Load pipeline
source("genomic_pipeline.R")

results <- run_genomic_pipeline_optimized(
  vcf_file = " ",
  coords_file = "",
  worldclim_path = " ",
  n_core = 50,
  output_dir = "lentil_results",
  methods = c("rrBLUP", "GBLUP-Gaussian", "GBLUP-Exponential", "BayesCpi"),
  n_cores = 10,  # Use all CPU cores
  use_all_tifs = TRUE  # ← Use all 36 files!
)
