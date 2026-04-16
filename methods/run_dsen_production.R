#!/usr/bin/env Rscript
# =============================================================================
# DSEN Production Run: Comprehensive EN vs DSEN Comparison
# =============================================================================
#
# Purpose: Run full analysis comparing Elastic Net and DSEN with:
#   - CV with alpha tuning for both methods
#   - Bootstrap for feature selection stability
#   - Comprehensive metrics saved for post-hoc analysis
#
# Output: results/dsen_results/dsen_production_results.Rda
#
# Usage:
#   # Test mode (2 drugs, 10 bootstrap reps):
#   Rscript methods/run_dsen_production.R
#
#   # Production mode (all drugs, 100 bootstrap reps):
#   TEST_MODE=FALSE Rscript methods/run_dsen_production.R
#
#   # Custom settings:
#   TEST_MODE=FALSE BOOTSTRAP_N=50 NCORES=8 Rscript methods/run_dsen_production.R
#
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "03_dsen_utils.R"))
source(here("methods", "05_parallel_utils.R"))

# =============================================================================
# CONFIGURATION
# =============================================================================
TEST_MODE <- Sys.getenv("TEST_MODE", "TRUE") == "TRUE"

max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) {
  as.integer(max_drugs_env)
} else if (TEST_MODE) {
  2L
} else {
  Inf
}

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "3"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) {
  CHECKPOINT_EVERY <- 3L
}

# Minimum samples
MIN_SAMPLES <- as.integer(Sys.getenv("MIN_SAMPLES", "50"))
MIN_SAMPLES_PER_TISSUE <- as.integer(Sys.getenv("MIN_SAMPLES_PER_TISSUE", "15"))
MIN_TISSUES <- as.integer(Sys.getenv("MIN_TISSUES", "3"))

# CV settings
CV_NFOLDS <- as.integer(Sys.getenv("CV_NFOLDS", if (TEST_MODE) "3" else "5"))
if (is.na(CV_NFOLDS) || CV_NFOLDS < 3) CV_NFOLDS <- 3L

# Bootstrap settings
SKIP_BOOTSTRAP <- Sys.getenv("SKIP_BOOTSTRAP", "FALSE") == "TRUE"
BOOTSTRAP_NREPS <- as.integer(Sys.getenv("BOOTSTRAP_N", if (TEST_MODE) "10" else "200"))
if (is.na(BOOTSTRAP_NREPS) || BOOTSTRAP_NREPS < 5) {
  BOOTSTRAP_NREPS <- if (TEST_MODE) 10L else 200L
}

# Feature screening settings
CORRELATION_THRESHOLD <- as.numeric(Sys.getenv("CORRELATION_THRESHOLD", "0.1"))
VARIANCE_THRESHOLD <- as.numeric(Sys.getenv("VARIANCE_THRESHOLD", "0.01"))

# Model tuning grids
# For fast runs: ALPHA_GRID=0.5 LAMRAT_GRID=2,4
alpha_env <- Sys.getenv("ALPHA_GRID", "")
ALPHA_GRID <- if (nzchar(alpha_env)) {
  as.numeric(strsplit(alpha_env, ",")[[1]])
} else {
  c(0.05, 0.1, 0.15, 0.2, seq(0.3, 1.0, 0.1))
}

lamrat_env <- Sys.getenv("LAMRAT_GRID", "")
LAMRAT_GRID <- if (nzchar(lamrat_env)) {
  as.numeric(strsplit(lamrat_env, ",")[[1]])
} else {
  c(0.125, 0.25, 0.5, 1, 2, 4)
}

# DSEN settings
# Note: union_threshold is best but very slow (~50k features)
# Use pooled for production feasibility (~2k features), can re-run select drugs later
SCREENING <- Sys.getenv("SCREENING", "pooled")  # pooled or union_threshold
TISSUE_PENALTY_MODE <- Sys.getenv("TISSUE_PENALTY_MODE", "sample_size_mild")

GLOBAL_SEED <- as.integer(Sys.getenv("GLOBAL_SEED", "42"))

NCORES <- as.integer(Sys.getenv("NCORES", NA))
if (is.na(NCORES)) {
  NCORES <- max(1L, parallel::detectCores() - 1L)
}

# Output directory
results_subdir <- Sys.getenv("RESULTS_SUBDIR", "")
results_dir <- if (nzchar(results_subdir)) {
  here("results", "dsen_results", results_subdir)
} else {
  here("results", "dsen_results")
}
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

cat(strrep("=", 80), "\n")
cat("DSEN PRODUCTION RUN: EN vs DSEN Comparison\n")
cat(strrep("=", 80), "\n")
cat("Mode:", if (TEST_MODE) "TEST" else "PRODUCTION", "\n")
cat("Max drugs:", if (is.infinite(MAX_DRUGS)) "all" else MAX_DRUGS, "\n")
cat("CV folds:", CV_NFOLDS, "\n")
cat("Bootstrap:", if (SKIP_BOOTSTRAP) "SKIPPED" else paste(BOOTSTRAP_NREPS, "reps"), "\n")
cat("Alpha grid:", paste(ALPHA_GRID, collapse = ", "), "\n")
cat("Lamrat grid:", paste(LAMRAT_GRID, collapse = ", "), "\n")
cat("Screening:", SCREENING, "\n")
cat("Tissue penalty mode:", TISSUE_PENALTY_MODE, "\n")
cat("Results dir:", results_dir, "\n")
cat("Checkpoint every:", CHECKPOINT_EVERY, "drugs\n")
cat("Cores:", NCORES, "\n")
cat(strrep("=", 80), "\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading data...\n")
data <- load_main_data("Sanger")
all_drugs <- get_drug_list(data)

if (!is.infinite(MAX_DRUGS)) {
  all_drugs <- all_drugs[1:min(MAX_DRUGS, length(all_drugs))]
  cat("Limited to", length(all_drugs), "drugs\n")
}

cat("Total drugs to process:", length(all_drugs), "\n\n")

# =============================================================================
# CHECKPOINT SYSTEM
# =============================================================================
checkpoint_file <- file.path(results_dir, "dsen_production_checkpoint.Rda")
results_file <- file.path(results_dir, "dsen_production_results.Rda")

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  cat("Resuming from checkpoint\n")
  cat("  Completed:", length(completed_drugs), "drugs\n")
  remaining_drugs <- setdiff(all_drugs, completed_drugs)
} else {
  completed_drugs <- character(0)
  all_results <- list()
  remaining_drugs <- all_drugs
  cat("Starting fresh\n")
}

cat("Remaining:", length(remaining_drugs), "drugs\n\n")

# =============================================================================
# HELPER: Jaccard similarity
# =============================================================================
jaccard <- function(set1, set2) {
  if (length(set1) == 0 && length(set2) == 0) return(NA)
  inter <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  if (union == 0) return(NA)
  return(inter / union)
}

# =============================================================================
# DRUG PROCESSING FUNCTION
# =============================================================================
process_drug <- function(drug_name, data, seed) {

  # Prepare data with per-tissue filtering
  drug_data <- prepare_drug_data_filtered(
    data,
    drug_name,
    min_per_tissue = MIN_SAMPLES_PER_TISSUE
  )

  if (is.null(drug_data)) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = NA, n_tissues = NA,
      mse_en = NA, mse_dsen = NA, pct_improvement = NA,
      alpha_en = NA, alpha_dsen = NA, lamrat_dsen = NA,
      lambda_en = NA, lambda_dsen = NA,
      en_bootstrap = NULL, dsen_bootstrap = NULL,
      error = NA
    ))
  }

  x <- drug_data$x
  y <- drug_data$y
  tissue <- drug_data$tissue
  n_tissues <- length(unique(tissue))

  # Check minimum tissues for DSEN
  if (n_tissues < MIN_TISSUES) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = length(y), n_tissues = n_tissues,
      mse_en = NA, mse_dsen = NA, pct_improvement = NA,
      alpha_en = NA, alpha_dsen = NA, lamrat_dsen = NA,
      lambda_en = NA, lambda_dsen = NA,
      en_bootstrap = NULL, dsen_bootstrap = NULL,
      error = NA
    ))
  }

  # -------------------------------------------------------------------------
  # Part 1: Elastic Net
  # -------------------------------------------------------------------------
  set.seed(seed)
  folds <- createFolds(tissue, k = CV_NFOLDS, list = FALSE)

  en_result <- cv_tune_alpha(
    x, y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed,
    foldid = folds
  )

  mse_en <- en_result$best_cvm
  alpha_en <- en_result$alpha.min
  lambda_en <- en_result$lambda.min

  # EN Bootstrap (optional)
  en_bs <- if (SKIP_BOOTSTRAP) {
    list(cnt = NULL)
  } else {
    bootstrap_features(
      x, y,
      alpha = alpha_en,
      lambda = lambda_en,
      nboot = BOOTSTRAP_NREPS,
      method = "right",
      seed = seed
    )
  }

  # -------------------------------------------------------------------------
  # Part 2: DSEN
  # -------------------------------------------------------------------------
  dsen_result <- tryCatch({
    dsen_tune_alpha_lamrat(
      x, y, tissue,
      alphas = ALPHA_GRID,
      lamrats = LAMRAT_GRID,
      nfolds = CV_NFOLDS,
      tissue_penalty_mode = TISSUE_PENALTY_MODE,
      screening = SCREENING,
      seed = seed,
      foldid = folds
    )
  }, error = function(e) {
    list(best_cvm = NA, alpha.min = NA, lamrat.min = NA, error = e$message)
  })

  if (is.na(dsen_result$best_cvm)) {
    return(list(
      drug = drug_name,
      status = "dsen_error",
      n_samples = length(y), n_tissues = n_tissues,
      mse_en = mse_en, mse_dsen = NA, pct_improvement = NA,
      alpha_en = alpha_en, alpha_dsen = NA, lamrat_dsen = NA,
      lambda_en = lambda_en, lambda_dsen = NA,
      en_bootstrap = en_bs$cnt,  # Raw counts
      dsen_bootstrap = NULL,
      error = dsen_result$error
    ))
  }

  mse_dsen <- dsen_result$best_cvm
  alpha_dsen <- dsen_result$alpha.min
  lamrat_dsen <- dsen_result$lamrat.min
  lambda_dsen <- dsen_result$lambda.min

  # DSEN Bootstrap (optional)
  dsen_bs <- if (SKIP_BOOTSTRAP) {
    list(cnt = NULL)
  } else {
    tryCatch({
      dsen_bootstrap(
        x, y, tissue,
        alpha = alpha_dsen,
        lambda = lambda_dsen,
        lamrat = lamrat_dsen,
        nboot = BOOTSTRAP_NREPS,
        tissue_penalty_mode = TISSUE_PENALTY_MODE,
        seed = seed
      )
    }, error = function(e) {
      list(cnt = NULL, mean_coefs = NULL, error = e$message)
    })
  }

  # -------------------------------------------------------------------------
  # Part 3: Compute metrics
  # -------------------------------------------------------------------------
  pct_improvement <- (mse_en - mse_dsen) / mse_en * 100

  return(list(
    drug = drug_name,
    status = "success",
    n_samples = length(y),
    n_tissues = n_tissues,

    # MSE
    mse_en = mse_en,
    mse_dsen = mse_dsen,
    pct_improvement = pct_improvement,

    # Hyperparameters
    alpha_en = alpha_en,
    alpha_dsen = alpha_dsen,
    lamrat_dsen = lamrat_dsen,
    lambda_en = lambda_en,
    lambda_dsen = lambda_dsen,

    # Bootstrap results (raw counts for flexible thresholding)
    en_bootstrap = en_bs$cnt,           # Named vector: feature -> count
    dsen_bootstrap = dsen_bs$cnt,       # Named matrix or list

    error = NA
  ))
}

# =============================================================================
# MAIN COMPUTATION
# =============================================================================
start_time <- Sys.time()

if (length(remaining_drugs) > 0) {
  cat("Setting up parallel processing...\n")
  cl <- setup_parallel(NCORES)
  cat("Using", NCORES, "cores\n\n")

  batches <- split(
    remaining_drugs,
    ceiling(seq_along(remaining_drugs) / CHECKPOINT_EVERY)
  )

  for (batch_idx in seq_along(batches)) {
    batch <- batches[[batch_idx]]
    cat(sprintf("\n[%s] Processing batch %d/%d (%d drugs)\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, length(batches), length(batch)))

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret", "Matrix"),
      .export = c(
        "process_drug",
        "prepare_drug_data_filtered",
        "preprocess_features", "as_numeric_matrix",
        "cv_tune_alpha", "bootstrap_features",
        "dsen_tune_alpha_lamrat", "dsen_bootstrap", "dsen_cv",
        "compute_dsen_penalties", "build_dsen_matrix", "prepare_tissue_data", "dsen_predict",
        "ALPHA_GRID", "LAMRAT_GRID", "CV_NFOLDS",
        "VARIANCE_THRESHOLD", "CORRELATION_THRESHOLD",
        "BOOTSTRAP_NREPS", "TISSUE_PENALTY_MODE",
        "MIN_SAMPLES_PER_TISSUE", "MIN_TISSUES",
        "SCREENING", "SKIP_BOOTSTRAP", "GLOBAL_SEED", "data"
      ),
      .combine = c,
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      tryCatch({
        list(process_drug(drug_name, data, seed = GLOBAL_SEED))
      }, error = function(e) {
        list(list(
          drug = drug_name,
          status = "error",
          n_samples = NA, n_tissues = NA,
          mse_en = NA, mse_dsen = NA, pct_improvement = NA,
          alpha_en = NA, alpha_dsen = NA, lamrat_dsen = NA,
          lambda_en = NA, lambda_dsen = NA,
          en_bootstrap = NULL, dsen_bootstrap = NULL,
          error = e$message
        ))
      })
    }

    for (res in batch_results) {
      all_results[[res$drug]] <- res
    }
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    batch_success <- sum(sapply(batch_results, function(r) r$status == "success"))
    batch_skipped <- sum(sapply(batch_results, function(r) grepl("^skipped_", r$status)))
    batch_errors <- sum(sapply(batch_results, function(r) r$status %in% c("error", "dsen_error")))

    cat(sprintf("[%s] Batch %d complete: %d success, %d skipped, %d errors\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, batch_success, batch_skipped, batch_errors))
    cat(sprintf("Overall progress: %d/%d drugs (%.1f%%)\n",
                length(completed_drugs), length(all_drugs),
                100 * length(completed_drugs) / length(all_drugs)))

    elapsed_mins <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate <- length(completed_drugs) / elapsed_mins
    remaining <- length(all_drugs) - length(completed_drugs)
    eta_mins <- if (rate > 0) remaining / rate else NA
    cat(sprintf("Rate: %.2f drugs/min | ETA: %.1f hours\n\n",
                rate, eta_mins / 60))
  }

  stopCluster(cl)
  registerDoSEQ()
}

# =============================================================================
# FINALIZE RESULTS
# =============================================================================
cat("\n", strrep("=", 80), "\n")
cat("FINALIZING RESULTS\n")
cat(strrep("=", 80), "\n\n")

# Build summary dataframe
summary_df <- do.call(rbind, lapply(all_results, function(res) {
  data.frame(
    drug = res$drug,
    status = res$status,
    n_samples = res$n_samples,
    n_tissues = res$n_tissues,
    mse_en = res$mse_en,
    mse_dsen = res$mse_dsen,
    pct_improvement = res$pct_improvement,
    alpha_en = res$alpha_en,
    alpha_dsen = res$alpha_dsen,
    lamrat_dsen = res$lamrat_dsen,
    lambda_en = res$lambda_en,
    lambda_dsen = res$lambda_dsen,
    stringsAsFactors = FALSE
  )
}))

# Extract bootstrap results
en_bootstrap <- lapply(all_results, function(res) res$en_bootstrap)
dsen_bootstrap <- lapply(all_results, function(res) res$dsen_bootstrap)

# Save configuration
config <- list(
  test_mode = TEST_MODE,
  cv_nfolds = CV_NFOLDS,
  bootstrap_nreps = BOOTSTRAP_NREPS,
  alpha_grid = ALPHA_GRID,
  lamrat_grid = LAMRAT_GRID,
  screening = SCREENING,
  tissue_penalty_mode = TISSUE_PENALTY_MODE,
  correlation_threshold = CORRELATION_THRESHOLD,
  variance_threshold = VARIANCE_THRESHOLD,
  min_samples = MIN_SAMPLES,
  min_samples_per_tissue = MIN_SAMPLES_PER_TISSUE,
  min_tissues = MIN_TISSUES,
  seed = GLOBAL_SEED,
  timestamp = Sys.time()
)

# Save comprehensive results
production_results <- list(
  summary_df = summary_df,
  en_bootstrap = en_bootstrap,
  dsen_bootstrap = dsen_bootstrap,
  config = config
)

save(production_results, file = results_file)
cat("Results saved to:", results_file, "\n")

# Also save summary CSV
csv_file <- file.path(results_dir, "dsen_summary.csv")
write.csv(summary_df, csv_file, row.names = FALSE)
cat("Summary CSV saved to:", csv_file, "\n")

# Clean up checkpoint
if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
  cat("Checkpoint cleaned up\n")
}

# Print summary
cat("\n=== SUMMARY ===\n")
cat("Total drugs:", nrow(summary_df), "\n")
cat("Successful:", sum(summary_df$status == "success"), "\n")
cat("Skipped (tissues):", sum(summary_df$status == "skipped_insufficient_tissues"), "\n")
cat("Errors:", sum(summary_df$status %in% c("error", "dsen_error")), "\n")

if (sum(summary_df$status == "success") > 0) {
  success_df <- summary_df[summary_df$status == "success", ]
  cat("\nDSEN vs EN (successful drugs):\n")
  cat("  DSEN wins:", sum(success_df$pct_improvement > 0, na.rm = TRUE), "\n")
  cat("  EN wins:", sum(success_df$pct_improvement <= 0, na.rm = TRUE), "\n")
  cat("  Mean improvement:", round(mean(success_df$pct_improvement, na.rm = TRUE), 2), "%\n")
  cat("  Median improvement:", round(median(success_df$pct_improvement, na.rm = TRUE), 2), "%\n")
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\nElapsed time:", round(elapsed, 1), "minutes\n")
cat("Finished at", format(Sys.time()), "\n")
