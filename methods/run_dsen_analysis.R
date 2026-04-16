#!/usr/bin/env Rscript
# Purpose: Run DSEN analysis across drugs with checkpointing.
# Saves to: paths$scratch/dsen_results.Rda

suppressPackageStartupMessages({
  library(here)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "03_dsen_utils.R"))
source(here("methods", "05_parallel_utils.R"))

# ==============================================================================
# CONFIGURATION
# ==============================================================================
TEST_MODE <- Sys.getenv("TEST_MODE", "TRUE") == "TRUE"

max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) {
  as.integer(max_drugs_env)
} else if (TEST_MODE) {
  2L
} else {
  Inf
}

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "5"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) {
  CHECKPOINT_EVERY <- 5L
}

MIN_SAMPLES_PER_TISSUE <- as.integer(Sys.getenv("MIN_SAMPLES_PER_TISSUE", "15"))
MIN_TISSUES <- as.integer(Sys.getenv("MIN_TISSUES", "3"))

CV_NFOLDS <- as.integer(Sys.getenv("CV_NFOLDS", if (TEST_MODE) "3" else "5"))
if (is.na(CV_NFOLDS) || CV_NFOLDS < 3) {
  CV_NFOLDS <- 3L
}
CORRELATION_THRESHOLD <- as.numeric(Sys.getenv("CORRELATION_THRESHOLD", "0.1"))
VARIANCE_THRESHOLD <- as.numeric(Sys.getenv("VARIANCE_THRESHOLD", "0.01"))
ALPHA_GRID <- seq(0.2, 1.0, by = 0.1)
GLOBAL_SEED <- as.integer(Sys.getenv("GLOBAL_SEED", "42"))

NCORES <- as.integer(Sys.getenv("NCORES", NA))
if (is.na(NCORES)) {
  NCORES <- max(1L, parallel::detectCores() - 1L)
}

if (!dir.exists(paths$scratch)) {
  dir.create(paths$scratch, recursive = TRUE)
}

cat(strrep("=", 80), "\n")
cat("ANALYSIS: DSEN vs Standard Elastic Net\n")
cat("Mode:", if (TEST_MODE) "TEST" else "PRODUCTION", "\n")
cat("Max drugs:", if (is.infinite(MAX_DRUGS)) "all" else MAX_DRUGS, "\n")
cat("CV folds:", CV_NFOLDS, "\n")
cat("Min tissues:", MIN_TISSUES, "\n")
cat(strrep("=", 80), "\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
data <- load_main_data("Sanger")
all_drugs <- get_drug_list(data)

if (!is.infinite(MAX_DRUGS)) {
  all_drugs <- all_drugs[1:min(MAX_DRUGS, length(all_drugs))]
  cat("TEST MODE: Using", length(all_drugs), "drugs\n")
}

cat("Total drugs to process:", length(all_drugs), "\n\n")

# ==============================================================================
# CHECKPOINT SYSTEM
# ==============================================================================
checkpoint_file <- file.path(paths$scratch, "dsen_checkpoint.Rda")
results_file <- file.path(paths$scratch, "dsen_results.Rda")

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

# ==============================================================================
# PARALLEL SETUP
# ==============================================================================
cat("Setting up parallel processing...\n")
cl <- setup_parallel(NCORES)
cat("Using", NCORES, "cores\n\n")

# ==============================================================================
# DRUG PROCESSING FUNCTION
# ==============================================================================
process_drug <- function(drug_name, data, seed) {
  drug_data <- prepare_drug_data_filtered(
    data,
    drug_name,
    min_per_tissue = MIN_SAMPLES_PER_TISSUE
  )

  if (is.null(drug_data) || drug_data$n_tissues < MIN_TISSUES) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = NA,
      n_tissues = NA,
      mse_standard = NA,
      mse_dsen = NA,
      mse_reduction_pct = NA,
      alpha_standard = NA,
      alpha_dsen = NA,
      lambda_standard = NA,
      lambda_dsen = NA,
      error = NA
    ))
  }

  x <- drug_data$x
  y <- drug_data$y
  tissue <- drug_data$tissue

  # Create shared tissue-stratified folds for both models
  set.seed(seed)
  folds <- createFolds(tissue, k = CV_NFOLDS, list = FALSE)

  standard_result <- cv_tune_alpha(
    x, y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed,
    foldid = folds
  )

  dsen_result <- dsen_tune_alpha(
    x, y, tissue,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed,
    foldid = folds
  )

  mse_reduction_pct <- (standard_result$best_cvm - dsen_result$best_cvm) /
    standard_result$best_cvm * 100

  list(
    drug = drug_name,
    status = "success",
    n_samples = drug_data$n_samples,
    n_tissues = drug_data$n_tissues,
    mse_standard = standard_result$best_cvm,
    mse_dsen = dsen_result$best_cvm,
    mse_reduction_pct = mse_reduction_pct,
    alpha_standard = standard_result$alpha.min,
    alpha_dsen = dsen_result$alpha.min,
    lambda_standard = standard_result$lambda.min,
    lambda_dsen = dsen_result$lambda.min,
    error = NA
  )
}

# ==============================================================================
# MAIN COMPUTATION
# ==============================================================================
start_time <- Sys.time()

if (length(remaining_drugs) > 0) {
  batches <- split(
    remaining_drugs,
    ceiling(seq_along(remaining_drugs) / CHECKPOINT_EVERY)
  )

  for (batch_idx in seq_along(batches)) {
    batch <- batches[[batch_idx]]
    cat(sprintf("\n[%s] Processing batch %d/%d (%d drugs)\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, length(batches), length(batch)))
    cat("Drugs in batch:", paste(batch, collapse=", "), "\n\n")

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret", "Matrix"),
      .export = c("process_drug", "prepare_drug_data_filtered",
                  "cv_tune_alpha", "dsen_tune_alpha",
                  "ALPHA_GRID", "CV_NFOLDS", "GLOBAL_SEED",
                  "CORRELATION_THRESHOLD", "VARIANCE_THRESHOLD",
                  "MIN_SAMPLES_PER_TISSUE", "MIN_TISSUES"),
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
          n_samples = NA,
          n_tissues = NA,
          mse_standard = NA,
          mse_dsen = NA,
          mse_reduction_pct = NA,
          alpha_standard = NA,
          alpha_dsen = NA,
          lambda_standard = NA,
          lambda_dsen = NA,
          error = e$message
        ))
      })
    }

    all_results <- c(all_results, batch_results)
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    # Report batch completion
    batch_success <- sum(sapply(batch_results, function(r) r$status == "success"))
    batch_skipped <- sum(sapply(batch_results, function(r) r$status == "skipped_insufficient_tissues"))
    batch_errors <- sum(sapply(batch_results, function(r) r$status == "error"))
    cat(sprintf("[%s] Batch %d complete: %d success, %d skipped, %d errors\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, batch_success, batch_skipped, batch_errors))
    cat(sprintf("Overall progress: %d/%d drugs (%.1f%%)\n",
                length(completed_drugs), length(all_drugs),
                100 * length(completed_drugs) / length(all_drugs)))

    # Estimate remaining time
    elapsed_mins <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate <- length(completed_drugs) / elapsed_mins
    remaining <- length(all_drugs) - length(completed_drugs)
    eta_mins <- remaining / rate
    cat(sprintf("Rate: %.2f drugs/min, ETA: %.1f hours\n\n",
                rate, eta_mins / 60))
  }
}

# ==============================================================================
# FINALIZE RESULTS
# ==============================================================================
if (length(completed_drugs) >= length(all_drugs)) {
  cat("All drugs processed. Building results dataframe...\n")

  dsen_results <- do.call(rbind, lapply(all_results, function(res) {
    data.frame(
      drug = res$drug,
      status = res$status,
      n_samples = res$n_samples,
      n_tissues = res$n_tissues,
      mse_standard = res$mse_standard,
      mse_dsen = res$mse_dsen,
      mse_reduction_pct = res$mse_reduction_pct,
      alpha_standard = res$alpha_standard,
      alpha_dsen = res$alpha_dsen,
      lambda_standard = res$lambda_standard,
      lambda_dsen = res$lambda_dsen,
      error = res$error,
      stringsAsFactors = FALSE
    )
  }))

  save(dsen_results, file = results_file)
  cat("Results saved to:", results_file, "\n")

  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint cleaned up\n")
  }

  cat("\nSummary:\n")
  cat("  Total drugs:", nrow(dsen_results), "\n")
  cat("  Successful:", sum(dsen_results$status == "success"), "\n")
  cat("  Failed:", sum(dsen_results$status == "error"), "\n")
  cat("  Skipped:", sum(dsen_results$status == "skipped_insufficient_tissues"), "\n")
}

# ==============================================================================
# CLEANUP
# ==============================================================================
stopCluster(cl)
registerDoSEQ()

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("Elapsed time:", round(elapsed, 1), "minutes\n")
cat("Finished at", format(Sys.time()), "\n")
