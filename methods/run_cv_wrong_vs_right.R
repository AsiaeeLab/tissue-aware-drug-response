#!/usr/bin/env Rscript
# Purpose: Run wrong vs right CV analysis with checkpointing.
# Saves to: paths$scratch/cv_wrong_vs_right_results.Rda

suppressPackageStartupMessages({
  library(here)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
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

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "10"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) {
  CHECKPOINT_EVERY <- 10L
}

MIN_SAMPLES <- as.integer(Sys.getenv("MIN_SAMPLES", "50"))

CV_NFOLDS <- as.integer(Sys.getenv("CV_NFOLDS", if (TEST_MODE) "3" else "5"))
if (is.na(CV_NFOLDS) || CV_NFOLDS < 3) {
  CV_NFOLDS <- 3L
}
BOOTSTRAP_NREPS <- as.integer(Sys.getenv(
  "BOOTSTRAP_NREPS",
  if (TEST_MODE) "5" else "200"
))

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
cat("ANALYSIS: CV Wrong vs Right\n")
cat("Mode:", if (TEST_MODE) "TEST" else "PRODUCTION", "\n")
cat("Max drugs:", if (is.infinite(MAX_DRUGS)) "all" else MAX_DRUGS, "\n")
cat("CV folds:", CV_NFOLDS, "\n")
cat("Bootstrap reps:", BOOTSTRAP_NREPS, "\n")
cat("Checkpoint every:", CHECKPOINT_EVERY, "drugs\n")
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
checkpoint_file <- file.path(paths$scratch, "cv_wrong_vs_right_checkpoint.Rda")
results_file <- file.path(paths$scratch, "cv_wrong_vs_right_results.Rda")

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
  drug_data <- prepare_drug_data(data, drug_name, min_samples = MIN_SAMPLES)

  if (is.null(drug_data)) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_samples",
      n_samples = NA,
      n_features_total = NA,
      mse_wrong = NA,
      mse_right = NA,
      inflation_pct = NA,
      n_features_wrong = NA,
      n_features_right = NA,
      jaccard = NA,
      alpha_wrong = NA,
      alpha_right = NA,
      lambda_wrong = NA,
      lambda_right = NA,
      error = NA
    ))
  }

  x <- drug_data$x
  y <- drug_data$y

  wrong_result <- wrong_cv_tune_alpha(
    x, y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed
  )

  right_result <- cv_tune_alpha(
    x, y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed
  )

  bs_wrong <- bootstrap_features(
    x, y,
    alpha = wrong_result$alpha.min,
    lambda = wrong_result$lambda.min,
    nboot = BOOTSTRAP_NREPS,
    method = "wrong",
    seed = seed
  )

  bs_right <- bootstrap_features(
    x, y,
    alpha = right_result$alpha.min,
    lambda = right_result$lambda.min,
    nboot = BOOTSTRAP_NREPS,
    method = "right",
    seed = seed
  )

  comparison <- compare_cv_results(
    wrong_result, right_result,
    bs_wrong, bs_right,
    threshold = 0.8
  )

  list(
    drug = drug_name,
    status = "success",
    n_samples = length(y),
    n_features_total = ncol(x),
    mse_wrong = comparison$mse_wrong,
    mse_right = comparison$mse_right,
    inflation_pct = comparison$mse_pct_increase,
    n_features_wrong = comparison$n_features_wrong,
    n_features_right = comparison$n_features_right,
    jaccard = comparison$jaccard,
    alpha_wrong = comparison$alpha_wrong,
    alpha_right = comparison$alpha_right,
    lambda_wrong = comparison$lambda_wrong,
    lambda_right = comparison$lambda_right,
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
      .export = c("process_drug", "prepare_drug_data", "wrong_cv_tune_alpha",
                  "cv_tune_alpha", "bootstrap_features", "compare_cv_results",
                  "ALPHA_GRID", "CV_NFOLDS", "BOOTSTRAP_NREPS",
                  "CORRELATION_THRESHOLD", "VARIANCE_THRESHOLD",
                  "GLOBAL_SEED", "MIN_SAMPLES"),
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
          n_features_total = NA,
          mse_wrong = NA,
          mse_right = NA,
          inflation_pct = NA,
          n_features_wrong = NA,
          n_features_right = NA,
          jaccard = NA,
          alpha_wrong = NA,
          alpha_right = NA,
          lambda_wrong = NA,
          lambda_right = NA,
          error = e$message
        ))
      })
    }

    all_results <- c(all_results, batch_results)
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    # Report batch completion
    batch_success <- sum(sapply(batch_results, function(r) r$status == "success"))
    batch_errors <- sum(sapply(batch_results, function(r) r$status == "error"))
    cat(sprintf("[%s] Batch %d complete: %d success, %d errors\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, batch_success, batch_errors))
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

  cv_results <- do.call(rbind, lapply(all_results, function(res) {
    data.frame(
      drug = res$drug,
      status = res$status,
      n_samples = res$n_samples,
      n_features_total = res$n_features_total,
      n_features_wrong = res$n_features_wrong,
      n_features_right = res$n_features_right,
      mse_wrong = res$mse_wrong,
      mse_right = res$mse_right,
      inflation_pct = res$inflation_pct,
      jaccard = res$jaccard,
      alpha_wrong = res$alpha_wrong,
      alpha_right = res$alpha_right,
      lambda_wrong = res$lambda_wrong,
      lambda_right = res$lambda_right,
      error = res$error,
      stringsAsFactors = FALSE
    )
  }))

  save(cv_results, file = results_file)
  cat("Results saved to:", results_file, "\n")

  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint cleaned up\n")
  }

  cat("\nSummary:\n")
  cat("  Total drugs:", nrow(cv_results), "\n")
  cat("  Successful:", sum(cv_results$status == "success"), "\n")
  cat("  Failed:", sum(cv_results$status == "error"), "\n")
  cat("  Skipped:", sum(cv_results$status == "skipped_insufficient_samples"), "\n")
}

# ==============================================================================
# CLEANUP
# ==============================================================================
stopCluster(cl)
registerDoSEQ()

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("Elapsed time:", round(elapsed, 1), "minutes\n")
cat("Finished at", format(Sys.time()), "\n")
