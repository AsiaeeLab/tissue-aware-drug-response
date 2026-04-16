#!/usr/bin/env Rscript
# Purpose: Run cross-tissue prediction analysis with checkpointing.
# Saves to: paths$scratch/cross_tissue_results.Rda

suppressPackageStartupMessages({
  library(here)
  library(dplyr)  # For bind_rows (safer than rbind for data frames)
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

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "5"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) {
  CHECKPOINT_EVERY <- 5L
}

MIN_SAMPLES_PER_TISSUE <- as.integer(Sys.getenv("MIN_SAMPLES_PER_TISSUE", "20"))
MIN_SAMPLES <- as.integer(Sys.getenv("MIN_SAMPLES", "10"))

CV_NFOLDS <- as.integer(Sys.getenv("CV_NFOLDS", if (TEST_MODE) "3" else "5"))
if (is.na(CV_NFOLDS) || CV_NFOLDS < 3) {
  CV_NFOLDS <- 3L
}
CORRELATION_THRESHOLD <- as.numeric(Sys.getenv("CORRELATION_THRESHOLD", "0.1"))
VARIANCE_THRESHOLD <- as.numeric(Sys.getenv("VARIANCE_THRESHOLD", "0.01"))
ALPHA <- as.numeric(Sys.getenv("ALPHA", "0.5"))
GLOBAL_SEED <- as.integer(Sys.getenv("GLOBAL_SEED", "42"))

NCORES <- as.integer(Sys.getenv("NCORES", NA))
if (is.na(NCORES)) {
  NCORES <- max(1L, parallel::detectCores() - 1L)
}

if (!dir.exists(paths$scratch)) {
  dir.create(paths$scratch, recursive = TRUE)
}

cat(strrep("=", 80), "\n")
cat("ANALYSIS: Cross-Tissue Prediction\n")
cat("Mode:", if (TEST_MODE) "TEST" else "PRODUCTION", "\n")
cat("Max drugs:", if (is.infinite(MAX_DRUGS)) "all" else MAX_DRUGS, "\n")
cat("CV folds:", CV_NFOLDS, "\n")
cat("Min samples per tissue:", MIN_SAMPLES_PER_TISSUE, "\n")
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
checkpoint_file <- file.path(paths$scratch, "cross_tissue_checkpoint.Rda")
results_file <- file.path(paths$scratch, "cross_tissue_results.Rda")

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
# HELPERS
# ==============================================================================
make_status_row <- function(drug_name, status, error_msg = NA) {
  data.frame(
    drug = drug_name,
    strategy = NA,
    train_tissue = NA,
    test_tissue = NA,
    n_train = NA,
    n_test = NA,
    mse = NA,
    correlation = NA,
    r2 = NA,
    status = status,
    error = error_msg,
    stringsAsFactors = FALSE
  )
}

fit_and_eval <- function(drug_data, train_idx, test_idx,
                         train_label, test_label, strategy_label) {
  if (sum(train_idx) < MIN_SAMPLES_PER_TISSUE || sum(test_idx) < 5) return(NULL)

  x_train <- drug_data$x[train_idx, , drop = FALSE]
  y_train <- drug_data$y[train_idx]

  nfolds_use <- min(CV_NFOLDS, length(y_train))
  if (nfolds_use < 3) return(NULL)

  preproc <- preprocess_features(
    x_train, y_train,
    var_thresh = VARIANCE_THRESHOLD,
    cor_thresh = CORRELATION_THRESHOLD,
    scale_features = TRUE
  )

  if (ncol(preproc$x) < 5) return(NULL)

  set.seed(GLOBAL_SEED)
  cv_fit <- cv.glmnet(preproc$x, y_train, alpha = ALPHA, nfolds = nfolds_use)

  x_test <- drug_data$x[test_idx, , drop = FALSE]
  y_test <- drug_data$y[test_idx]
  x_test_proc <- apply_preprocessing(x_test, preproc)

  common_features <- intersect(colnames(preproc$x), colnames(x_test_proc))
  if (length(common_features) < 5) return(NULL)

  y_pred <- predict(cv_fit,
                    newx = x_test_proc[, common_features, drop = FALSE],
                    s = "lambda.min")

  mse <- mean((y_test - y_pred)^2)
  # Convert y_pred to vector (glmnet returns matrix with lambda name as column)
  y_pred_vec <- as.vector(y_pred)
  cor_val <- if (length(unique(y_test)) > 1) {
    suppressWarnings(cor(y_test, y_pred_vec))
  } else {
    NA
  }

  y_var <- var(y_test)
  r2 <- if (is.na(y_var) || y_var == 0) NA else 1 - mse / y_var

  data.frame(
    drug = drug_data$drug,
    strategy = strategy_label,
    train_tissue = train_label,
    test_tissue = test_label,
    n_train = sum(train_idx),
    n_test = sum(test_idx),
    mse = mse,
    correlation = cor_val,
    r2 = r2,
    status = "success",
    error = NA,
    stringsAsFactors = FALSE
  )
}

process_drug <- function(drug_name, data) {
  drug_data <- prepare_drug_data(data, drug_name, min_samples = MIN_SAMPLES)
  if (is.null(drug_data)) {
    return(make_status_row(drug_name, "skipped_insufficient_samples"))
  }

  drug_data$drug <- drug_name
  tissue_counts <- table(drug_data$tissue)
  valid_tissues <- names(tissue_counts)[tissue_counts >= MIN_SAMPLES_PER_TISSUE]

  if (length(valid_tissues) < 2) {
    return(make_status_row(drug_name, "skipped_insufficient_tissues"))
  }

  results <- list()
  tissue_labels <- drug_data$tissue

  # Within-tissue CV
  for (tis in valid_tissues) {
    train_idx <- tissue_labels == tis
    if (sum(train_idx) < MIN_SAMPLES_PER_TISSUE) next

    x_train <- drug_data$x[train_idx, , drop = FALSE]
    y_train <- drug_data$y[train_idx]
    nfolds_use <- min(CV_NFOLDS, length(y_train))
    if (nfolds_use < 3) next

    preproc <- preprocess_features(
      x_train, y_train,
      var_thresh = VARIANCE_THRESHOLD,
      cor_thresh = CORRELATION_THRESHOLD,
      scale_features = TRUE
    )

    if (ncol(preproc$x) < 5) next

    set.seed(GLOBAL_SEED)
    cv_fit <- cv.glmnet(preproc$x, y_train, alpha = ALPHA, nfolds = nfolds_use)

    results[[length(results) + 1]] <- data.frame(
      drug = drug_name,
      strategy = "within",
      train_tissue = tis,
      test_tissue = tis,
      n_train = sum(train_idx),
      n_test = sum(train_idx),
      mse = min(cv_fit$cvm),
      correlation = NA,
      r2 = NA,
      status = "success",
      error = NA,
      stringsAsFactors = FALSE
    )
  }

  # Cross-tissue: train on all other tissues, test on held-out tissue
  for (tis in valid_tissues) {
    test_idx <- tissue_labels == tis
    train_idx <- tissue_labels %in% setdiff(valid_tissues, tis)

    res <- fit_and_eval(
      drug_data,
      train_idx,
      test_idx,
      train_label = "all_other",
      test_label = tis,
      strategy_label = "cross"
    )

    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  }

  # Reverse: train on one tissue, test on all others
  for (tis in valid_tissues) {
    train_idx <- tissue_labels == tis
    test_idx <- tissue_labels %in% setdiff(valid_tissues, tis)

    res <- fit_and_eval(
      drug_data,
      train_idx,
      test_idx,
      train_label = tis,
      test_label = "all_other",
      strategy_label = "reverse"
    )

    if (!is.null(res)) {
      results[[length(results) + 1]] <- res
    }
  }

  if (length(results) == 0) {
    return(make_status_row(drug_name, "skipped_no_results"))
  }

  bind_rows(results)
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
    cat("Drugs in batch:", paste(batch, collapse=", "), "\n")
    cat("Three-way evaluation: within-tissue + cross-tissue + reverse\n\n")

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret", "Matrix", "dplyr"),
      .export = c("process_drug", "prepare_drug_data", "preprocess_features",
                  "apply_preprocessing", "fit_and_eval", "make_status_row",
                  "VARIANCE_THRESHOLD", "CORRELATION_THRESHOLD", "CV_NFOLDS",
                  "MIN_SAMPLES_PER_TISSUE", "MIN_SAMPLES", "ALPHA",
                  "GLOBAL_SEED"),
      .combine = c,
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      tryCatch({
        list(process_drug(drug_name, data))
      }, error = function(e) {
        list(make_status_row(drug_name, "error", e$message))
      })
    }

    all_results <- c(all_results, batch_results)
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    # Report batch completion
    batch_df <- bind_rows(batch_results)
    batch_success <- sum(batch_df$status == "success", na.rm = TRUE)
    cat(sprintf("[%s] Batch %d complete: %d result rows generated\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, nrow(batch_df)))
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
  cross_tissue_results <- bind_rows(all_results)

  save(cross_tissue_results, file = results_file)
  cat("Results saved to:", results_file, "\n")

  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint cleaned up\n")
  }

  cat("\nSummary:\n")
  cat("  Total rows:", nrow(cross_tissue_results), "\n")
  cat("  Drugs processed:", length(unique(cross_tissue_results$drug)), "\n")
}

# ==============================================================================
# CLEANUP
# ==============================================================================
stopCluster(cl)
registerDoSEQ()

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("Elapsed time:", round(elapsed, 1), "minutes\n")
cat("Finished at", format(Sys.time()), "\n")
