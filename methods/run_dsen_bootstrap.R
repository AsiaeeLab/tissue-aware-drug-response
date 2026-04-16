#!/usr/bin/env Rscript
# Purpose: Run bootstrap feature selection for EN and DSEN
# Uses optimal parameters from Stage 1 results
# Saves raw bootstrap counts for post-hoc threshold analysis

suppressPackageStartupMessages({
  library(here)
  library(glmnet)
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

# Max drugs to process (for testing)
max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) {
  as.integer(max_drugs_env)
} else {
  Inf
}

BOOTSTRAP_NREPS <- as.integer(Sys.getenv("BOOTSTRAP_N", if (TEST_MODE) "10" else "200"))
if (is.na(BOOTSTRAP_NREPS) || BOOTSTRAP_NREPS < 5) BOOTSTRAP_NREPS <- 200L

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "20"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) CHECKPOINT_EVERY <- 20L

MIN_SAMPLES <- as.integer(Sys.getenv("MIN_SAMPLES", "50"))
MIN_SAMPLES_PER_TISSUE <- as.integer(Sys.getenv("MIN_SAMPLES_PER_TISSUE", "15"))
if (is.na(MIN_SAMPLES_PER_TISSUE) || MIN_SAMPLES_PER_TISSUE < 1) {
  MIN_SAMPLES_PER_TISSUE <- 15L
}
MIN_TISSUES <- as.integer(Sys.getenv("MIN_TISSUES", "3"))
if (is.na(MIN_TISSUES) || MIN_TISSUES < 2) {
  MIN_TISSUES <- 3L
}
VARIANCE_THRESHOLD <- as.numeric(Sys.getenv("VARIANCE_THRESHOLD", "0.01"))
CORRELATION_THRESHOLD <- as.numeric(Sys.getenv("CORRELATION_THRESHOLD", "0.1"))
TISSUE_PENALTY_MODE <- Sys.getenv("TISSUE_PENALTY_MODE", "sample_size_mild")
GLOBAL_SEED <- as.integer(Sys.getenv("GLOBAL_SEED", "42"))

NCORES <- as.integer(Sys.getenv("NCORES", NA))
if (is.na(NCORES)) {
  NCORES <- max(1L, parallel::detectCores() - 1L)
}

# Paths
stage1_csv <- Sys.getenv(
  "STAGE1_SUMMARY_CSV",
  here("results", "dsen_results", "dsen_summary.csv")
)
results_file <- Sys.getenv(
  "BOOTSTRAP_RESULTS_FILE",
  here("results", "dsen_results", "dsen_bootstrap_results.Rda")
)
checkpoint_file <- Sys.getenv(
  "BOOTSTRAP_CHECKPOINT_FILE",
  here("results", "dsen_results", "dsen_bootstrap_checkpoint.Rda")
)

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 80), "\n")
cat("DSEN BOOTSTRAP: Feature Selection Analysis\n")
cat(strrep("=", 80), "\n")
cat("Mode:", if (TEST_MODE) "TEST" else "PRODUCTION", "\n")
cat("Max drugs:", if (is.infinite(MAX_DRUGS)) "all" else MAX_DRUGS, "\n")
cat("Bootstrap reps:", BOOTSTRAP_NREPS, "\n")
cat("Checkpoint every:", CHECKPOINT_EVERY, "drugs\n")
cat("Min samples per tissue:", MIN_SAMPLES_PER_TISSUE, "\n")
cat("Min tissues:", MIN_TISSUES, "\n")
cat("Stage 1 summary:", stage1_csv, "\n")
cat("Results file:", results_file, "\n")
cat("Checkpoint file:", checkpoint_file, "\n")
cat("Cores:", NCORES, "\n")
cat(strrep("=", 80), "\n\n")

# ==============================================================================
# LOAD STAGE 1 RESULTS
# ==============================================================================
cat("Loading Stage 1 results...\n")
if (!file.exists(stage1_csv)) {
  stop("Stage 1 results not found: ", stage1_csv)
}

stage1_df <- read.csv(stage1_csv)
cat("Found", nrow(stage1_df), "drugs from Stage 1\n")

if (!"status" %in% names(stage1_df)) {
  stop("Stage 1 CSV must include a 'status' column.")
}

# Filter successful drugs
stage1_df <- stage1_df[stage1_df$status == "success", ]
cat("Successful drugs:", nrow(stage1_df), "\n")

# Harmonize Stage 1 column names across runner variants
if (!"alpha_en" %in% names(stage1_df)) {
  if ("alpha_standard" %in% names(stage1_df)) {
    stage1_df$alpha_en <- stage1_df$alpha_standard
  } else {
    stop("Stage 1 CSV missing alpha for EN ('alpha_en' or 'alpha_standard').")
  }
}

if (!"alpha_dsen" %in% names(stage1_df)) {
  stop("Stage 1 CSV missing 'alpha_dsen'.")
}

if (!"lamrat_dsen" %in% names(stage1_df)) {
  stage1_df$lamrat_dsen <- NA_real_
}

if (!"lambda_en" %in% names(stage1_df)) {
  if ("lambda_standard" %in% names(stage1_df)) {
    stage1_df$lambda_en <- stage1_df$lambda_standard
  } else {
    stage1_df$lambda_en <- NA_real_
  }
}

if (!"lambda_dsen" %in% names(stage1_df)) {
  stage1_df$lambda_dsen <- NA_real_
}

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("\nLoading data...\n")
data <- load_main_data("Sanger")

# ==============================================================================
# CHECK FOR EXISTING CHECKPOINT
# ==============================================================================
all_drugs <- stage1_df$drug
if (!is.infinite(MAX_DRUGS) && MAX_DRUGS < length(all_drugs)) {
  all_drugs <- all_drugs[1:MAX_DRUGS]
  cat("Limited to first", MAX_DRUGS, "drugs for testing\n")
}

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  remaining_drugs <- setdiff(all_drugs, completed_drugs)
  cat("Resuming from checkpoint:", length(completed_drugs), "completed\n")
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
# BOOTSTRAP FUNCTION
# ==============================================================================
process_drug_bootstrap <- function(drug_name, data, stage1_params, seed) {
  # Get Stage 1 parameters for this drug
  params <- stage1_params[stage1_params$drug == drug_name, ]
  if (nrow(params) == 0) {
    return(list(
      drug = drug_name,
      status = "error",
      error = "No Stage 1 parameters found"
    ))
  }

  # Prepare data with the same per-tissue filtering used in MSE analysis
  drug_data <- prepare_drug_data_filtered(
    data,
    drug_name,
    min_per_tissue = MIN_SAMPLES_PER_TISSUE
  )
  if (is.null(drug_data) || drug_data$n_tissues < MIN_TISSUES) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues"
    ))
  }

  x <- drug_data$x
  y <- drug_data$y
  tissue <- drug_data$tissue

  # Get optimal parameters from Stage 1
  alpha_en <- as.numeric(params$alpha_en[1])
  alpha_dsen <- as.numeric(params$alpha_dsen[1])
  lamrat_dsen <- as.numeric(params$lamrat_dsen[1])
  if (!is.finite(lamrat_dsen) || lamrat_dsen <= 0) {
    lamrat_dsen <- NULL
  }

  lambda_en <- as.numeric(params$lambda_en[1])
  if (!is.finite(lambda_en) || lambda_en <= 0) {
    lambda_en <- NA_real_
  }

  lambda_dsen <- as.numeric(params$lambda_dsen[1])
  if (!is.finite(lambda_dsen) || lambda_dsen <= 0) {
    lambda_dsen <- NA_real_
  }

  # Re-fit only when Stage 1 summary lacks lambda values

  # EN: fallback CV to get lambda when Stage 1 does not provide one.
  if (!is.finite(lambda_en)) {
    set.seed(seed)
    en_fit <- tryCatch({
      my_cv_glmnet(
        x, y,
        alpha = alpha_en,
        nfolds = 5,
        var_thresh = VARIANCE_THRESHOLD,
        cor_thresh = CORRELATION_THRESHOLD,
        seed = seed
      )
    }, error = function(e) NULL)

    if (is.null(en_fit)) {
      return(list(
        drug = drug_name,
        status = "error",
        error = "EN CV fit failed"
      ))
    }

    lambda_en <- en_fit$lambda.min
  }

  # EN Bootstrap
  en_bs <- tryCatch({
    bootstrap_features(
      x, y,
      alpha = alpha_en,
      lambda = lambda_en,
      nboot = BOOTSTRAP_NREPS,
      method = "right",
      seed = seed
    )
  }, error = function(e) {
    list(cnt = NULL, error = e$message)
  })

  # DSEN: Get lambda via quick CV
  if (!is.finite(lambda_dsen)) {
    dsen_fit <- tryCatch({
      dsen_cv(
        x, y, tissue,
        alpha = alpha_dsen,
        lamrat = lamrat_dsen,
        nfolds = 5,
        tissue_penalty_mode = TISSUE_PENALTY_MODE,
        seed = seed
      )
    }, error = function(e) NULL)

    if (is.null(dsen_fit)) {
      lambda_dsen <- lambda_en  # Fallback
    } else {
      lambda_dsen <- dsen_fit$lambda.min
    }
  }

  # DSEN Bootstrap
  dsen_bs <- tryCatch({
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
    list(cnt = NULL, error = e$message)
  })

  list(
    drug = drug_name,
    status = "success",
    n_samples = length(y),
    n_features = ncol(x),
    n_tissues = length(unique(tissue)),
    alpha_en = alpha_en,
    alpha_dsen = alpha_dsen,
    lamrat_dsen = if (is.null(lamrat_dsen)) NA_real_ else lamrat_dsen,
    lambda_en = lambda_en,
    lambda_dsen = lambda_dsen,
    en_bootstrap = en_bs,
    dsen_bootstrap = dsen_bs
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

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret", "Matrix"),
      .export = c(
        # Main processing function
        "process_drug_bootstrap",
        # From 01_data_utils.R
        "prepare_drug_data_filtered", "preprocess_features", "apply_preprocessing", "as_numeric_matrix",
        # From 02_cv_utils.R
        "bootstrap_features", "my_cv_glmnet",
        # From 03_dsen_utils.R
        "dsen_cv", "dsen_bootstrap", "compute_dsen_penalties",
        "build_dsen_matrix", "prepare_tissue_data", "dsen_predict",
        # Constants
        "VARIANCE_THRESHOLD", "CORRELATION_THRESHOLD",
        "BOOTSTRAP_NREPS", "TISSUE_PENALTY_MODE",
        "MIN_SAMPLES_PER_TISSUE", "MIN_TISSUES",
        "GLOBAL_SEED", "stage1_df", "data"
      ),
      .combine = c,
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      tryCatch({
        list(process_drug_bootstrap(drug_name, data, stage1_df, seed = GLOBAL_SEED))
      }, error = function(e) {
        list(list(
          drug = drug_name,
          status = "error",
          error = e$message
        ))
      })
    }

    # Store results
    for (i in seq_along(batch_results)) {
      res <- batch_results[[i]]
      all_results[[res$drug]] <- res
    }
    completed_drugs <- c(completed_drugs, batch)

    # Save checkpoint
    save(all_results, completed_drugs, file = checkpoint_file)

    # Report progress
    batch_success <- sum(sapply(batch_results, function(r) r$status == "success"))
    cat(sprintf("[%s] Batch complete: %d success\n",
                format(Sys.time(), "%H:%M:%S"), batch_success))
    cat(sprintf("Overall: %d/%d drugs (%.1f%%)\n",
                length(completed_drugs), length(all_drugs),
                100 * length(completed_drugs) / length(all_drugs)))

    # ETA
    elapsed_mins <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate <- length(completed_drugs) / elapsed_mins
    remaining <- length(all_drugs) - length(completed_drugs)
    eta_mins <- remaining / rate
    cat(sprintf("Rate: %.2f drugs/min, ETA: %.1f min\n\n", rate, eta_mins))
  }
}

# ==============================================================================
# FINALIZE
# ==============================================================================
if (length(completed_drugs) >= length(all_drugs)) {
  save(all_results, file = results_file)
  cat("\nResults saved to:", results_file, "\n")

  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint cleaned up\n")
  }

  # Summary
  statuses <- sapply(all_results, function(r) r$status)
  cat("\nSummary:\n")
  cat("  Total:", length(all_results), "\n")
  cat("  Success:", sum(statuses == "success"), "\n")
  cat("  Errors:", sum(statuses == "error"), "\n")
}

# ==============================================================================
# CLEANUP
# ==============================================================================
stopCluster(cl)
registerDoSEQ()

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat("\nElapsed time:", round(elapsed, 1), "minutes\n")
cat("Finished at", format(Sys.time()), "\n")
