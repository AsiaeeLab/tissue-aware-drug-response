#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "03_dsen_utils.R"))
source(here("methods", "05_parallel_utils.R"))
source(here("methods", "comparison_common.R"))

# =============================================================================
# CONFIG
# =============================================================================
RUN_ID <- Sys.getenv("RUN_ID", "dsen_all")
FEATURE_SUBSET <- Sys.getenv("FEATURE_SUBSET", "all_features")

FEATURE_SUBSETS_FILE <- Sys.getenv(
  "FEATURE_SUBSETS_FILE",
  here("results", "comparison", "feature_subsets.rds")
)
FOLD_ASSIGNMENTS_FILE <- Sys.getenv(
  "FOLD_ASSIGNMENTS_FILE",
  here("results", "comparison", "fold_assignments.rds")
)

RESULTS_DIR <- Sys.getenv(
  "RESULTS_DIR",
  here("results", "comparison", "ablation", RUN_ID)
)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

TEST_MODE <- Sys.getenv("TEST_MODE", "FALSE") == "TRUE"

max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) {
  as.integer(max_drugs_env)
} else if (TEST_MODE) {
  2L
} else {
  Inf
}

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "3"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) CHECKPOINT_EVERY <- 3L

MIN_SAMPLES <- 50L
MIN_SAMPLES_PER_TISSUE <- 25L
MIN_TISSUES <- 3L

CV_NFOLDS <- 5L
GLOBAL_SEED <- 42L
BOOTSTRAP_NREPS <- as.integer(Sys.getenv("BOOTSTRAP_N", "200"))
if (is.na(BOOTSTRAP_NREPS) || BOOTSTRAP_NREPS < 5) BOOTSTRAP_NREPS <- 200L
SKIP_BOOTSTRAP <- Sys.getenv("SKIP_BOOTSTRAP", "TRUE") == "TRUE"

CORRELATION_THRESHOLD <- 0.1
VARIANCE_THRESHOLD <- 0.01
ALPHA_GRID <- c(0.05, 0.1, 0.15, 0.2, seq(0.3, 1.0, 0.1))
LAMRAT_GRID <- c(0.125, 0.25, 0.5, 1, 2, 4)
SCREENING <- "pooled"
TISSUE_PENALTY_MODE <- "sample_size_mild"

NCORES <- as.integer(Sys.getenv("NCORES", "10"))
if (is.na(NCORES) || NCORES < 1) NCORES <- 10L
NCORES <- min(NCORES, 10L)

checkpoint_file <- file.path(RESULTS_DIR, "dsen_ablation_checkpoint.Rda")
results_file <- file.path(RESULTS_DIR, "dsen_production_results.Rda")
bootstrap_file <- file.path(RESULTS_DIR, "dsen_bootstrap_results.Rda")
summary_csv <- file.path(RESULTS_DIR, "dsen_summary.csv")
summary_alias_csv <- file.path(RESULTS_DIR, "summary.csv")

if (!file.exists(FEATURE_SUBSETS_FILE)) {
  stop("Missing feature subsets file: ", FEATURE_SUBSETS_FILE)
}
if (!file.exists(FOLD_ASSIGNMENTS_FILE)) {
  stop("Missing fold assignment file: ", FOLD_ASSIGNMENTS_FILE)
}

feature_bundle <- readRDS(FEATURE_SUBSETS_FILE)
feature_subsets <- if (is.list(feature_bundle) && "subsets" %in% names(feature_bundle)) {
  feature_bundle$subsets
} else {
  feature_bundle
}

if (!FEATURE_SUBSET %in% names(feature_subsets)) {
  stop("FEATURE_SUBSET not found in feature_subsets.rds: ", FEATURE_SUBSET)
}

selected_feature_names <- as.character(feature_subsets[[FEATURE_SUBSET]])
fold_assignments <- readRDS(FOLD_ASSIGNMENTS_FILE)

cat(strrep("=", 80), "\n")
cat("DSEN ABLATION RUN\n")
cat(strrep("=", 80), "\n")
cat("Run ID:", RUN_ID, "\n")
cat("Feature subset:", FEATURE_SUBSET, "\n")
cat("Features requested:", length(selected_feature_names), "\n")
cat("Output dir:", RESULTS_DIR, "\n")
cat("Bootstrap:", if (SKIP_BOOTSTRAP) "SKIPPED" else paste(BOOTSTRAP_NREPS, "reps"), "\n")
cat("Cores:", NCORES, "\n")
cat(strrep("=", 80), "\n\n")

# =============================================================================
# DATA
# =============================================================================
cat("[", format(Sys.time(), "%H:%M:%S"), "] Loading data...\n")
data <- load_main_data("Sanger")

available_features <- intersect(selected_feature_names, colnames(data$predictors))
if (length(available_features) < 2) {
  stop("Selected feature subset has fewer than 2 features available in predictors.")
}

if (length(available_features) != length(selected_feature_names)) {
  warning(sprintf(
    "Feature subset '%s': %d/%d requested features available.",
    FEATURE_SUBSET, length(available_features), length(selected_feature_names)
  ))
}

data$predictors <- data$predictors[, available_features, drop = FALSE]
all_drugs <- get_drug_list(data)

if (!is.infinite(MAX_DRUGS)) {
  all_drugs <- all_drugs[1:min(MAX_DRUGS, length(all_drugs))]
}

cat("Drugs to process:", length(all_drugs), "\n")
cat("Features available in run:", ncol(data$predictors), "\n\n")

# =============================================================================
# CHECKPOINT
# =============================================================================
if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  cat("Resuming from checkpoint:\n")
  cat("  Completed:", length(completed_drugs), "\n")
  remaining_drugs <- setdiff(all_drugs, completed_drugs)
} else {
  completed_drugs <- character(0)
  all_results <- list()
  remaining_drugs <- all_drugs
  cat("Starting fresh\n")
}
cat("Remaining:", length(remaining_drugs), "\n\n")

compute_screened_feature_counts <- function(x, y, tissue, foldid) {
  folds <- sort(unique(foldid))
  en_counts <- numeric(0)
  dsen_counts <- numeric(0)

  for (k in folds) {
    train_idx <- which(foldid != k)

    en_pre <- tryCatch(
      preprocess_features(
        x[train_idx, , drop = FALSE], y[train_idx],
        var_thresh = VARIANCE_THRESHOLD,
        cor_thresh = CORRELATION_THRESHOLD,
        scale_features = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(en_pre)) {
      en_counts <- c(en_counts, en_pre$n_features)
    }

    dsen_pre <- tryCatch(
      prepare_tissue_data(
        x, y, tissue,
        train_idx = train_idx,
        var_thresh = VARIANCE_THRESHOLD,
        cor_thresh = CORRELATION_THRESHOLD,
        screening = SCREENING
      ),
      error = function(e) NULL
    )

    if (!is.null(dsen_pre)) {
      dsen_counts <- c(dsen_counts, length(dsen_pre$selected_features))
    }
  }

  list(
    mean_en = if (length(en_counts) == 0) NA_real_ else mean(en_counts),
    mean_dsen = if (length(dsen_counts) == 0) NA_real_ else mean(dsen_counts)
  )
}

process_drug <- function(drug_name, data, fold_assignments, seed) {
  fold_entry <- fold_assignments[[drug_name]]
  if (is.null(fold_entry) || !identical(fold_entry$status, "ok")) {
    return(list(
      drug = drug_name,
      status = "skipped_no_fold",
      n_samples = NA_integer_,
      n_tissues = NA_integer_,
      n_features_available = ncol(data$predictors),
      mean_features_after_screening_en = NA_real_,
      mean_features_after_screening_dsen = NA_real_,
      mse_en = NA_real_,
      mse_dsen = NA_real_,
      pct_improvement = NA_real_,
      alpha_en = NA_real_,
      alpha_dsen = NA_real_,
      lamrat_dsen = NA_real_,
      lambda_en = NA_real_,
      lambda_dsen = NA_real_,
      en_bootstrap = NULL,
      dsen_bootstrap = NULL,
      error = NA_character_
    ))
  }

  dd <- prepare_drug_data_filtered(
    data,
    drug_name,
    min_per_tissue = MIN_SAMPLES_PER_TISSUE
  )

  if (is.null(dd)) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = NA_integer_,
      n_tissues = NA_integer_,
      n_features_available = ncol(data$predictors),
      mean_features_after_screening_en = NA_real_,
      mean_features_after_screening_dsen = NA_real_,
      mse_en = NA_real_,
      mse_dsen = NA_real_,
      pct_improvement = NA_real_,
      alpha_en = NA_real_,
      alpha_dsen = NA_real_,
      lamrat_dsen = NA_real_,
      lambda_en = NA_real_,
      lambda_dsen = NA_real_,
      en_bootstrap = NULL,
      dsen_bootstrap = NULL,
      error = NA_character_
    ))
  }

  if (length(unique(dd$tissue)) < MIN_TISSUES) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = length(dd$y),
      n_tissues = length(unique(dd$tissue)),
      n_features_available = ncol(dd$x),
      mean_features_after_screening_en = NA_real_,
      mean_features_after_screening_dsen = NA_real_,
      mse_en = NA_real_,
      mse_dsen = NA_real_,
      pct_improvement = NA_real_,
      alpha_en = NA_real_,
      alpha_dsen = NA_real_,
      lamrat_dsen = NA_real_,
      lambda_en = NA_real_,
      lambda_dsen = NA_real_,
      en_bootstrap = NULL,
      dsen_bootstrap = NULL,
      error = NA_character_
    ))
  }

  x <- dd$x
  y <- dd$y
  tissue <- as.character(dd$tissue)

  if (!identical(rownames(x), fold_entry$cell_line_ids)) {
    idx <- match(fold_entry$cell_line_ids, rownames(x))
    if (any(is.na(idx))) {
      return(list(
        drug = drug_name,
        status = "error",
        n_samples = nrow(x),
        n_tissues = length(unique(tissue)),
        n_features_available = ncol(x),
        mean_features_after_screening_en = NA_real_,
        mean_features_after_screening_dsen = NA_real_,
        mse_en = NA_real_,
        mse_dsen = NA_real_,
        pct_improvement = NA_real_,
        alpha_en = NA_real_,
        alpha_dsen = NA_real_,
        lamrat_dsen = NA_real_,
        lambda_en = NA_real_,
        lambda_dsen = NA_real_,
        en_bootstrap = NULL,
        dsen_bootstrap = NULL,
        error = "Fold cell_line_ids not found in prepared matrix"
      ))
    }

    x <- x[idx, , drop = FALSE]
    y <- y[idx]
    tissue <- tissue[idx]
  }

  if (!identical(as.character(tissue), as.character(fold_entry$tissue))) {
    return(list(
      drug = drug_name,
      status = "error",
      n_samples = nrow(x),
      n_tissues = length(unique(tissue)),
      n_features_available = ncol(x),
      mean_features_after_screening_en = NA_real_,
      mean_features_after_screening_dsen = NA_real_,
      mse_en = NA_real_,
      mse_dsen = NA_real_,
      pct_improvement = NA_real_,
      alpha_en = NA_real_,
      alpha_dsen = NA_real_,
      lamrat_dsen = NA_real_,
      lambda_en = NA_real_,
      lambda_dsen = NA_real_,
      en_bootstrap = NULL,
      dsen_bootstrap = NULL,
      error = "Fold tissue vector mismatch"
    ))
  }

  foldid <- as.integer(fold_entry$foldid)

  screened_counts <- compute_screened_feature_counts(x, y, tissue, foldid)

  en_result <- cv_tune_alpha(
    x, y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = seed,
    foldid = foldid
  )

  mse_en <- en_result$best_cvm
  alpha_en <- en_result$alpha.min
  lambda_en <- en_result$lambda.min

  en_bs <- if (SKIP_BOOTSTRAP) {
    list(cnt = NULL, nboot = BOOTSTRAP_NREPS)
  } else {
    bootstrap_features(
      x, y,
      alpha = alpha_en,
      lambda = lambda_en,
      nboot = BOOTSTRAP_NREPS,
      method = "right",
      seed = seed,
      var_thresh = VARIANCE_THRESHOLD,
      cor_thresh = CORRELATION_THRESHOLD
    )
  }

  dsen_result <- tryCatch({
    dsen_tune_alpha_lamrat(
      x, y, tissue,
      alphas = ALPHA_GRID,
      lamrats = LAMRAT_GRID,
      nfolds = CV_NFOLDS,
      tissue_penalty_mode = TISSUE_PENALTY_MODE,
      screening = SCREENING,
      seed = seed,
      foldid = foldid
    )
  }, error = function(e) {
    list(best_cvm = NA_real_, error = e$message)
  })

  if (is.na(dsen_result$best_cvm)) {
    return(list(
      drug = drug_name,
      status = "dsen_error",
      n_samples = length(y),
      n_tissues = length(unique(tissue)),
      n_features_available = ncol(x),
      mean_features_after_screening_en = screened_counts$mean_en,
      mean_features_after_screening_dsen = screened_counts$mean_dsen,
      mse_en = mse_en,
      mse_dsen = NA_real_,
      pct_improvement = NA_real_,
      alpha_en = alpha_en,
      alpha_dsen = NA_real_,
      lamrat_dsen = NA_real_,
      lambda_en = lambda_en,
      lambda_dsen = NA_real_,
      en_bootstrap = en_bs,
      dsen_bootstrap = NULL,
      error = dsen_result$error
    ))
  }

  mse_dsen <- dsen_result$best_cvm
  alpha_dsen <- dsen_result$alpha.min
  lamrat_dsen <- dsen_result$lamrat.min
  lambda_dsen <- dsen_result$lambda.min

  dsen_bs <- if (SKIP_BOOTSTRAP) {
    list(cnt = NULL, tissue_cnt = NULL, nboot = BOOTSTRAP_NREPS)
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
      list(cnt = NULL, tissue_cnt = NULL, nboot = BOOTSTRAP_NREPS, error = e$message)
    })
  }

  pct_improvement <- (mse_en - mse_dsen) / mse_en * 100

  list(
    drug = drug_name,
    status = "success",
    n_samples = length(y),
    n_tissues = length(unique(tissue)),
    n_features_available = ncol(x),
    mean_features_after_screening_en = screened_counts$mean_en,
    mean_features_after_screening_dsen = screened_counts$mean_dsen,
    mse_en = mse_en,
    mse_dsen = mse_dsen,
    pct_improvement = pct_improvement,
    alpha_en = alpha_en,
    alpha_dsen = alpha_dsen,
    lamrat_dsen = lamrat_dsen,
    lambda_en = lambda_en,
    lambda_dsen = lambda_dsen,
    en_bootstrap = en_bs,
    dsen_bootstrap = dsen_bs,
    error = NA_character_
  )
}

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
    cat(sprintf("[%s] Processing batch %d/%d (%d drugs)\n",
                format(Sys.time(), "%H:%M:%S"),
                batch_idx,
                length(batches),
                length(batch)))

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret", "Matrix"),
      .export = c(
        "process_drug",
        "compute_screened_feature_counts",
        "prepare_drug_data_filtered",
        "preprocess_features",
        "as_numeric_matrix",
        "apply_preprocessing",
        "cv_tune_alpha",
        "bootstrap_features",
        "dsen_tune_alpha_lamrat",
        "dsen_bootstrap",
        "dsen_cv",
        "compute_dsen_penalties",
        "build_dsen_matrix",
        "prepare_tissue_data",
        "dsen_predict",
        "ALPHA_GRID",
        "LAMRAT_GRID",
        "CV_NFOLDS",
        "VARIANCE_THRESHOLD",
        "CORRELATION_THRESHOLD",
        "BOOTSTRAP_NREPS",
        "TISSUE_PENALTY_MODE",
        "SCREENING",
        "MIN_SAMPLES_PER_TISSUE",
        "MIN_TISSUES",
        "SKIP_BOOTSTRAP",
        "GLOBAL_SEED",
        "data",
        "fold_assignments"
      ),
      .combine = c,
      .errorhandling = "pass",
      .verbose = FALSE
    ) %dopar% {
      tryCatch({
        list(process_drug(drug_name, data, fold_assignments, seed = GLOBAL_SEED))
      }, error = function(e) {
        list(list(
          drug = drug_name,
          status = "error",
          n_samples = NA_integer_,
          n_tissues = NA_integer_,
          n_features_available = NA_integer_,
          mean_features_after_screening_en = NA_real_,
          mean_features_after_screening_dsen = NA_real_,
          mse_en = NA_real_,
          mse_dsen = NA_real_,
          pct_improvement = NA_real_,
          alpha_en = NA_real_,
          alpha_dsen = NA_real_,
          lamrat_dsen = NA_real_,
          lambda_en = NA_real_,
          lambda_dsen = NA_real_,
          en_bootstrap = NULL,
          dsen_bootstrap = NULL,
          error = e$message
        ))
      })
    }

    for (res in batch_results) {
      all_results[[res$drug]] <- res
    }
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    batch_success <- sum(sapply(batch_results, function(r) identical(r$status, "success")))
    batch_skipped <- sum(sapply(batch_results, function(r) grepl("^skipped_", r$status)))
    batch_errors <- sum(sapply(batch_results, function(r) r$status %in% c("error", "dsen_error")))

    cat(sprintf("[%s] Batch done: %d success | %d skipped | %d errors\n",
                format(Sys.time(), "%H:%M:%S"),
                batch_success,
                batch_skipped,
                batch_errors))
    cat(sprintf("Progress: %d/%d (%.1f%%)\n\n",
                length(completed_drugs),
                length(all_drugs),
                100 * length(completed_drugs) / length(all_drugs)))

    append_progress(
      phase = sprintf("Phase 1 - %s", RUN_ID),
      status = "in_progress",
      details = c(
        sprintf("Batch %d/%d", batch_idx, length(batches)),
        sprintf("Progress: %d/%d", length(completed_drugs), length(all_drugs)),
        sprintf("Batch success/skipped/errors: %d/%d/%d", batch_success, batch_skipped, batch_errors)
      )
    )
  }

  stopCluster(cl)
  registerDoSEQ()
}

summary_df <- do.call(rbind, lapply(all_results, function(res) {
  data.frame(
    drug = res$drug,
    status = res$status,
    run_id = RUN_ID,
    feature_subset = FEATURE_SUBSET,
    n_samples = res$n_samples,
    n_tissues = res$n_tissues,
    n_features_available = res$n_features_available,
    mean_features_after_screening_en = res$mean_features_after_screening_en,
    mean_features_after_screening_dsen = res$mean_features_after_screening_dsen,
    mse_en = res$mse_en,
    mse_dsen = res$mse_dsen,
    pct_improvement = res$pct_improvement,
    dsen_wins = ifelse(is.finite(res$pct_improvement), res$pct_improvement > 0, NA),
    alpha_en = res$alpha_en,
    alpha_dsen = res$alpha_dsen,
    lamrat_dsen = res$lamrat_dsen,
    lambda_en = res$lambda_en,
    lambda_dsen = res$lambda_dsen,
    error = res$error,
    stringsAsFactors = FALSE
  )
}))

# Keep canonical drug order for readability
summary_df <- summary_df[match(all_drugs, summary_df$drug), , drop = FALSE]

en_bootstrap <- lapply(all_results, function(res) res$en_bootstrap)
dsen_bootstrap <- lapply(all_results, function(res) res$dsen_bootstrap)

config <- list(
  run_id = RUN_ID,
  feature_subset = FEATURE_SUBSET,
  n_features_requested = length(selected_feature_names),
  n_features_available = length(available_features),
  test_mode = TEST_MODE,
  cv_nfolds = CV_NFOLDS,
  bootstrap_nreps = BOOTSTRAP_NREPS,
  skip_bootstrap = SKIP_BOOTSTRAP,
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

production_results <- list(
  summary_df = summary_df,
  en_bootstrap = en_bootstrap,
  dsen_bootstrap = dsen_bootstrap,
  config = config
)

save(production_results, file = results_file)
save(all_results, file = bootstrap_file)
write.csv(summary_df, summary_csv, row.names = FALSE)
write.csv(summary_df, summary_alias_csv, row.names = FALSE)

if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

success_df <- summary_df %>% filter(status == "success")
win_rate <- if (nrow(success_df) == 0) NA_real_ else mean(success_df$pct_improvement > 0, na.rm = TRUE) * 100

cat("\n", strrep("=", 80), "\n", sep = "")
cat("Ablation run complete\n")
cat(strrep("=", 80), "\n")
cat("Run ID:", RUN_ID, "\n")
cat("Feature subset:", FEATURE_SUBSET, "\n")
cat("Successful drugs:", nrow(success_df), "/", nrow(summary_df), "\n")
cat("DSEN win rate:", sprintf("%.2f%%", win_rate), "\n")
cat("Mean improvement:", sprintf("%.3f%%", mean(success_df$pct_improvement, na.rm = TRUE)), "\n")
cat("Median improvement:", sprintf("%.3f%%", median(success_df$pct_improvement, na.rm = TRUE)), "\n")
cat("Elapsed minutes:", round(elapsed, 1), "\n")
cat("Saved:", summary_csv, "\n")
cat(strrep("=", 80), "\n")

append_progress(
  phase = sprintf("Phase 1 - %s", RUN_ID),
  status = "completed",
  runtime = sprintf("%.1f minutes", elapsed),
  details = c(
    sprintf("Feature subset: %s", FEATURE_SUBSET),
    sprintf("Features available: %d", length(available_features)),
    sprintf("Successful drugs: %d/%d", nrow(success_df), nrow(summary_df)),
    sprintf("DSEN win rate: %.2f%%", win_rate),
    sprintf("Mean improvement: %.3f%%", mean(success_df$pct_improvement, na.rm = TRUE)),
    sprintf("Median improvement: %.3f%%", median(success_df$pct_improvement, na.rm = TRUE)),
    sprintf("Bootstrap skipped: %s", ifelse(SKIP_BOOTSTRAP, "TRUE", "FALSE"))
  )
)
