#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "04_tglasso_utils.R"))
source(here("methods", "05_parallel_utils.R"))
source(here("methods", "comparison_common.R"))

# =============================================================================
# CONFIG
# =============================================================================
GLOBAL_SEED <- 42L
CV_NFOLDS <- 5L
MIN_SAMPLES_PER_TISSUE <- 25L
MIN_TISSUES <- 3L
CORRELATION_THRESHOLD <- 0.1
VARIANCE_THRESHOLD <- 0.01

SKIP_BOOTSTRAP <- Sys.getenv("SKIP_BOOTSTRAP", "TRUE") == "TRUE"
BOOTSTRAP_NREPS <- as.integer(Sys.getenv("BOOTSTRAP_N", "200"))
if (is.na(BOOTSTRAP_NREPS) || BOOTSTRAP_NREPS < 5) BOOTSTRAP_NREPS <- 200L

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "3"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) CHECKPOINT_EVERY <- 3L

NCORES <- as.integer(Sys.getenv("NCORES", "10"))
if (is.na(NCORES) || NCORES < 1) NCORES <- 10L
NCORES <- min(NCORES, 10L)

max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) as.integer(max_drugs_env) else Inf

FOLD_ASSIGNMENTS_FILE <- Sys.getenv(
  "FOLD_ASSIGNMENTS_FILE",
  here("results", "comparison", "fold_assignments.rds")
)

FEATURE_SUBSETS_FILE <- Sys.getenv(
  "FEATURE_SUBSETS_FILE",
  here("results", "comparison", "feature_subsets.rds")
)
FEATURE_SUBSET <- Sys.getenv("FEATURE_SUBSET", "all_features")

DSEN_ALL_SUMMARY <- Sys.getenv(
  "DSEN_ALL_SUMMARY",
  here("results", "comparison", "ablation", "dsen_all", "dsen_summary.csv")
)

OUT_DIR <- Sys.getenv("OUT_DIR", here("results", "comparison", "tglasso"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results_rds <- file.path(OUT_DIR, "tglasso_results.rds")
summary_csv <- file.path(OUT_DIR, "tglasso_summary.csv")
bootstrap_rds <- file.path(OUT_DIR, "tglasso_bootstrap_features.rds")
checkpoint_file <- file.path(OUT_DIR, "tglasso_checkpoint.Rda")

if (!file.exists(FOLD_ASSIGNMENTS_FILE)) {
  stop("Missing fold assignment file: ", FOLD_ASSIGNMENTS_FILE)
}
if (!file.exists(FEATURE_SUBSETS_FILE)) {
  stop("Missing feature subset file: ", FEATURE_SUBSETS_FILE)
}

feature_bundle <- readRDS(FEATURE_SUBSETS_FILE)
feature_subsets <- if (is.list(feature_bundle) && "subsets" %in% names(feature_bundle)) {
  feature_bundle$subsets
} else {
  feature_bundle
}
if (!FEATURE_SUBSET %in% names(feature_subsets)) {
  stop("FEATURE_SUBSET not found: ", FEATURE_SUBSET)
}

selected_feature_names <- feature_subsets[[FEATURE_SUBSET]]
fold_assignments <- readRDS(FOLD_ASSIGNMENTS_FILE)

cat(strrep("=", 80), "\n")
cat("TG-LASSO RUNNER\n")
cat(strrep("=", 80), "\n")
cat("Feature subset:", FEATURE_SUBSET, "\n")
cat("Bootstrap:", if (SKIP_BOOTSTRAP) "SKIPPED" else paste(BOOTSTRAP_NREPS, "reps"), "\n")
cat("Output dir:", OUT_DIR, "\n")
cat("Cores:", NCORES, "\n")
cat(strrep("=", 80), "\n\n")

cat("Loading data...\n")
data <- load_main_data("Sanger")
available_features <- intersect(selected_feature_names, colnames(data$predictors))
data$predictors <- data$predictors[, available_features, drop = FALSE]

all_drugs <- get_drug_list(data)
if (!is.infinite(MAX_DRUGS)) {
  all_drugs <- all_drugs[1:min(MAX_DRUGS, length(all_drugs))]
}

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  remaining_drugs <- setdiff(all_drugs, completed_drugs)
  cat("Resuming from checkpoint. Completed:", length(completed_drugs), "\n")
} else {
  completed_drugs <- character(0)
  all_results <- list()
  remaining_drugs <- all_drugs
  cat("Starting fresh\n")
}
cat("Remaining:", length(remaining_drugs), "\n\n")

process_drug_tglasso <- function(drug_name, data, fold_assignments) {
  fold_entry <- fold_assignments[[drug_name]]
  if (is.null(fold_entry) || !identical(fold_entry$status, "ok")) {
    return(list(
      drug = drug_name,
      status = "skipped_no_fold",
      n_samples = NA_integer_,
      n_tissues = NA_integer_,
      mse_tglasso = NA_real_,
      lambda_by_fold_tissue = NULL,
      selected_features = NULL,
      bootstrap = NULL,
      error = NA_character_
    ))
  }

  dd <- prepare_drug_data_filtered(data, drug_name, min_per_tissue = MIN_SAMPLES_PER_TISSUE)
  if (is.null(dd) || length(unique(dd$tissue)) < MIN_TISSUES) {
    return(list(
      drug = drug_name,
      status = "skipped_insufficient_tissues",
      n_samples = if (is.null(dd)) NA_integer_ else length(dd$y),
      n_tissues = if (is.null(dd)) NA_integer_ else length(unique(dd$tissue)),
      mse_tglasso = NA_real_,
      lambda_by_fold_tissue = NULL,
      selected_features = NULL,
      bootstrap = NULL,
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
        mse_tglasso = NA_real_,
        lambda_by_fold_tissue = NULL,
        selected_features = NULL,
        bootstrap = NULL,
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
      mse_tglasso = NA_real_,
      lambda_by_fold_tissue = NULL,
      selected_features = NULL,
      bootstrap = NULL,
      error = "Fold tissue vector mismatch"
    ))
  }

  foldid <- as.integer(fold_entry$foldid)

  cv_res <- tryCatch({
    tglasso_cv(
      x, y, tissue,
      nfolds = CV_NFOLDS,
      foldid = foldid,
      seed = GLOBAL_SEED,
      cor_threshold = CORRELATION_THRESHOLD,
      var_threshold = VARIANCE_THRESHOLD
    )
  }, error = function(e) {
    list(error = e$message)
  })

  if (!is.null(cv_res$error)) {
    return(list(
      drug = drug_name,
      status = "error",
      n_samples = nrow(x),
      n_tissues = length(unique(tissue)),
      mse_tglasso = NA_real_,
      lambda_by_fold_tissue = NULL,
      selected_features = NULL,
      bootstrap = NULL,
      error = cv_res$error
    ))
  }

  bs <- if (SKIP_BOOTSTRAP) {
    NULL
  } else {
    tryCatch(
      tglasso_bootstrap_features(
        x, y, tissue,
        nboot = BOOTSTRAP_NREPS,
        seed = GLOBAL_SEED,
        var_threshold = VARIANCE_THRESHOLD,
        cor_threshold = CORRELATION_THRESHOLD
      ),
      error = function(e) list(error = e$message)
    )
  }

  list(
    drug = drug_name,
    status = "success",
    n_samples = nrow(x),
    n_tissues = length(unique(tissue)),
    mse_tglasso = cv_res$cvm,
    lambda_by_fold_tissue = cv_res$lambda_by_fold_tissue,
    selected_features = cv_res$selected_features,
    bootstrap = bs,
    error = NA_character_
  )
}

start_time <- Sys.time()

if (length(remaining_drugs) > 0) {
  cl <- setup_parallel(NCORES)

  batches <- split(
    remaining_drugs,
    ceiling(seq_along(remaining_drugs) / CHECKPOINT_EVERY)
  )

  for (batch_idx in seq_along(batches)) {
    batch <- batches[[batch_idx]]
    cat(sprintf("[%s] Batch %d/%d (%d drugs)\n",
                format(Sys.time(), "%H:%M:%S"),
                batch_idx,
                length(batches),
                length(batch)))

    batch_results <- foreach(
      drug_name = batch,
      .packages = c("glmnet", "caret"),
      .export = c(
        "process_drug_tglasso",
        "prepare_drug_data_filtered",
        "tglasso_cv",
        "tglasso_bootstrap_features",
        "fit_tglasso_models_for_training",
        "fit_tglasso_tissue_model",
        "predict_tglasso_model",
        "align_to_features",
        "preprocess_features",
        "apply_preprocessing",
        "as_numeric_matrix",
        "MIN_SAMPLES_PER_TISSUE",
        "MIN_TISSUES",
        "CV_NFOLDS",
        "GLOBAL_SEED",
        "CORRELATION_THRESHOLD",
        "VARIANCE_THRESHOLD",
        "SKIP_BOOTSTRAP",
        "BOOTSTRAP_NREPS",
        "data",
        "fold_assignments"
      ),
      .combine = c,
      .errorhandling = "pass"
    ) %dopar% {
      tryCatch(
        list(process_drug_tglasso(drug_name, data, fold_assignments)),
        error = function(e) {
          list(list(
            drug = drug_name,
            status = "error",
            n_samples = NA_integer_,
            n_tissues = NA_integer_,
            mse_tglasso = NA_real_,
            lambda_by_fold_tissue = NULL,
            selected_features = NULL,
            bootstrap = NULL,
            error = e$message
          ))
        }
      )
    }

    for (res in batch_results) {
      all_results[[res$drug]] <- res
    }
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    cat(sprintf("Progress: %d/%d\n\n", length(completed_drugs), length(all_drugs)))

    append_progress(
      phase = "Phase 2 - TG-LASSO",
      status = "in_progress",
      details = c(
        sprintf("Batch %d/%d", batch_idx, length(batches)),
        sprintf("Progress: %d/%d", length(completed_drugs), length(all_drugs))
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
    n_samples = res$n_samples,
    n_tissues = res$n_tissues,
    mse_tglasso = res$mse_tglasso,
    error = res$error,
    stringsAsFactors = FALSE
  )
}))
summary_df <- summary_df[match(all_drugs, summary_df$drug), , drop = FALSE]

# Merge EN baseline from dsen_all if available.
if (file.exists(DSEN_ALL_SUMMARY)) {
  dsen_all <- read.csv(DSEN_ALL_SUMMARY, stringsAsFactors = FALSE)
  keep <- dsen_all[, intersect(c("drug", "mse_en", "mse_dsen"), names(dsen_all)), drop = FALSE]
  summary_df <- summary_df %>% left_join(keep, by = "drug")
}

payload <- list(
  summary_df = summary_df,
  all_results = all_results,
  config = list(
    seed = GLOBAL_SEED,
    cv_nfolds = CV_NFOLDS,
    min_samples_per_tissue = MIN_SAMPLES_PER_TISSUE,
    min_tissues = MIN_TISSUES,
    feature_subset = FEATURE_SUBSET,
    n_features = length(available_features),
    skip_bootstrap = SKIP_BOOTSTRAP,
    bootstrap_nreps = BOOTSTRAP_NREPS,
    correlation_threshold = CORRELATION_THRESHOLD,
    variance_threshold = VARIANCE_THRESHOLD,
    timestamp = Sys.time()
  )
)

saveRDS(payload, results_rds)
write.csv(summary_df, summary_csv, row.names = FALSE)

bs_out <- lapply(all_results, function(res) {
  list(
    drug = res$drug,
    status = res$status,
    bootstrap = res$bootstrap
  )
})
saveRDS(bs_out, bootstrap_rds)

if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
success_df <- summary_df %>% filter(status == "success")

cat(strrep("=", 80), "\n")
cat("TG-LASSO run complete\n")
cat("Successful drugs:", nrow(success_df), "/", nrow(summary_df), "\n")
cat("Mean MSE:", sprintf("%.6f", mean(success_df$mse_tglasso, na.rm = TRUE)), "\n")
cat("Elapsed (min):", sprintf("%.1f", elapsed), "\n")
cat("Saved summary:", summary_csv, "\n")
cat(strrep("=", 80), "\n")

append_progress(
  phase = "Phase 2 - TG-LASSO",
  status = "completed",
  runtime = sprintf("%.1f minutes", elapsed),
  details = c(
    sprintf("Successful drugs: %d/%d", nrow(success_df), nrow(summary_df)),
    sprintf("Feature subset: %s", FEATURE_SUBSET),
    sprintf("Bootstrap skipped: %s", ifelse(SKIP_BOOTSTRAP, "TRUE", "FALSE")),
    sprintf("Mean TG-LASSO MSE: %.6f", mean(success_df$mse_tglasso, na.rm = TRUE))
  )
)
