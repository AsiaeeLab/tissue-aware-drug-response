#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(caret)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "03_dsen_utils.R"))
source(here("methods", "comparison_common.R"))

# =============================================================================
# CONFIG
# =============================================================================
GLOBAL_SEED <- 42L
CV_NFOLDS <- 5L
MIN_SAMPLES_PER_TISSUE <- 25L
MIN_SAMPLES_TOTAL <- 100L
MIN_TISSUES_LTO <- 5L
MIN_TISSUES_MODEL <- 3L

CORRELATION_THRESHOLD <- 0.1
VARIANCE_THRESHOLD <- 0.01
ALPHA_GRID <- c(0.05, 0.1, 0.15, 0.2, seq(0.3, 1.0, 0.1))
LAMRAT_GRID <- c(0.125, 0.25, 0.5, 1, 2, 4)
SCREENING <- "pooled"
TISSUE_PENALTY_MODE <- "sample_size_mild"

NCORES <- as.integer(Sys.getenv("NCORES", "10"))
if (is.na(NCORES) || NCORES < 1) NCORES <- 10L
NCORES <- min(NCORES, 10L)

CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "2"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1) CHECKPOINT_EVERY <- 2L

max_drugs_env <- Sys.getenv("MAX_DRUGS", "")
MAX_DRUGS <- if (nzchar(max_drugs_env)) as.integer(max_drugs_env) else Inf

TUNE_PER_SPLIT <- Sys.getenv("TUNE_PER_SPLIT", "TRUE") == "TRUE"

DSEN_ALL_SUMMARY <- Sys.getenv(
  "DSEN_ALL_SUMMARY",
  here("results", "comparison", "ablation", "dsen_all", "dsen_summary.csv")
)

OUT_DIR <- Sys.getenv("OUT_DIR", here("results", "comparison", "leave_tissue_out"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results_csv <- file.path(OUT_DIR, "lto_results.csv")
summary_csv <- file.path(OUT_DIR, "lto_summary.csv")
results_rds <- file.path(OUT_DIR, "lto_results.rds")
checkpoint_file <- file.path(OUT_DIR, "lto_checkpoint.Rda")

cat(strrep("=", 80), "\n")
cat("LEAVE-TISSUE-OUT VALIDATION\n")
cat(strrep("=", 80), "\n")
cat("Tune per split:", TUNE_PER_SPLIT, "\n")
cat("Cores:", NCORES, "\n")
cat("Output dir:", OUT_DIR, "\n")
cat(strrep("=", 80), "\n\n")

if (!TUNE_PER_SPLIT && !file.exists(DSEN_ALL_SUMMARY)) {
  stop("TUNE_PER_SPLIT=FALSE requires DSEN_ALL_SUMMARY: ", DSEN_ALL_SUMMARY)
}

cached_params <- NULL
if (file.exists(DSEN_ALL_SUMMARY)) {
  cached_params <- read.csv(DSEN_ALL_SUMMARY, stringsAsFactors = FALSE)
}

get_cached_params <- function(drug_name, cache_df) {
  if (is.null(cache_df)) return(NULL)
  row <- cache_df[cache_df$drug == drug_name & cache_df$status == "success", , drop = FALSE]
  if (nrow(row) == 0) return(NULL)

  out <- list(
    alpha_en = as.numeric(row$alpha_en[1]),
    lambda_en = as.numeric(row$lambda_en[1]),
    alpha_dsen = as.numeric(row$alpha_dsen[1]),
    lamrat_dsen = as.numeric(row$lamrat_dsen[1]),
    lambda_dsen = as.numeric(row$lambda_dsen[1])
  )

  if (any(!is.finite(unlist(out)))) return(NULL)
  out
}

fit_en_on_split <- function(x_train, y_train, x_test, alpha_en, lambda_en) {
  pre <- preprocess_features(
    x_train, y_train,
    var_thresh = VARIANCE_THRESHOLD,
    cor_thresh = CORRELATION_THRESHOLD,
    scale_features = TRUE
  )
  if (is.null(pre) || ncol(pre$x) < 2) {
    return(rep(mean(y_train, na.rm = TRUE), nrow(x_test)))
  }

  x_train_proc <- pre$x
  x_test_proc <- apply_preprocessing(x_test, pre)
  common <- intersect(colnames(x_train_proc), colnames(x_test_proc))
  if (length(common) < 2) {
    return(rep(mean(y_train, na.rm = TRUE), nrow(x_test)))
  }

  x_train_proc <- x_train_proc[, common, drop = FALSE]
  x_test_proc <- x_test_proc[, common, drop = FALSE]

  fit <- glmnet::glmnet(x_train_proc, y_train, alpha = alpha_en, lambda = lambda_en)
  as.numeric(predict(fit, newx = x_test_proc, s = lambda_en))
}

fit_dsen_on_split <- function(x_train, y_train, tissue_train, x_test,
                              alpha_dsen, lambda_dsen, lamrat_dsen) {
  train_idx <- seq_len(nrow(x_train))

  prep <- prepare_tissue_data(
    x_train, y_train, tissue_train,
    train_idx = train_idx,
    var_thresh = VARIANCE_THRESHOLD,
    cor_thresh = CORRELATION_THRESHOLD,
    screening = SCREENING
  )

  data_list <- prep$data_list
  if (length(data_list) < MIN_TISSUES_MODEL) {
    pred <- rep(mean(y_train, na.rm = TRUE), nrow(x_test))
    return(list(shared = pred, full = pred))
  }

  p <- ncol(data_list[[1]]$x)
  Z <- build_dsen_matrix(data_list)
  yz <- unlist(lapply(data_list, function(item) item$y))
  pf <- compute_dsen_penalties(data_list, lamrat_dsen, TISSUE_PENALTY_MODE)

  fit <- glmnet::glmnet(
    Z, yz,
    alpha = alpha_dsen,
    lambda = lambda_dsen,
    penalty.factor = pf,
    intercept = TRUE,
    standardize = FALSE
  )

  coef_vec <- as.vector(coef(fit, s = lambda_dsen))
  intercept <- coef_vec[1]
  beta <- coef_vec[-1]
  beta_shared <- beta[1:p]

  x_test_proc <- as_numeric_matrix(x_test[, prep$selected_features, drop = FALSE])
  x_test_proc[is.na(x_test_proc)] <- 0

  if (!is.null(prep$scale_center)) {
    is_mut <- grepl("_mut(\\.|$)|^mut_", colnames(x_test_proc), ignore.case = TRUE)
    for (feat in names(prep$scale_center)) {
      if (feat %in% colnames(x_test_proc) && !is_mut[colnames(x_test_proc) == feat]) {
        x_test_proc[, feat] <- (x_test_proc[, feat] - prep$scale_center[feat]) / prep$scale_scale[feat]
      }
    }
  }

  shared_pred <- as.numeric(x_test_proc %*% beta_shared + intercept)

  # Held-out tissue is not in training data, so tissue-specific block is zero.
  full_pred <- shared_pred

  list(shared = shared_pred, full = full_pred)
}

process_drug_lto <- function(drug_name, data, cache_df) {
  dd <- prepare_drug_data_filtered(data, drug_name, min_per_tissue = MIN_SAMPLES_PER_TISSUE)
  if (is.null(dd)) {
    return(data.frame(
      drug = drug_name,
      held_out_tissue = NA_character_,
      n_test = NA_integer_,
      mse_en = NA_real_,
      mse_dsen_shared = NA_real_,
      mse_dsen_full = NA_real_,
      status = "skipped_insufficient_tissues",
      stringsAsFactors = FALSE
    ))
  }

  tissue <- as.character(dd$tissue)
  tissue_counts <- table(tissue)
  valid_holdout_tissues <- names(tissue_counts)[tissue_counts >= MIN_SAMPLES_PER_TISSUE]

  if (length(valid_holdout_tissues) < MIN_TISSUES_LTO || length(y <- dd$y) < MIN_SAMPLES_TOTAL) {
    return(data.frame(
      drug = drug_name,
      held_out_tissue = NA_character_,
      n_test = NA_integer_,
      mse_en = NA_real_,
      mse_dsen_shared = NA_real_,
      mse_dsen_full = NA_real_,
      status = "skipped_lto_criteria",
      stringsAsFactors = FALSE
    ))
  }

  x <- dd$x

  if (!TUNE_PER_SPLIT) {
    cached <- get_cached_params(drug_name, cache_df)
    if (is.null(cached)) {
      return(data.frame(
        drug = drug_name,
        held_out_tissue = NA_character_,
        n_test = NA_integer_,
        mse_en = NA_real_,
        mse_dsen_shared = NA_real_,
        mse_dsen_full = NA_real_,
        status = "skipped_no_cached_params",
        stringsAsFactors = FALSE
      ))
    }
  }

  out_rows <- list()

  for (tis in valid_holdout_tissues) {
    test_idx <- which(tissue == tis)
    train_idx <- which(tissue != tis)

    x_train <- x[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    t_train <- tissue[train_idx]

    x_test <- x[test_idx, , drop = FALSE]
    y_test <- y[test_idx]

    if (length(unique(t_train)) < MIN_TISSUES_MODEL) {
      out_rows[[length(out_rows) + 1]] <- data.frame(
        drug = drug_name,
        held_out_tissue = tis,
        n_test = length(y_test),
        mse_en = NA_real_,
        mse_dsen_shared = NA_real_,
        mse_dsen_full = NA_real_,
        status = "skipped_training_tissues",
        stringsAsFactors = FALSE
      )
      next
    }

    params <- NULL
    if (TUNE_PER_SPLIT) {
      foldid_train <- caret::createFolds(t_train, k = CV_NFOLDS, list = FALSE)

      en_tuned <- tryCatch(
        cv_tune_alpha(
          x_train,
          y_train,
          alphas = ALPHA_GRID,
          nfolds = CV_NFOLDS,
          seed = GLOBAL_SEED,
          foldid = foldid_train
        ),
        error = function(e) NULL
      )

      dsen_tuned <- tryCatch(
        dsen_tune_alpha_lamrat(
          x_train,
          y_train,
          t_train,
          alphas = ALPHA_GRID,
          lamrats = LAMRAT_GRID,
          nfolds = CV_NFOLDS,
          seed = GLOBAL_SEED,
          foldid = foldid_train,
          tissue_penalty_mode = TISSUE_PENALTY_MODE,
          screening = SCREENING
        ),
        error = function(e) NULL
      )

      if (is.null(en_tuned) || is.null(dsen_tuned)) {
        out_rows[[length(out_rows) + 1]] <- data.frame(
          drug = drug_name,
          held_out_tissue = tis,
          n_test = length(y_test),
          mse_en = NA_real_,
          mse_dsen_shared = NA_real_,
          mse_dsen_full = NA_real_,
          status = "error_tuning",
          stringsAsFactors = FALSE
        )
        next
      }

      params <- list(
        alpha_en = en_tuned$alpha.min,
        lambda_en = en_tuned$lambda.min,
        alpha_dsen = dsen_tuned$alpha.min,
        lamrat_dsen = dsen_tuned$lamrat.min,
        lambda_dsen = dsen_tuned$lambda.min
      )
    } else {
      params <- get_cached_params(drug_name, cache_df)
    }

    en_pred <- tryCatch(
      fit_en_on_split(x_train, y_train, x_test, params$alpha_en, params$lambda_en),
      error = function(e) rep(NA_real_, nrow(x_test))
    )

    dsen_pred <- tryCatch(
      fit_dsen_on_split(
        x_train, y_train, t_train, x_test,
        params$alpha_dsen,
        params$lambda_dsen,
        params$lamrat_dsen
      ),
      error = function(e) list(shared = rep(NA_real_, nrow(x_test)), full = rep(NA_real_, nrow(x_test)))
    )

    out_rows[[length(out_rows) + 1]] <- data.frame(
      drug = drug_name,
      held_out_tissue = tis,
      n_test = length(y_test),
      mse_en = mean((y_test - en_pred)^2, na.rm = TRUE),
      mse_dsen_shared = mean((y_test - dsen_pred$shared)^2, na.rm = TRUE),
      mse_dsen_full = mean((y_test - dsen_pred$full)^2, na.rm = TRUE),
      status = "success",
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(out_rows)
}

cat("Loading data...\n")
data <- load_main_data("Sanger")
all_drugs <- get_drug_list(data)
if (!is.infinite(MAX_DRUGS)) {
  all_drugs <- all_drugs[1:min(MAX_DRUGS, length(all_drugs))]
}

if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  remaining_drugs <- setdiff(all_drugs, completed_drugs)
  cat("Resuming checkpoint. Completed drugs:", length(completed_drugs), "\n")
} else {
  completed_drugs <- character(0)
  all_results <- list()
  remaining_drugs <- all_drugs
  cat("Starting fresh\n")
}

cat("Remaining drugs:", length(remaining_drugs), "\n\n")

start_time <- Sys.time()

if (length(remaining_drugs) > 0) {
  batches <- split(
    remaining_drugs,
    ceiling(seq_along(remaining_drugs) / CHECKPOINT_EVERY)
  )

  for (batch_idx in seq_along(batches)) {
    batch <- batches[[batch_idx]]
    cat(sprintf("[%s] Batch %d/%d (%d drugs)\n",
                format(Sys.time(), "%H:%M:%S"), batch_idx, length(batches), length(batch)))

    batch_results <- parallel::mclapply(
      batch,
      function(drug_name) {
        tryCatch(
          process_drug_lto(drug_name, data, cached_params),
          error = function(e) data.frame(
            drug = drug_name,
            held_out_tissue = NA_character_,
            n_test = NA_integer_,
            mse_en = NA_real_,
            mse_dsen_shared = NA_real_,
            mse_dsen_full = NA_real_,
            status = paste0("error: ", e$message),
            stringsAsFactors = FALSE
          )
        )
      },
      mc.cores = NCORES
    )

    names(batch_results) <- batch
    all_results <- c(all_results, batch_results)
    completed_drugs <- c(completed_drugs, batch)

    save(all_results, completed_drugs, file = checkpoint_file)

    cat(sprintf("Progress: %d/%d\n\n", length(completed_drugs), length(all_drugs)))

    append_progress(
      phase = "Phase 3 - Leave tissue out",
      status = "in_progress",
      details = c(
        sprintf("Batch %d/%d", batch_idx, length(batches)),
        sprintf("Progress: %d/%d", length(completed_drugs), length(all_drugs))
      )
    )
  }
}

lto_df <- bind_rows(all_results)
write.csv(lto_df, results_csv, row.names = FALSE)

lto_summary <- lto_df %>%
  filter(status == "success") %>%
  group_by(drug) %>%
  summarize(
    n_tissues_tested = n(),
    mean_mse_en = mean(mse_en, na.rm = TRUE),
    mean_mse_dsen_shared = mean(mse_dsen_shared, na.rm = TRUE),
    mean_mse_dsen_full = mean(mse_dsen_full, na.rm = TRUE),
    dsen_shared_win_rate = mean(mse_dsen_shared < mse_en, na.rm = TRUE),
    dsen_full_win_rate = mean(mse_dsen_full < mse_en, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(lto_summary, summary_csv, row.names = FALSE)
saveRDS(list(results = lto_df, summary = lto_summary), results_rds)

if (file.exists(checkpoint_file)) file.remove(checkpoint_file)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat(strrep("=", 80), "\n")
cat("LTO run complete\n")
cat("Drugs with successful LTO rows:", n_distinct(lto_df$drug[lto_df$status == "success"]), "\n")
cat("Total successful tissue-holdout rows:", sum(lto_df$status == "success"), "\n")
cat("Elapsed (min):", sprintf("%.1f", elapsed), "\n")
cat("Saved:", results_csv, "\n")
cat(strrep("=", 80), "\n")

append_progress(
  phase = "Phase 3 - Leave tissue out",
  status = "completed",
  runtime = sprintf("%.1f minutes", elapsed),
  details = c(
    sprintf("Tune per split: %s", ifelse(TUNE_PER_SPLIT, "TRUE", "FALSE")),
    sprintf("Successful tissue rows: %d", sum(lto_df$status == "success")),
    sprintf("Drugs represented: %d", n_distinct(lto_df$drug[lto_df$status == "success"]))
  )
)
