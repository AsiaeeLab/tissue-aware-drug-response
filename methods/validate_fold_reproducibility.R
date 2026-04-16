#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("methods", "02_cv_utils.R"))
source(here("methods", "comparison_common.R"))

GLOBAL_SEED <- 42L
CV_NFOLDS <- 5L
MIN_SAMPLES_PER_TISSUE <- 25L
VARIANCE_THRESHOLD <- 0.01
CORRELATION_THRESHOLD <- 0.1
ALPHA_GRID <- c(0.05, 0.1, 0.15, 0.2, seq(0.3, 1.0, 0.1))

out_dir <- ensure_comparison_dir()
fold_file <- file.path(out_dir, "fold_assignments.rds")
out_file <- file.path(out_dir, "fold_reproducibility_check.csv")

if (!file.exists(fold_file)) {
  stop("Missing fold assignment file: ", fold_file)
}

prod_summary_file <- here("results", "dsen_results", "threshold_25_v2", "dsen_summary.csv")
if (!file.exists(prod_summary_file)) {
  stop("Missing production summary file: ", prod_summary_file)
}

fold_assignments <- readRDS(fold_file)
prod <- read.csv(prod_summary_file, stringsAsFactors = FALSE)
prod <- prod %>% filter(status == "success")

valid_drugs <- intersect(prod$drug, names(fold_assignments))
valid_drugs <- valid_drugs[vapply(valid_drugs, function(d) {
  is.list(fold_assignments[[d]]) && identical(fold_assignments[[d]]$status, "ok")
}, logical(1))]

if (length(valid_drugs) < 5) {
  stop("Need at least 5 valid drugs for reproducibility check.")
}

set.seed(GLOBAL_SEED)
check_drugs <- sample(valid_drugs, 5)

message("Loading data...")
data <- load_main_data("Sanger")

rows <- lapply(check_drugs, function(drug) {
  message(sprintf("Checking %s", drug))

  fd <- fold_assignments[[drug]]
  dd <- prepare_drug_data_filtered(data, drug, min_per_tissue = MIN_SAMPLES_PER_TISSUE)

  if (is.null(dd)) {
    return(data.frame(
      drug = drug,
      mse_en_recomputed = NA_real_,
      mse_en_production = NA_real_,
      abs_diff = NA_real_,
      match = FALSE,
      note = "prepare_drug_data_filtered returned NULL",
      stringsAsFactors = FALSE
    ))
  }

  if (!identical(rownames(dd$x), fd$cell_line_ids)) {
    return(data.frame(
      drug = drug,
      mse_en_recomputed = NA_real_,
      mse_en_production = NA_real_,
      abs_diff = NA_real_,
      match = FALSE,
      note = "cell_line_ids mismatch",
      stringsAsFactors = FALSE
    ))
  }

  en <- cv_tune_alpha(
    dd$x,
    dd$y,
    alphas = ALPHA_GRID,
    nfolds = CV_NFOLDS,
    seed = GLOBAL_SEED,
    foldid = as.integer(fd$foldid)
  )

  ref <- prod$mse_en[prod$drug == drug][1]
  est <- en$best_cvm
  diff <- abs(est - ref)

  data.frame(
    drug = drug,
    mse_en_recomputed = est,
    mse_en_production = ref,
    abs_diff = diff,
    match = diff < 1e-8,
    note = ifelse(diff < 1e-8, "ok", "mismatch"),
    stringsAsFactors = FALSE
  )
})

check_df <- bind_rows(rows)
write.csv(check_df, out_file, row.names = FALSE)

message("Fold reproducibility check complete")
message(sprintf("  Matches: %d/%d", sum(check_df$match), nrow(check_df)))
message(sprintf("  Output: %s", out_file))

append_progress(
  phase = "Phase 0.3 - Fold reproducibility",
  status = if (all(check_df$match)) "completed" else "completed_with_warnings",
  details = c(
    sprintf("Drugs checked: %d", nrow(check_df)),
    sprintf("Matches (abs diff < 1e-8): %d", sum(check_df$match)),
    sprintf("Mismatches: %d", sum(!check_df$match))
  )
)
