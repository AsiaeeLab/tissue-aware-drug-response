#!/usr/bin/env Rscript
# Purpose: Build threshold-specific Stage 1 summary tables from run_dsen_production output.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
})

threshold_label <- Sys.getenv("THRESHOLD_LABEL", "threshold_25")
input_rda <- Sys.getenv(
  "INPUT_RDA",
  here("results", "dsen_results", threshold_label, "dsen_production_results.Rda")
)
out_table_dir <- Sys.getenv(
  "OUT_TABLE_DIR",
  here("results", "tables", threshold_label)
)

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rda)) {
  stop("Missing input Rda: ", input_rda)
}

env <- new.env(parent = emptyenv())
load(input_rda, envir = env)

if (!exists("production_results", envir = env, inherits = FALSE)) {
  stop("Expected object 'production_results' in ", input_rda)
}
production_results <- get("production_results", envir = env)

if (!"summary_df" %in% names(production_results)) {
  stop("'production_results' missing summary_df in ", input_rda)
}

summary_df <- production_results$summary_df
if (!"lambda_en" %in% names(summary_df)) summary_df$lambda_en <- NA_real_
if (!"lambda_dsen" %in% names(summary_df)) summary_df$lambda_dsen <- NA_real_

full_df <- summary_df %>%
  select(
    drug, status, n_samples, n_tissues,
    mse_en, mse_dsen, pct_improvement,
    alpha_en, alpha_dsen, lamrat_dsen,
    lambda_en, lambda_dsen
  )

success_df <- full_df %>% filter(status == "success")
if (nrow(success_df) == 0) {
  stop("No successful drugs found in Stage 1 results.")
}

max_idx <- which.max(success_df$pct_improvement)
max_drug <- success_df$drug[max_idx]
max_val <- success_df$pct_improvement[max_idx]

summary_stats <- tribble(
  ~Metric, ~Value,
  "Total drugs analyzed", as.character(nrow(success_df)),
  "Mean samples per drug", sprintf("%.0f", mean(success_df$n_samples, na.rm = TRUE)),
  "Mean tissues per drug", sprintf("%.1f", mean(success_df$n_tissues, na.rm = TRUE)),
  "---", "---",
  "DSEN win rate", sprintf("%.1f%%", 100 * mean(success_df$pct_improvement > 0, na.rm = TRUE)),
  "Mean MSE improvement", sprintf("%.2f%%", mean(success_df$pct_improvement, na.rm = TRUE)),
  "Median MSE improvement", sprintf("%.2f%%", median(success_df$pct_improvement, na.rm = TRUE)),
  "Max improvement", sprintf("%.1f%% (%s)", max_val, max_drug)
)

out_full_csv <- file.path(out_table_dir, "dsen_full_results.csv")
out_summary_csv <- file.path(out_table_dir, "dsen_analysis_summary.csv")

write_csv(full_df, out_full_csv)
write_csv(summary_stats, out_summary_csv)

cat("Wrote:\n")
cat("  ", out_full_csv, "\n", sep = "")
cat("  ", out_summary_csv, "\n", sep = "")
