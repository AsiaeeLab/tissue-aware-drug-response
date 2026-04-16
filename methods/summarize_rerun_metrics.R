#!/usr/bin/env Rscript
# Purpose: Summarize key manuscript metrics from rerun outputs.
# Output: results/tables/rerun_metric_summary.csv

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
})

source(here("methods", "00_paths.R"))

cv_file <- file.path(paths$scratch, "cv_wrong_vs_right_results.Rda")
dsen_file <- file.path(paths$scratch, "dsen_results.Rda")
out_file <- file.path(paths$results, "tables", "rerun_metric_summary.csv")

if (!file.exists(cv_file)) {
  stop("Missing CV rerun output: ", cv_file)
}
if (!file.exists(dsen_file)) {
  stop("Missing DSEN rerun output: ", dsen_file)
}

load(cv_file)   # expects cv_results
load(dsen_file) # expects dsen_results

if (!exists("cv_results")) stop("Object 'cv_results' not found in ", cv_file)
if (!exists("dsen_results")) stop("Object 'dsen_results' not found in ", dsen_file)

cv_ok <- cv_results %>% filter(status == "success")
dsen_ok <- dsen_results %>% filter(status == "success")

summary_tbl <- tibble(
  metric = c(
    "cv_n_success",
    "cv_mean_jaccard",
    "cv_median_jaccard",
    "cv_pct_zero_overlap",
    "cv_mean_n_features_wrong",
    "cv_mean_n_features_right",
    "dsen_n_success",
    "dsen_win_rate_pct",
    "dsen_mean_mse_improvement_pct"
  ),
  value = c(
    nrow(cv_ok),
    mean(cv_ok$jaccard, na.rm = TRUE),
    median(cv_ok$jaccard, na.rm = TRUE),
    mean(cv_ok$jaccard == 0, na.rm = TRUE) * 100,
    mean(cv_ok$n_features_wrong, na.rm = TRUE),
    mean(cv_ok$n_features_right, na.rm = TRUE),
    nrow(dsen_ok),
    mean(dsen_ok$mse_reduction_pct > 0, na.rm = TRUE) * 100,
    mean(dsen_ok$mse_reduction_pct, na.rm = TRUE)
  )
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, out_file)
cat("Wrote:", out_file, "\n")
