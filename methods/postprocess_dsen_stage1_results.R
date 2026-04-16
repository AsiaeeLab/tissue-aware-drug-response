#!/usr/bin/env Rscript
# Purpose: Convert run_dsen_analysis.R output into threshold-specific Stage 1 artifacts.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here("methods", "00_paths.R"))

INPUT_RDA <- Sys.getenv("INPUT_RDA", file.path(paths$scratch, "dsen_results.Rda"))
THRESHOLD_LABEL <- Sys.getenv("THRESHOLD_LABEL", "threshold_50")
OUT_DSEn_DIR <- Sys.getenv(
  "OUT_DSEN_DIR",
  here("results", "dsen_results", THRESHOLD_LABEL)
)
OUT_TABLE_DIR <- Sys.getenv(
  "OUT_TABLE_DIR",
  here("results", "tables", THRESHOLD_LABEL)
)

dir.create(OUT_DSEn_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT_RDA)) {
  stop("Missing Stage 1 input Rda: ", INPUT_RDA)
}

load(INPUT_RDA)
if (!exists("dsen_results")) {
  stop("Expected object 'dsen_results' in ", INPUT_RDA)
}

stage1_df <- dsen_results %>%
  mutate(
    mse_en = mse_standard,
    pct_improvement = mse_reduction_pct,
    alpha_en = alpha_standard,
    lamrat_dsen = NA_real_,
    lambda_en = lambda_standard
  ) %>%
  select(
    drug,
    status,
    n_samples,
    n_tissues,
    mse_en,
    mse_dsen,
    pct_improvement,
    alpha_en,
    alpha_dsen,
    lamrat_dsen,
    lambda_en,
    lambda_dsen
  )

success_df <- stage1_df %>% filter(status == "success")
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

out_rda <- file.path(OUT_DSEn_DIR, "dsen_production_results.Rda")
out_summary_csv <- file.path(OUT_DSEn_DIR, "dsen_summary.csv")
out_full_csv <- file.path(OUT_TABLE_DIR, "dsen_full_results.csv")
out_analysis_csv <- file.path(OUT_TABLE_DIR, "dsen_analysis_summary.csv")

save(dsen_results, file = out_rda)
write_csv(stage1_df, out_summary_csv)
write_csv(success_df, out_full_csv)
write_csv(summary_stats, out_analysis_csv)

cat("Wrote Stage 1 threshold artifacts:\n")
cat("  ", out_rda, "\n", sep = "")
cat("  ", out_summary_csv, "\n", sep = "")
cat("  ", out_full_csv, "\n", sep = "")
cat("  ", out_analysis_csv, "\n", sep = "")
