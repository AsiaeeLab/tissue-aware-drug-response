#!/usr/bin/env Rscript
# Purpose: Compare threshold_25_v2 (sample_size_mild) vs threshold_25_uniform Stage-1
# and select winner mode for final bootstrap.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
})

mild_label <- Sys.getenv("MILD_LABEL", "threshold_25_v2")
uniform_label <- Sys.getenv("UNIFORM_LABEL", "threshold_25_uniform")

mild_csv <- Sys.getenv("MILD_SUMMARY_CSV", here("results", "dsen_results", mild_label, "dsen_summary.csv"))
uniform_csv <- Sys.getenv("UNIFORM_SUMMARY_CSV", here("results", "dsen_results", uniform_label, "dsen_summary.csv"))

out_dir <- Sys.getenv("OUT_DIR", here("results", "tables"))
out_csv <- Sys.getenv("OUT_COMPARISON_CSV", file.path(out_dir, "threshold_25_penalty_mode_comparison.csv"))
out_winner <- Sys.getenv("OUT_WINNER_FILE", file.path(out_dir, "threshold_25_penalty_mode_winner.txt"))

for (f in c(mild_csv, uniform_csv)) {
  if (!file.exists(f)) stop("Missing summary CSV: ", f)
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_winner), recursive = TRUE, showWarnings = FALSE)

read_success <- function(csv_file) {
  df <- read_csv(csv_file, show_col_types = FALSE)
  if (!"status" %in% names(df)) df$status <- "success"
  if (!"pct_improvement" %in% names(df) && all(c("mse_en", "mse_dsen") %in% names(df))) {
    df$pct_improvement <- (df$mse_en - df$mse_dsen) / df$mse_en * 100
  }
  df %>%
    filter(status == "success")
}

summarize_mode <- function(df, mode_name) {
  tibble(
    mode = mode_name,
    n_drugs = nrow(df),
    dsen_win_rate = 100 * mean(df$pct_improvement > 0, na.rm = TRUE),
    dsen_mean_improvement = mean(df$pct_improvement, na.rm = TRUE),
    dsen_median_improvement = median(df$pct_improvement, na.rm = TRUE),
    dsen_max_improvement = max(df$pct_improvement, na.rm = TRUE),
    dsen_max_improvement_drug = df$drug[which.max(df$pct_improvement)]
  )
}

mild_df <- read_success(mild_csv)
uniform_df <- read_success(uniform_csv)

if (nrow(mild_df) == 0 || nrow(uniform_df) == 0) {
  stop("One of the summary files has zero successful drugs; cannot compare modes.")
}

summary_tbl <- bind_rows(
  summarize_mode(mild_df, "sample_size_mild"),
  summarize_mode(uniform_df, "uniform")
)

m <- summary_tbl %>% filter(mode == "sample_size_mild")
u <- summary_tbl %>% filter(mode == "uniform")

winner_mode <- if (
  u$dsen_win_rate > m$dsen_win_rate + 1e-9 ||
    (abs(u$dsen_win_rate - m$dsen_win_rate) < 1e-9 && u$dsen_mean_improvement > m$dsen_mean_improvement)
) {
  "uniform"
} else {
  "sample_size_mild"
}

summary_with_delta <- summary_tbl %>%
  mutate(
    win_rate_delta_vs_mild = dsen_win_rate - m$dsen_win_rate,
    mean_improvement_delta_vs_mild = dsen_mean_improvement - m$dsen_mean_improvement
  )

write_csv(summary_with_delta, out_csv)
writeLines(winner_mode, out_winner)

cat("Wrote:\n")
cat("  ", out_csv, "\n", sep = "")
cat("  ", out_winner, "\n", sep = "")
cat("Winner mode: ", winner_mode, "\n", sep = "")
