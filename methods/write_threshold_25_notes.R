#!/usr/bin/env Rscript
# Purpose: Write threshold-25 comparison note and append Issue 9 to decision log.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
})

comparison_csv <- Sys.getenv(
  "COMPARISON_CSV",
  here("results", "tables", "threshold_comparison_15_25_50.csv")
)
note_file <- Sys.getenv(
  "NOTE_FILE",
  here("doc", "notes", "2026-03-18_threshold_25_comparison.md")
)
decision_log_file <- Sys.getenv(
  "DECISION_LOG_FILE",
  here("doc", "notes", "2026-03-15_tissue_recovery_decisions.md")
)

if (!file.exists(comparison_csv)) {
  stop("Missing comparison CSV: ", comparison_csv)
}

cmp <- read_csv(comparison_csv, show_col_types = FALSE)
if (!all(c("metric", "threshold_15", "threshold_25", "threshold_50") %in% names(cmp))) {
  stop("Comparison CSV missing expected columns.")
}

get_val <- function(metric, col) {
  row <- cmp %>% filter(.data$metric == .env$metric)
  if (nrow(row) == 0) return(NA_character_)
  as.character(row[[col]][1])
}

to_num <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  token <- str_extract(x, "-?[0-9]+\\.?[0-9]*")
  suppressWarnings(as.numeric(token))
}

num_val <- function(metric, col) to_num(get_val(metric, col))

win15 <- num_val("dsen_win_rate", "threshold_15")
win25 <- num_val("dsen_win_rate", "threshold_25")
win50 <- num_val("dsen_win_rate", "threshold_50")

mean15 <- num_val("dsen_mean_improvement", "threshold_15")
mean25 <- num_val("dsen_mean_improvement", "threshold_25")
mean50 <- num_val("dsen_mean_improvement", "threshold_50")

comb50_15 <- num_val("dsen_combined_recovery_50pct", "threshold_15")
comb50_25 <- num_val("dsen_combined_recovery_50pct", "threshold_25")
comb50_50 <- num_val("dsen_combined_recovery_50pct", "threshold_50")

shared_genes_15 <- num_val("mean_dsen_shared_genes", "threshold_15")
shared_genes_25 <- num_val("mean_dsen_shared_genes", "threshold_25")
shared_genes_50 <- num_val("mean_dsen_shared_genes", "threshold_50")

lamrat_lt1_25 <- get_val("lamrat_lt1_pct", "threshold_25")
lamrat_dist_25 <- get_val("lamrat_distribution", "threshold_25")

q1_25 <- get_val("difficulty_q1_win_rate", "threshold_25")
q2_25 <- get_val("difficulty_q2_win_rate", "threshold_25")
q3_25 <- get_val("difficulty_q3_win_rate", "threshold_25")
q4_25 <- get_val("difficulty_q4_win_rate", "threshold_25")

recommendation <- if (
  is.finite(win25) && is.finite(win15) && is.finite(comb50_25) && is.finite(comb50_15) &&
  win25 >= win15 && comb50_25 >= comb50_15
) {
  "Use threshold 25 as primary, with thresholds 15 and 50 as sensitivities."
} else {
  "Keep threshold 15 as primary; present threshold 25 as a targeted sensitivity analysis with corrected lamrat tuning."
}

note_lines <- c(
  "# Threshold Comparison: 15 vs 25 (and 50 reference)",
  "",
  sprintf("Generated on %s from `%s`.", Sys.Date(), comparison_csv),
  "",
  "## Headline Results",
  sprintf("- DSEN win rate: 15=%s, 25=%s, 50=%s",
          get_val("dsen_win_rate", "threshold_15"),
          get_val("dsen_win_rate", "threshold_25"),
          get_val("dsen_win_rate", "threshold_50")),
  sprintf("- Mean DSEN improvement: 15=%s, 25=%s, 50=%s",
          get_val("dsen_mean_improvement", "threshold_15"),
          get_val("dsen_mean_improvement", "threshold_25"),
          get_val("dsen_mean_improvement", "threshold_50")),
  sprintf("- DSEN combined target recovery @50%%: 15=%s, 25=%s, 50=%s",
          get_val("dsen_combined_recovery_50pct", "threshold_15"),
          get_val("dsen_combined_recovery_50pct", "threshold_25"),
          get_val("dsen_combined_recovery_50pct", "threshold_50")),
  sprintf("- Mean DSEN shared genes: 15=%s, 25=%s, 50=%s",
          get_val("mean_dsen_shared_genes", "threshold_15"),
          get_val("mean_dsen_shared_genes", "threshold_25"),
          get_val("mean_dsen_shared_genes", "threshold_50")),
  "",
  "## Lamrat Behavior at Threshold 25",
  sprintf("- lamrat < 1 selected for %s of successful drugs.", lamrat_lt1_25),
  sprintf("- Full lamrat distribution: %s", lamrat_dist_25),
  "",
  "## Difficulty Quartiles (Threshold 25)",
  sprintf("- Q1 win rate: %s", q1_25),
  sprintf("- Q2 win rate: %s", q2_25),
  sprintf("- Q3 win rate: %s", q3_25),
  sprintf("- Q4 win rate: %s", q4_25),
  "",
  "## Interpretation",
  sprintf("- Win-rate deltas: 25-15 = %.2f pp, 25-50 = %.2f pp.", win25 - win15, win25 - win50),
  sprintf("- Mean-improvement deltas: 25-15 = %.2f pp, 25-50 = %.2f pp.", mean25 - mean15, mean25 - mean50),
  sprintf("- Combined-recovery deltas @50%%: 25-15 = %.2f pp, 25-50 = %.2f pp.", comb50_25 - comb50_15, comb50_25 - comb50_50),
  sprintf("- Shared-gene count deltas: 25-15 = %.2f, 25-50 = %.2f.", shared_genes_25 - shared_genes_15, shared_genes_25 - shared_genes_50),
  "",
  "## Recommendation",
  paste0("- ", recommendation)
)

dir.create(dirname(note_file), recursive = TRUE, showWarnings = FALSE)
writeLines(note_lines, note_file)

issue_header <- "## Issue 9: Threshold 25 with corrected lamrat grid"
issue_exists <- FALSE
if (file.exists(decision_log_file)) {
  issue_exists <- any(str_detect(readLines(decision_log_file, warn = FALSE), fixed(issue_header)))
}

if (!issue_exists) {
  issue_block <- c(
    "",
    "---",
    "",
    issue_header,
    "",
    sprintf("**Date:** %s", Sys.Date()),
    "",
    "**What changed vs threshold-50 attempt:**",
    "- Used `run_dsen_production.R` (alpha + lamrat tuning) instead of the non-lamrat analysis runner.",
    "- Extended default lamrat grid to include values < 1 (`0.25, 0.5, 1, 2, 4`).",
    "- Applied threshold 25 tissue filtering (min 25 samples per tissue).",
    "",
    "**Results snapshot (15 / 25 / 50):**",
    sprintf("- DSEN win rate: %s / %s / %s",
            get_val("dsen_win_rate", "threshold_15"),
            get_val("dsen_win_rate", "threshold_25"),
            get_val("dsen_win_rate", "threshold_50")),
    sprintf("- Mean DSEN improvement: %s / %s / %s",
            get_val("dsen_mean_improvement", "threshold_15"),
            get_val("dsen_mean_improvement", "threshold_25"),
            get_val("dsen_mean_improvement", "threshold_50")),
    sprintf("- DSEN combined recovery @50%%: %s / %s / %s",
            get_val("dsen_combined_recovery_50pct", "threshold_15"),
            get_val("dsen_combined_recovery_50pct", "threshold_25"),
            get_val("dsen_combined_recovery_50pct", "threshold_50")),
    "",
    sprintf("**Lamrat distribution at threshold 25:** %s", lamrat_dist_25),
    sprintf("**Lamrat < 1 frequency:** %s", lamrat_lt1_25),
    "",
    paste0("**Recommendation:** ", recommendation)
  )

  dir.create(dirname(decision_log_file), recursive = TRUE, showWarnings = FALSE)
  cat(paste(issue_block, collapse = "\n"), file = decision_log_file, append = TRUE)
}

cat("Wrote note:", note_file, "\n")
if (issue_exists) {
  cat("Issue 9 already exists in decision log; left unchanged:", decision_log_file, "\n")
} else {
  cat("Appended Issue 9 to decision log:", decision_log_file, "\n")
}
