#!/usr/bin/env Rscript
# Purpose: Write threshold-15 vs threshold-50 narrative notes and append decision log entry.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
})

comparison_csv <- Sys.getenv(
  "COMPARISON_CSV",
  here("results", "tables", "threshold_comparison_15_vs_50.csv")
)
note_file <- Sys.getenv(
  "NOTE_FILE",
  here("doc", "notes", "2026-03-18_threshold_comparison.md")
)
decision_log_file <- Sys.getenv(
  "DECISION_LOG_FILE",
  here("doc", "notes", "2026-03-15_tissue_recovery_decisions.md")
)

if (!file.exists(comparison_csv)) {
  stop("Missing comparison CSV: ", comparison_csv)
}

cmp <- read_csv(comparison_csv, show_col_types = FALSE)
if (!all(c("metric", "threshold_15", "threshold_50") %in% names(cmp))) {
  stop("Comparison CSV missing expected columns.")
}

get_val <- function(metric, col) {
  row <- cmp %>% filter(.data$metric == .env$metric)
  if (nrow(row) == 0) return(NA_character_)
  as.character(row[[col]][1])
}

to_num <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  # Keep only the first numeric token (handles strings like "23.6% (Imatinib)")
  token <- str_extract(x, "-?[0-9]+\\.?[0-9]*")
  suppressWarnings(as.numeric(token))
}

metric_num <- function(metric, col) to_num(get_val(metric, col))
metric_chr <- function(metric, col) get_val(metric, col)

win15 <- metric_num("dsen_win_rate", "threshold_15")
win50 <- metric_num("dsen_win_rate", "threshold_50")
mean15 <- metric_num("dsen_mean_improvement", "threshold_15")
mean50 <- metric_num("dsen_mean_improvement", "threshold_50")
med15 <- metric_num("dsen_median_improvement", "threshold_15")
med50 <- metric_num("dsen_median_improvement", "threshold_50")

tr50_15 <- metric_num("dsen_tissue_recovery_50pct", "threshold_15")
tr50_50 <- metric_num("dsen_tissue_recovery_50pct", "threshold_50")
tc50_15 <- metric_num("dsen_combined_recovery_50pct", "threshold_15")
tc50_50 <- metric_num("dsen_combined_recovery_50pct", "threshold_50")

win_delta <- win50 - win15
mean_delta <- mean50 - mean15
med_delta <- med50 - med15
tr_delta <- tr50_50 - tr50_15
tc_delta <- tc50_50 - tc50_15

primary_recommendation <- if (
  is.finite(win_delta) && is.finite(mean_delta) &&
  win_delta >= 0 && mean_delta >= 0
) {
  "Use threshold 50 as the primary analysis and keep threshold 15 as sensitivity."
} else if (
  is.finite(win_delta) && is.finite(mean_delta) &&
  win_delta <= 0 && mean_delta <= 0
) {
  "Keep threshold 15 as primary and treat threshold 50 as sensitivity."
} else {
  "Mixed signal: keep threshold 15 as primary and report threshold 50 as a prespecified sensitivity analysis."
}

improved <- character(0)
worsened <- character(0)

add_change <- function(metric_label, old_val, new_val) {
  if (!is.finite(old_val) || !is.finite(new_val)) return(invisible(NULL))
  line <- sprintf("%s: %.2f -> %.2f (delta %.2f)", metric_label, old_val, new_val, new_val - old_val)
  if (new_val > old_val) {
    improved <<- c(improved, line)
  } else if (new_val < old_val) {
    worsened <<- c(worsened, line)
  }
}

add_change("DSEN win rate (%)", win15, win50)
add_change("Mean DSEN improvement (%)", mean15, mean50)
add_change("Median DSEN improvement (%)", med15, med50)
add_change("DSEN tissue recovery @50% threshold (%)", tr50_15, tr50_50)
add_change("DSEN combined recovery @50% threshold (%)", tc50_15, tc50_50)

if (length(improved) == 0) improved <- "- None"
if (length(worsened) == 0) worsened <- "- None"

note_lines <- c(
  "# Threshold Sensitivity: 15 vs 50",
  "",
  sprintf("Generated on %s from `%s`.", Sys.Date(), comparison_csv),
  "",
  "## Summary",
  sprintf("- DSEN win rate: %s (15) vs %s (50)", metric_chr("dsen_win_rate", "threshold_15"), metric_chr("dsen_win_rate", "threshold_50")),
  sprintf("- Mean DSEN improvement: %s (15) vs %s (50)", metric_chr("dsen_mean_improvement", "threshold_15"), metric_chr("dsen_mean_improvement", "threshold_50")),
  sprintf("- Median DSEN improvement: %s (15) vs %s (50)", metric_chr("dsen_median_improvement", "threshold_15"), metric_chr("dsen_median_improvement", "threshold_50")),
  sprintf("- DSEN tissue target recovery @50%%: %s (15) vs %s (50)", metric_chr("dsen_tissue_recovery_50pct", "threshold_15"), metric_chr("dsen_tissue_recovery_50pct", "threshold_50")),
  sprintf("- DSEN combined target recovery @50%%: %s (15) vs %s (50)", metric_chr("dsen_combined_recovery_50pct", "threshold_15"), metric_chr("dsen_combined_recovery_50pct", "threshold_50")),
  "",
  "## Metrics That Improved",
  improved,
  "",
  "## Metrics That Worsened",
  worsened,
  "",
  "## Interpretation",
  sprintf("- DSEN advantage (win-rate delta): %.2f percentage points.", win_delta),
  sprintf("- DSEN advantage (mean-improvement delta): %.2f percentage points.", mean_delta),
  sprintf("- Tissue-specific recovery change at 50%% threshold: %.2f percentage points.", tr_delta),
  sprintf("- Combined recovery change at 50%% threshold: %.2f percentage points.", tc_delta),
  "",
  "## Recommendation",
  paste0("- ", primary_recommendation),
  "- Report both thresholds in the manuscript/supplement as robustness evidence."
)

dir.create(dirname(note_file), recursive = TRUE, showWarnings = FALSE)
writeLines(note_lines, note_file)

decision_header <- "## Issue 8: Tissue threshold sensitivity analysis"
decision_exists <- FALSE
if (file.exists(decision_log_file)) {
  decision_exists <- any(str_detect(readLines(decision_log_file, warn = FALSE), fixed(decision_header)))
}

if (!decision_exists) {
  decision_block <- c(
    "",
    "---",
    "",
    decision_header,
    "",
    sprintf("**Date:** %s", Sys.Date()),
    "",
    "**Rationale:** Small tissues can destabilize tissue-specific coefficient estimates in DSEN. Raising the minimum tissue sample threshold from 15 to 50 was evaluated as a sensitivity check.",
    "",
    "**Results snapshot (15 -> 50):**",
    sprintf("- DSEN win rate: %s -> %s", metric_chr("dsen_win_rate", "threshold_15"), metric_chr("dsen_win_rate", "threshold_50")),
    sprintf("- Mean DSEN improvement: %s -> %s", metric_chr("dsen_mean_improvement", "threshold_15"), metric_chr("dsen_mean_improvement", "threshold_50")),
    sprintf("- Median DSEN improvement: %s -> %s", metric_chr("dsen_median_improvement", "threshold_15"), metric_chr("dsen_median_improvement", "threshold_50")),
    sprintf("- DSEN tissue recovery @50%%: %s -> %s", metric_chr("dsen_tissue_recovery_50pct", "threshold_15"), metric_chr("dsen_tissue_recovery_50pct", "threshold_50")),
    sprintf("- DSEN combined recovery @50%%: %s -> %s", metric_chr("dsen_combined_recovery_50pct", "threshold_15"), metric_chr("dsen_combined_recovery_50pct", "threshold_50")),
    "",
    paste0("**Recommendation:** ", primary_recommendation)
  )

  dir.create(dirname(decision_log_file), recursive = TRUE, showWarnings = FALSE)
  cat(paste(decision_block, collapse = "\n"), file = decision_log_file, append = TRUE)
}

cat("Wrote note:", note_file, "\n")
if (decision_exists) {
  cat("Decision log already contains Issue 8; left unchanged:", decision_log_file, "\n")
} else {
  cat("Appended Issue 8 to decision log:", decision_log_file, "\n")
}
