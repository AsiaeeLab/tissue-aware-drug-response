#!/usr/bin/env Rscript
# Purpose: Write DSEN improvement note and append Issue 10 decision-log entry.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
})

comparison_csv <- Sys.getenv(
  "COMPARISON_CSV",
  here("results", "tables", "threshold_comparison_25_variants.csv")
)
penalty_cmp_csv <- Sys.getenv(
  "PENALTY_COMPARISON_CSV",
  here("results", "tables", "threshold_25_penalty_mode_comparison.csv")
)
winner_file <- Sys.getenv(
  "WINNER_FILE",
  here("results", "tables", "threshold_25_penalty_mode_winner.txt")
)
note_file <- Sys.getenv(
  "NOTE_FILE",
  here("doc", "notes", "2026-03-19_dsen_improvement_v2.md")
)
decision_log_file <- Sys.getenv(
  "DECISION_LOG_FILE",
  here("doc", "notes", "2026-03-15_tissue_recovery_decisions.md")
)

if (!file.exists(comparison_csv)) stop("Missing comparison CSV: ", comparison_csv)
if (!file.exists(penalty_cmp_csv)) stop("Missing penalty comparison CSV: ", penalty_cmp_csv)
if (!file.exists(winner_file)) stop("Missing winner file: ", winner_file)

cmp <- read_csv(comparison_csv, show_col_types = FALSE)
penalty_cmp <- read_csv(penalty_cmp_csv, show_col_types = FALSE)
winner_mode <- trimws(readLines(winner_file, warn = FALSE)[1])
if (!winner_mode %in% c("sample_size_mild", "uniform")) {
  stop("Winner mode file contained invalid value: ", winner_mode)
}

get_val <- function(metric, col) {
  row <- cmp %>% filter(.data$metric == .env$metric)
  if (nrow(row) == 0 || !col %in% names(row)) return(NA_character_)
  as.character(row[[col]][1])
}

to_num <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  token <- str_extract(x, "-?[0-9]+\\.?[0-9]*")
  suppressWarnings(as.numeric(token))
}

num_val <- function(metric, col) to_num(get_val(metric, col))

v1_col <- "threshold_25_v1"
v2_col <- "threshold_25_v2"
u_col <- "threshold_25_uniform"

win_v1 <- num_val("dsen_win_rate", v1_col)
win_v2 <- num_val("dsen_win_rate", v2_col)
win_u <- num_val("dsen_win_rate", u_col)

mean_v1 <- num_val("dsen_mean_improvement", v1_col)
mean_v2 <- num_val("dsen_mean_improvement", v2_col)
mean_u <- num_val("dsen_mean_improvement", u_col)

comb50_v1 <- num_val("dsen_combined_recovery_50pct", v1_col)
comb50_v2 <- num_val("dsen_combined_recovery_50pct", v2_col)

alpha_lt02_v2 <- get_val("alpha_lt_0_2_pct", v2_col)
alpha_eq005_v2 <- get_val("alpha_eq_0_05_pct", v2_col)
alpha_eq1_v2 <- get_val("alpha_eq_1_pct", v2_col)
lamrat_eq0125_v2 <- get_val("lamrat_eq_0_125_pct", v2_col)
lamrat_lt1_v2 <- get_val("lamrat_lt1_pct", v2_col)
lamrat_dist_v2 <- get_val("lamrat_distribution", v2_col)

pen_mild <- penalty_cmp %>% filter(mode == "sample_size_mild")
pen_uniform <- penalty_cmp %>% filter(mode == "uniform")

if (nrow(pen_mild) == 0 || nrow(pen_uniform) == 0) {
  stop("Penalty comparison table missing sample_size_mild or uniform rows.")
}

note_lines <- c(
  "# DSEN Improvement Rerun at Threshold 25 (Audit Fixes)",
  "",
  sprintf("Generated on %s.", Sys.Date()),
  "",
  "## Implemented Fixes",
  "- Fixed mutation regex so mutation subtype features (`_mut.mis/.tr/.hs/.hsrec`) are excluded from scaling.",
  "- Enforced paired tissue-stratified folds for EN and DSEN in production tuning.",
  "- Expanded alpha grid to `0.05, 0.1, 0.15, 0.2, 0.3, ..., 1.0`.",
  "- Extended lamrat grid low end to `0.125, 0.25, 0.5, 1, 2, 4`.",
  "- Removed global preprocessing leakage in bootstrap runner (raw `x` now passed into per-replicate preprocessing paths).",
  "- Passed tuned `lamrat_dsen` into production DSEN bootstrap calls.",
  "",
  "## Penalty-Mode Selection (Stage 1)",
  sprintf("- sample_size_mild: win rate %.1f%%, mean improvement %.2f%%", pen_mild$dsen_win_rate, pen_mild$dsen_mean_improvement),
  sprintf("- uniform: win rate %.1f%%, mean improvement %.2f%%", pen_uniform$dsen_win_rate, pen_uniform$dsen_mean_improvement),
  sprintf("- Winner for final bootstrap: `%s`", winner_mode),
  "",
  "## v1 vs v2 vs uniform",
  sprintf("- DSEN win rate: v1=%s, v2=%s, uniform=%s", get_val("dsen_win_rate", v1_col), get_val("dsen_win_rate", v2_col), get_val("dsen_win_rate", u_col)),
  sprintf("- Mean DSEN improvement: v1=%s, v2=%s, uniform=%s", get_val("dsen_mean_improvement", v1_col), get_val("dsen_mean_improvement", v2_col), get_val("dsen_mean_improvement", u_col)),
  sprintf("- DSEN combined target recovery @50%%: v1=%s, v2=%s", get_val("dsen_combined_recovery_50pct", v1_col), get_val("dsen_combined_recovery_50pct", v2_col)),
  "",
  "## Grid-Boundary Diagnostics (v2)",
  sprintf("- alpha < 0.2 frequency: %s", alpha_lt02_v2),
  sprintf("- alpha = 0.05 frequency: %s", alpha_eq005_v2),
  sprintf("- alpha = 1.0 frequency: %s", alpha_eq1_v2),
  sprintf("- lamrat = 0.125 frequency: %s", lamrat_eq0125_v2),
  sprintf("- lamrat < 1 frequency: %s", lamrat_lt1_v2),
  sprintf("- lamrat distribution: %s", lamrat_dist_v2),
  "",
  "## Interpretation",
  sprintf("- Win-rate delta v2-v1: %.2f pp", win_v2 - win_v1),
  sprintf("- Mean-improvement delta v2-v1: %.2f pp", mean_v2 - mean_v1),
  sprintf("- Combined-recovery delta v2-v1 (@50%%): %.2f pp", comb50_v2 - comb50_v1),
  sprintf("- Uniform-v2 Stage1 win-rate delta: %.2f pp", win_u - win_v2),
  sprintf("- Uniform-v2 Stage1 mean-improvement delta: %.2f pp", mean_u - mean_v2)
)

dir.create(dirname(note_file), recursive = TRUE, showWarnings = FALSE)
writeLines(note_lines, note_file)

issue_header <- "## Issue 10: DSEN audit fixes + threshold-25 rerun"
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
    "**What changed:**",
    "- Fixed mutation scaling regex to correctly exclude mutation subtype features.",
    "- Forced paired tissue-stratified folds for EN/DSEN tuning in production.",
    "- Expanded alpha grid (`0.05..1.0`) and lamrat grid (`0.125..4`).",
    "- Removed bootstrap global-preprocessing leakage by passing raw `x` into per-replicate preprocessing.",
    "- Added tuned `lamrat_dsen` passthrough into DSEN production bootstrap.",
    "",
    "**Penalty-mode test (Stage 1):**",
    sprintf("- sample_size_mild win rate / mean improvement: %.1f%% / %.2f%%", pen_mild$dsen_win_rate, pen_mild$dsen_mean_improvement),
    sprintf("- uniform win rate / mean improvement: %.1f%% / %.2f%%", pen_uniform$dsen_win_rate, pen_uniform$dsen_mean_improvement),
    sprintf("- Winner: `%s`", winner_mode),
    "",
    "**v1 vs v2 headline:**",
    sprintf("- DSEN win rate: %s -> %s", get_val("dsen_win_rate", v1_col), get_val("dsen_win_rate", v2_col)),
    sprintf("- Mean DSEN improvement: %s -> %s", get_val("dsen_mean_improvement", v1_col), get_val("dsen_mean_improvement", v2_col)),
    sprintf("- DSEN combined recovery @50%%: %s -> %s", get_val("dsen_combined_recovery_50pct", v1_col), get_val("dsen_combined_recovery_50pct", v2_col)),
    "",
    "**Grid-boundary behavior in v2:**",
    sprintf("- alpha < 0.2: %s", alpha_lt02_v2),
    sprintf("- lamrat = 0.125: %s", lamrat_eq0125_v2),
    sprintf("- lamrat < 1: %s", lamrat_lt1_v2),
    "",
    "**Recommendation:**",
    if (winner_mode == "uniform") {
      "- Use threshold_25_v2 with `uniform` tissue penalty and keep `sample_size_mild` as sensitivity."
    } else {
      "- Keep `sample_size_mild` as default at threshold_25_v2; report `uniform` as sensitivity."
    }
  )

  dir.create(dirname(decision_log_file), recursive = TRUE, showWarnings = FALSE)
  cat(paste(issue_block, collapse = "\n"), file = decision_log_file, append = TRUE)
}

cat("Wrote note:", note_file, "\n")
if (issue_exists) {
  cat("Issue 10 already exists in decision log; left unchanged:", decision_log_file, "\n")
} else {
  cat("Appended Issue 10 to decision log:", decision_log_file, "\n")
}
