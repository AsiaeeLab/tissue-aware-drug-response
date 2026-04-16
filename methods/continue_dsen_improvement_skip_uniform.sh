#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

NCORES="${NCORES:-15}"
CHECKPOINT_EVERY_BOOTSTRAP="${CHECKPOINT_EVERY_BOOTSTRAP:-20}"
MIN_SAMPLES_PER_TISSUE="${MIN_SAMPLES_PER_TISSUE:-25}"
MIN_TISSUES="${MIN_TISSUES:-3}"
GLOBAL_SEED="${GLOBAL_SEED:-42}"
BOOTSTRAP_N="${BOOTSTRAP_N:-200}"

V1_LABEL="${V1_LABEL:-threshold_25}"
V2_LABEL="${V2_LABEL:-threshold_25_v2}"

WAIT_FOR_STAGE1="${WAIT_FOR_STAGE1:-TRUE}"
WAIT_SECONDS="${WAIT_SECONDS:-300}"

LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "${LOG_DIR}"
mkdir -p "results/dsen_results/${V2_LABEL}" "results/tables/${V2_LABEL}"

stage1_rda="results/dsen_results/${V2_LABEL}/dsen_production_results.Rda"
stage1_csv="results/dsen_results/${V2_LABEL}/dsen_summary.csv"

if [[ "${WAIT_FOR_STAGE1}" == "TRUE" ]]; then
  echo "[Wait] Waiting for Phase 2 Stage-1 output: ${stage1_rda}"
  while [[ ! -f "${stage1_rda}" ]]; do
    echo "[Wait] $(date): Stage-1 not complete yet; sleeping ${WAIT_SECONDS}s"
    sleep "${WAIT_SECONDS}"
  done
fi

if [[ ! -f "${stage1_rda}" ]]; then
  echo "[Error] Missing Stage-1 RDA: ${stage1_rda}" >&2
  exit 1
fi

if [[ ! -f "${stage1_csv}" ]]; then
  echo "[Error] Missing Stage-1 summary CSV: ${stage1_csv}" >&2
  exit 1
fi

echo "[Phase 2.5] Postprocess Stage-1 tables for ${V2_LABEL}"
THRESHOLD_LABEL="${V2_LABEL}" \
Rscript methods/postprocess_dsen_production_stage1.R \
  2>&1 | tee "${LOG_DIR}/${V2_LABEL}_postprocess.log"

echo "[Phase 4] Bootstrap for ${V2_LABEL} (sample_size_mild only)"
rm -f "results/dsen_results/${V2_LABEL}/dsen_bootstrap_checkpoint.Rda"

TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE="${MIN_SAMPLES_PER_TISSUE}" MIN_TISSUES="${MIN_TISSUES}" \
TISSUE_PENALTY_MODE="sample_size_mild" GLOBAL_SEED="${GLOBAL_SEED}" \
BOOTSTRAP_N="${BOOTSTRAP_N}" \
NCORES="${NCORES}" CHECKPOINT_EVERY="${CHECKPOINT_EVERY_BOOTSTRAP}" \
STAGE1_SUMMARY_CSV="${stage1_csv}" \
BOOTSTRAP_RESULTS_FILE="results/dsen_results/${V2_LABEL}/dsen_bootstrap_results.Rda" \
BOOTSTRAP_CHECKPOINT_FILE="results/dsen_results/${V2_LABEL}/dsen_bootstrap_checkpoint.Rda" \
Rscript methods/run_dsen_bootstrap.R \
  2>&1 | tee "${LOG_DIR}/${V2_LABEL}_bootstrap.log"

echo "[Phase 5] Target recovery + v1 vs v2 comparison (uniform skipped)"
V1_LABEL="${V1_LABEL}" V2_LABEL="${V2_LABEL}" INCLUDE_UNIFORM=FALSE \
OUT_V2_TABLE_DIR="results/tables/${V2_LABEL}" \
OUT_V2_RECOVERY_RDA="results/dsen_results/${V2_LABEL}/tissue_target_recovery.Rda" \
OUT_COMPARISON_CSV="results/tables/threshold_comparison_25_v1_v2.csv" \
Rscript methods/run_threshold_comparison_25_variants.R \
  2>&1 | tee "${LOG_DIR}/threshold_25_v1_v2_comparison.log"

echo "[Phase 5] Quick v2-v1 Stage-1 gate metrics"
Rscript -e '
  suppressPackageStartupMessages(library(readr));
  v1 <- read_csv("results/dsen_results/'"${V1_LABEL}"'/dsen_summary.csv", show_col_types = FALSE);
  v2 <- read_csv("results/dsen_results/'"${V2_LABEL}"'/dsen_summary.csv", show_col_types = FALSE);
  v1 <- v1[v1$status == "success", ];
  v2 <- v2[v2$status == "success", ];
  win1 <- 100 * mean(v1$pct_improvement > 0, na.rm = TRUE);
  win2 <- 100 * mean(v2$pct_improvement > 0, na.rm = TRUE);
  mean1 <- mean(v1$pct_improvement, na.rm = TRUE);
  mean2 <- mean(v2$pct_improvement, na.rm = TRUE);
  out <- data.frame(
    metric = c("dsen_win_rate", "dsen_mean_improvement", "delta_win_rate_pp", "delta_mean_improvement_pp"),
    threshold_25_v1 = c(sprintf("%.1f%%", win1), sprintf("%.2f%%", mean1), NA, NA),
    threshold_25_v2 = c(sprintf("%.1f%%", win2), sprintf("%.2f%%", mean2), NA, NA),
    delta_v2_minus_v1 = c(NA, NA, sprintf("%.2f", win2 - win1), sprintf("%.2f", mean2 - mean1))
  );
  readr::write_csv(out, "results/tables/threshold_25_v2_vs_v1_gate_metrics.csv");
  cat("Wrote results/tables/threshold_25_v2_vs_v1_gate_metrics.csv\n");
'

echo "Completed skip-uniform continuation: $(date)"
