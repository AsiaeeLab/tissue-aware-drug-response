#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

NCORES="${NCORES:-15}"
CHECKPOINT_EVERY_STAGE1="${CHECKPOINT_EVERY_STAGE1:-60}"
CHECKPOINT_EVERY_BOOTSTRAP="${CHECKPOINT_EVERY_BOOTSTRAP:-20}"
MIN_SAMPLES_PER_TISSUE="${MIN_SAMPLES_PER_TISSUE:-25}"
MIN_TISSUES="${MIN_TISSUES:-3}"
GLOBAL_SEED="${GLOBAL_SEED:-42}"
BOOTSTRAP_N="${BOOTSTRAP_N:-200}"

V1_LABEL="${V1_LABEL:-threshold_25}"
V2_LABEL="${V2_LABEL:-threshold_25_v2}"
UNIFORM_LABEL="${UNIFORM_LABEL:-threshold_25_uniform}"
MILD_BACKUP_LABEL="${MILD_BACKUP_LABEL:-threshold_25_v2_sample_size_mild}"

LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "${LOG_DIR}"
mkdir -p "results/dsen_results/${V2_LABEL}" "results/tables/${V2_LABEL}"
mkdir -p "results/dsen_results/${UNIFORM_LABEL}" "results/tables/${UNIFORM_LABEL}"

run_stage1() {
  local label="$1"
  local mode="$2"

  echo "[Stage1] ${label} (${mode})"
  rm -f "results/dsen_results/${label}/dsen_production_checkpoint.Rda"

  TEST_MODE=FALSE SKIP_BOOTSTRAP=TRUE \
  MIN_SAMPLES_PER_TISSUE="${MIN_SAMPLES_PER_TISSUE}" MIN_TISSUES="${MIN_TISSUES}" \
  TISSUE_PENALTY_MODE="${mode}" \
  NCORES="${NCORES}" CHECKPOINT_EVERY="${CHECKPOINT_EVERY_STAGE1}" \
  GLOBAL_SEED="${GLOBAL_SEED}" RESULTS_SUBDIR="${label}" \
  Rscript methods/run_dsen_production.R \
    2>&1 | tee "${LOG_DIR}/${label}_production.log"

  THRESHOLD_LABEL="${label}" \
  Rscript methods/postprocess_dsen_production_stage1.R \
    2>&1 | tee "${LOG_DIR}/${label}_postprocess.log"
}

echo "======================================================================"
echo "DSEN Improvement Pipeline"
echo "Started: $(date)"
echo "Repo: ${REPO_DIR}"
echo "Cores: ${NCORES}"
echo "======================================================================"

run_stage1 "${V2_LABEL}" "sample_size_mild"
run_stage1 "${UNIFORM_LABEL}" "uniform"

echo "[Select] Compare penalty modes"
MILD_LABEL="${V2_LABEL}" \
UNIFORM_LABEL="${UNIFORM_LABEL}" \
OUT_DIR="results/tables" \
OUT_COMPARISON_CSV="results/tables/threshold_25_penalty_mode_comparison.csv" \
OUT_WINNER_FILE="results/tables/threshold_25_penalty_mode_winner.txt" \
Rscript methods/select_threshold_25_penalty_winner.R \
  2>&1 | tee "${LOG_DIR}/threshold_25_penalty_selection.log"

winner_mode="$(head -n 1 results/tables/threshold_25_penalty_mode_winner.txt | tr -d '[:space:]')"
echo "Winner mode: ${winner_mode}"

if [[ "${winner_mode}" == "uniform" ]]; then
  echo "[Switch] Uniform won: preserving sample_size_mild outputs under ${MILD_BACKUP_LABEL}"
  mkdir -p "results/dsen_results/${MILD_BACKUP_LABEL}" "results/tables/${MILD_BACKUP_LABEL}"
  cp -a "results/dsen_results/${V2_LABEL}/." "results/dsen_results/${MILD_BACKUP_LABEL}/"
  cp -a "results/tables/${V2_LABEL}/." "results/tables/${MILD_BACKUP_LABEL}/"

  echo "[Switch] Copying uniform Stage-1 artifacts into ${V2_LABEL}"
  cp -f "results/dsen_results/${UNIFORM_LABEL}/dsen_production_results.Rda" "results/dsen_results/${V2_LABEL}/dsen_production_results.Rda"
  cp -f "results/dsen_results/${UNIFORM_LABEL}/dsen_summary.csv" "results/dsen_results/${V2_LABEL}/dsen_summary.csv"
  cp -f "results/tables/${UNIFORM_LABEL}/dsen_full_results.csv" "results/tables/${V2_LABEL}/dsen_full_results.csv"
  cp -f "results/tables/${UNIFORM_LABEL}/dsen_analysis_summary.csv" "results/tables/${V2_LABEL}/dsen_analysis_summary.csv"
fi

echo "[Bootstrap] Winner mode into ${V2_LABEL}"
rm -f "results/dsen_results/${V2_LABEL}/dsen_bootstrap_checkpoint.Rda"

TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE="${MIN_SAMPLES_PER_TISSUE}" MIN_TISSUES="${MIN_TISSUES}" \
TISSUE_PENALTY_MODE="${winner_mode}" GLOBAL_SEED="${GLOBAL_SEED}" \
BOOTSTRAP_N="${BOOTSTRAP_N}" \
NCORES="${NCORES}" CHECKPOINT_EVERY="${CHECKPOINT_EVERY_BOOTSTRAP}" \
STAGE1_SUMMARY_CSV="results/dsen_results/${V2_LABEL}/dsen_summary.csv" \
BOOTSTRAP_RESULTS_FILE="results/dsen_results/${V2_LABEL}/dsen_bootstrap_results.Rda" \
BOOTSTRAP_CHECKPOINT_FILE="results/dsen_results/${V2_LABEL}/dsen_bootstrap_checkpoint.Rda" \
Rscript methods/run_dsen_bootstrap.R \
  2>&1 | tee "${LOG_DIR}/${V2_LABEL}_bootstrap.log"

echo "[Compare] Build v1-v2-uniform tables"
V1_LABEL="${V1_LABEL}" V2_LABEL="${V2_LABEL}" UNIFORM_LABEL="${UNIFORM_LABEL}" \
OUT_V2_TABLE_DIR="results/tables/${V2_LABEL}" \
OUT_V2_RECOVERY_RDA="results/dsen_results/${V2_LABEL}/tissue_target_recovery.Rda" \
OUT_COMPARISON_CSV="results/tables/threshold_comparison_25_variants.csv" \
Rscript methods/run_threshold_comparison_25_variants.R \
  2>&1 | tee "${LOG_DIR}/threshold_25_variant_comparison.log"

echo "[Notes] Write Issue 10 and summary note"
COMPARISON_CSV="results/tables/threshold_comparison_25_variants.csv" \
PENALTY_COMPARISON_CSV="results/tables/threshold_25_penalty_mode_comparison.csv" \
WINNER_FILE="results/tables/threshold_25_penalty_mode_winner.txt" \
NOTE_FILE="doc/notes/2026-03-19_dsen_improvement_v2.md" \
DECISION_LOG_FILE="doc/notes/2026-03-15_tissue_recovery_decisions.md" \
Rscript methods/write_dsen_improvement_notes.R \
  2>&1 | tee "${LOG_DIR}/threshold_25_issue10_notes.log"

echo "======================================================================"
echo "DSEN improvement pipeline complete: $(date)"
echo "======================================================================"
