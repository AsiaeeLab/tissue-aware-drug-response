#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

echo "======================================================================"
echo "Threshold-25 DSEN Pipeline (Corrected Lamrat Grid)"
echo "Started: $(date)"
echo "Repo: ${REPO_DIR}"
echo "======================================================================"

echo "[0/7] Ensure archive/output directories exist"
mkdir -p results/dsen_results/threshold_15 results/tables/threshold_15
mkdir -p results/dsen_results/threshold_50 results/tables/threshold_50
mkdir -p results/dsen_results/threshold_25 results/tables/threshold_25

echo "[1/7] Clean threshold-25 checkpoints"
rm -f results/dsen_results/threshold_25/dsen_production_checkpoint.Rda
rm -f results/dsen_results/threshold_25/dsen_bootstrap_checkpoint.Rda

echo "[2/7] Run DSEN production at threshold 25 (alpha+lamrat tuning)"
TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE=25 MIN_TISSUES=3 \
NCORES=15 CHECKPOINT_EVERY=15 SKIP_BOOTSTRAP=TRUE \
RESULTS_SUBDIR=threshold_25 \
Rscript methods/run_dsen_production.R

echo "[3/7] Build threshold-25 Stage 1 summary tables"
THRESHOLD_LABEL=threshold_25 \
Rscript methods/postprocess_dsen_production_stage1.R

echo "[4/7] Run DSEN bootstrap at threshold 25"
TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE=25 MIN_TISSUES=3 \
NCORES=15 CHECKPOINT_EVERY=20 \
STAGE1_SUMMARY_CSV=results/dsen_results/threshold_25/dsen_summary.csv \
BOOTSTRAP_RESULTS_FILE=results/dsen_results/threshold_25/dsen_bootstrap_results.Rda \
BOOTSTRAP_CHECKPOINT_FILE=results/dsen_results/threshold_25/dsen_bootstrap_checkpoint.Rda \
Rscript methods/run_dsen_bootstrap.R

echo "[5/7] Run threshold-25 recovery/Jaccard/difficulty + 15/25/50 comparison"
THRESHOLD15_STAGE1_RDA=results/dsen_results/threshold_15/dsen_production_results.Rda \
THRESHOLD15_BOOTSTRAP_RDA=results/dsen_results/threshold_15/dsen_bootstrap_results.Rda \
THRESHOLD25_STAGE1_RDA=results/dsen_results/threshold_25/dsen_production_results.Rda \
THRESHOLD25_BOOTSTRAP_RDA=results/dsen_results/threshold_25/dsen_bootstrap_results.Rda \
THRESHOLD50_STAGE1_RDA=results/dsen_results/threshold_50/dsen_production_results.Rda \
THRESHOLD50_BOOTSTRAP_RDA=results/dsen_results/threshold_50/dsen_bootstrap_results.Rda \
OUT_THRESHOLD25_TABLE_DIR=results/tables/threshold_25 \
OUT_COMPARISON_CSV=results/tables/threshold_comparison_15_25_50.csv \
OUT_THRESHOLD25_RECOVERY_RDA=results/dsen_results/threshold_25/tissue_target_recovery.Rda \
Rscript methods/run_threshold_comparison_15_25_50.R

echo "[6/7] Write threshold-25 narrative notes and decision-log update"
COMPARISON_CSV=results/tables/threshold_comparison_15_25_50.csv \
NOTE_FILE=doc/notes/2026-03-18_threshold_25_comparison.md \
DECISION_LOG_FILE=doc/notes/2026-03-15_tissue_recovery_decisions.md \
Rscript methods/write_threshold_25_notes.R

echo "======================================================================"
echo "Threshold-25 pipeline complete: $(date)"
echo "======================================================================"
