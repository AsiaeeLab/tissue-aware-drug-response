#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_DIR}"

mkdir -p results/dsen_results/threshold_50 results/tables/threshold_50

echo "======================================================================"
echo "Threshold-50 DSEN Pipeline"
echo "Started: $(date)"
echo "Repo: ${REPO_DIR}"
echo "======================================================================"

echo "[1/6] Archive current threshold-15 baseline files"
mkdir -p results/dsen_results/threshold_15 results/tables/threshold_15
cp -f results/dsen_results/dsen_production_results.Rda results/dsen_results/threshold_15/
cp -f results/dsen_results/dsen_bootstrap_results.Rda results/dsen_results/threshold_15/
cp -f results/dsen_results/dsen_summary.csv results/dsen_results/threshold_15/
cp -f results/tables/dsen_full_results.csv results/tables/threshold_15/
cp -f results/tables/dsen_analysis_summary.csv results/tables/threshold_15/
cp -f results/tables/feature_stability_summary.csv results/tables/threshold_15/
cp -f results/tables/target_recovery_extended.csv results/tables/threshold_15/
cp -f results/tables/target_recovery_extended_summary.csv results/tables/threshold_15/
cp -f results/tables/rerun_metric_summary.csv results/tables/threshold_15/
cp -f results/tables/publication_summary_stats.csv results/tables/threshold_15/

echo "[2/6] Run Stage 1 DSEN analysis (threshold 50) with checkpoint resume"
TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE=50 NCORES=15 CHECKPOINT_EVERY=15 \
  Rscript methods/run_dsen_analysis.R

echo "[3/6] Postprocess Stage 1 outputs into threshold_50 folders"
THRESHOLD_LABEL=threshold_50 \
  Rscript methods/postprocess_dsen_stage1_results.R

echo "[4/6] Run DSEN bootstrap with threshold-50 tissue filtering"
TEST_MODE=FALSE MIN_SAMPLES_PER_TISSUE=50 MIN_TISSUES=3 NCORES=15 \
  STAGE1_SUMMARY_CSV=results/dsen_results/threshold_50/dsen_summary.csv \
  BOOTSTRAP_RESULTS_FILE=results/dsen_results/threshold_50/dsen_bootstrap_results.Rda \
  BOOTSTRAP_CHECKPOINT_FILE=results/dsen_results/threshold_50/dsen_bootstrap_checkpoint.Rda \
  Rscript methods/run_dsen_bootstrap.R

echo "[5/6] Compute target recovery, Jaccard, and threshold comparison tables"
THRESHOLD50_STAGE1_RDA=results/dsen_results/threshold_50/dsen_production_results.Rda \
THRESHOLD50_BOOTSTRAP_RDA=results/dsen_results/threshold_50/dsen_bootstrap_results.Rda \
THRESHOLD15_STAGE1_RDA=results/dsen_results/threshold_15/dsen_production_results.Rda \
THRESHOLD15_BOOTSTRAP_RDA=results/dsen_results/threshold_15/dsen_bootstrap_results.Rda \
OUT_THRESHOLD50_TABLE_DIR=results/tables/threshold_50 \
OUT_COMPARISON_CSV=results/tables/threshold_comparison_15_vs_50.csv \
OUT_THRESHOLD50_RECOVERY_RDA=results/dsen_results/threshold_50/tissue_target_recovery.Rda \
  Rscript methods/run_threshold_comparison_analysis.R

echo "[6/6] Write narrative summary and decision-log update"
COMPARISON_CSV=results/tables/threshold_comparison_15_vs_50.csv \
NOTE_FILE=doc/notes/2026-03-18_threshold_comparison.md \
DECISION_LOG_FILE=doc/notes/2026-03-15_tissue_recovery_decisions.md \
  Rscript methods/write_threshold_comparison_notes.R

echo "======================================================================"
echo "Pipeline complete: $(date)"
echo "======================================================================"
