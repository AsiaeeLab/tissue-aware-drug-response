#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

TOTAL_CORES="${TOTAL_CORES:-20}"
JOB_CORES="${JOB_CORES:-10}"
SKIP_BOOTSTRAP="${SKIP_BOOTSTRAP:-TRUE}"
BOOTSTRAP_N="${BOOTSTRAP_N:-200}"
CHECKPOINT_EVERY="${CHECKPOINT_EVERY:-3}"
MAX_DRUGS="${MAX_DRUGS:-}"
TUNE_PER_SPLIT="${TUNE_PER_SPLIT:-TRUE}"
SKIP_PHASE0="${SKIP_PHASE0:-FALSE}"
SKIP_DSEN_ALL="${SKIP_DSEN_ALL:-FALSE}"

LOG_DIR="${LOG_DIR:-$ROOT_DIR/logs/comparison}"
mkdir -p "$LOG_DIR"

if [[ "$JOB_CORES" -gt 10 ]]; then
  echo "JOB_CORES must be <= 10 to follow experiment constraints."
  exit 1
fi

if [[ $((2 * JOB_CORES)) -gt "$TOTAL_CORES" ]]; then
  echo "With two concurrent jobs, 2*JOB_CORES exceeds TOTAL_CORES."
  exit 1
fi

COMMON_ENV=(
  "NCORES=$JOB_CORES"
  "SKIP_BOOTSTRAP=$SKIP_BOOTSTRAP"
  "BOOTSTRAP_N=$BOOTSTRAP_N"
  "CHECKPOINT_EVERY=$CHECKPOINT_EVERY"
)
if [[ -n "$MAX_DRUGS" ]]; then
  COMMON_ENV+=("MAX_DRUGS=$MAX_DRUGS")
fi

run_phase0() {
  echo "[Suite] Phase 0.1 save folds"
  Rscript methods/save_fold_assignments.R | tee "$LOG_DIR/phase0_save_folds.log"

  echo "[Suite] Phase 0.2 feature subsets"
  Rscript methods/create_feature_subsets.R | tee "$LOG_DIR/phase0_feature_subsets.log"

  echo "[Suite] Phase 0.3 fold reproducibility"
  Rscript methods/validate_fold_reproducibility.R | tee "$LOG_DIR/phase0_fold_reproducibility.log"
}

run_ablation_sync() {
  local run_id="$1"
  local feature_subset="$2"
  local log_file="$LOG_DIR/${run_id}.log"

  echo "[Suite] Running ablation $run_id ($feature_subset)"
  env "${COMMON_ENV[@]}" \
    RUN_ID="$run_id" \
    FEATURE_SUBSET="$feature_subset" \
    Rscript methods/run_dsen_ablation.R | tee "$log_file"
}

run_bg_job() {
  local job_name="$1"
  shift
  local log_file="$LOG_DIR/${job_name}.log"
  echo "[Suite] Launching $job_name in background"
  (
    "$@"
  ) >"$log_file" 2>&1 &
  RUN_BG_PID=$!
}

wait_for_any() {
  local -n pids_ref=$1
  local -n names_ref=$2

  while true; do
    for i in "${!pids_ref[@]}"; do
      local pid="${pids_ref[$i]}"
      if ! kill -0 "$pid" 2>/dev/null; then
        set +e
        wait "$pid"
        local rc=$?
        set -e
        if [[ "$rc" -ne 0 ]]; then
          echo "[Suite] Job failed: ${names_ref[$i]} (exit $rc)"
          return "$rc"
        fi
        echo "[Suite] Job complete: ${names_ref[$i]}"
        unset 'pids_ref[i]'
        unset 'names_ref[i]'
        pids_ref=("${pids_ref[@]}")
        names_ref=("${names_ref[@]}")
        return 0
      fi
    done
    sleep 20
  done
}

run_parallel_queue() {
  local -a names=()
  local -a pids=()

  local -a run_ids=(
    dsen_no_rppa
    dsen_expr_mut_cn
    dsen_expr_only
    dsen_mut_only
    dsen_rppa_only
    dsen_rppa_phospho
    dsen_rppa_total
    dsen_all_no_phospho
    dsen_all_no_total_rppa
    dsen_rppa_genes_expr
    dsen_rppa_genes_all_modalities
  )

  local -A subsets=(
    [dsen_no_rppa]=no_rppa
    [dsen_expr_mut_cn]=expr_mut_cn
    [dsen_expr_only]=expr_only
    [dsen_mut_only]=mut_only
    [dsen_rppa_only]=rppa_only
    [dsen_rppa_phospho]=rppa_phospho_only
    [dsen_rppa_total]=rppa_total_only
    [dsen_all_no_phospho]=all_no_phospho
    [dsen_all_no_total_rppa]=all_no_total_rppa
    [dsen_rppa_genes_expr]=rppa_genes_expr
    [dsen_rppa_genes_all_modalities]=rppa_genes_all_modalities
  )

  # Queue ablations, keeping max 2 concurrent jobs.
  for run_id in "${run_ids[@]}"; do
    while [[ "${#pids[@]}" -ge 2 ]]; do
      wait_for_any pids names
    done

    local subset="${subsets[$run_id]}"
    run_bg_job "$run_id" env "${COMMON_ENV[@]}" RUN_ID="$run_id" FEATURE_SUBSET="$subset" Rscript methods/run_dsen_ablation.R
    pids+=("$RUN_BG_PID")
    names+=("$run_id")
  done

  # TG-LASSO and LTO are queued under the same core budget.
  while [[ "${#pids[@]}" -ge 2 ]]; do
    wait_for_any pids names
  done
  run_bg_job "tglasso" env "${COMMON_ENV[@]}" OUT_DIR="results/comparison/tglasso" FEATURE_SUBSET="all_features" Rscript methods/run_tglasso.R
  pids+=("$RUN_BG_PID")
  names+=("tglasso")

  while [[ "${#pids[@]}" -ge 2 ]]; do
    wait_for_any pids names
  done
  run_bg_job "leave_tissue_out" env "${COMMON_ENV[@]}" "NCORES=$JOB_CORES" "TUNE_PER_SPLIT=$TUNE_PER_SPLIT" Rscript methods/run_leave_tissue_out.R
  pids+=("$RUN_BG_PID")
  names+=("leave_tissue_out")

  # Wait all remaining jobs.
  while [[ "${#pids[@]}" -gt 0 ]]; do
    wait_for_any pids names
  done
}

run_aggregation() {
  echo "[Suite] Phase 5 aggregation"
  Rscript methods/aggregate_comparison_results.R | tee "$LOG_DIR/phase5_aggregation.log"
}

main() {
  echo "[Suite] Starting DSEN comparison suite"
  echo "[Suite] TOTAL_CORES=$TOTAL_CORES JOB_CORES=$JOB_CORES SKIP_BOOTSTRAP=$SKIP_BOOTSTRAP"

  if [[ "$SKIP_PHASE0" != "TRUE" ]]; then
    run_phase0
  else
    echo "[Suite] Skipping Phase 0 (SKIP_PHASE0=TRUE)"
  fi

  # Required validation baseline run first unless explicitly skipped.
  if [[ "$SKIP_DSEN_ALL" != "TRUE" ]]; then
    run_ablation_sync dsen_all all_features
  else
    echo "[Suite] Skipping dsen_all baseline (SKIP_DSEN_ALL=TRUE)"
  fi

  # Remaining ablations + TG-LASSO + LTO under 20-core queue.
  run_parallel_queue

  run_aggregation
  echo "[Suite] Complete"
}

main "$@"
