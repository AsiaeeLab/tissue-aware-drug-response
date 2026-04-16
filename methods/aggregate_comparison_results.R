#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(purrr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "comparison_common.R"))

# =============================================================================
# CONFIG
# =============================================================================
out_dir <- ensure_comparison_dir()
ablation_dir <- file.path(out_dir, "ablation")
tglasso_dir <- file.path(out_dir, "tglasso")

feature_subset_file <- file.path(out_dir, "feature_subsets.rds")
if (!file.exists(feature_subset_file)) {
  stop("Missing feature subsets file: ", feature_subset_file)
}
feature_bundle <- readRDS(feature_subset_file)
feature_subsets <- if (is.list(feature_bundle) && "subsets" %in% names(feature_bundle)) {
  feature_bundle$subsets
} else {
  feature_bundle
}

thresholds <- c(0.005, 0.05, 0.1, 0.25, 0.5, 0.75)

run_to_subset <- c(
  dsen_all = "all_features",
  dsen_no_rppa = "no_rppa",
  dsen_expr_mut_cn = "expr_mut_cn",
  dsen_expr_only = "expr_only",
  dsen_mut_only = "mut_only",
  dsen_rppa_only = "rppa_only",
  dsen_rppa_phospho = "rppa_phospho_only",
  dsen_rppa_total = "rppa_total_only",
  dsen_all_no_phospho = "all_no_phospho",
  dsen_all_no_total_rppa = "all_no_total_rppa",
  dsen_rppa_genes_expr = "rppa_genes_expr",
  dsen_rppa_genes_all_modalities = "rppa_genes_all_modalities"
)
run_ids <- names(run_to_subset)

# =============================================================================
# HELPERS
# =============================================================================
load_run_summary <- function(run_id) {
  csv <- file.path(ablation_dir, run_id, "dsen_summary.csv")
  if (!file.exists(csv)) return(NULL)
  read.csv(csv, stringsAsFactors = FALSE)
}

load_run_bootstrap <- function(run_id) {
  rda <- file.path(ablation_dir, run_id, "dsen_bootstrap_results.Rda")
  if (!file.exists(rda)) return(NULL)
  env <- new.env(parent = emptyenv())
  load(rda, envir = env)
  if (exists("all_results", envir = env, inherits = FALSE)) {
    return(get("all_results", envir = env))
  }
  NULL
}

feature_info <- build_feature_gene_map(paths)
rppa_map <- parse_rppa_mapping()
rppa_genes <- unique(rppa_map$hgnc_symbol)
rppa_phospho <- rppa_map$rppa_feature[rppa_map$is_phospho]
rppa_features <- rppa_map$rppa_feature

feature_info <- feature_info %>%
  mutate(
    is_expr = grepl("_Expr$", feature),
    is_mut = grepl("_mut\\.", feature),
    is_cn = grepl("_CN$", feature),
    is_rppa = feature %in% rppa_features,
    is_rppa_phospho = feature %in% rppa_phospho
  )

feature_to_gene <- setNames(feature_info$gene, feature_info$feature)

targets <- standardize_targets(feature_info)

extract_dsen_selected <- function(res, threshold) {
  if (is.null(res) || is.null(res$dsen_bootstrap)) return(character(0))
  bs <- res$dsen_bootstrap
  if (is.null(bs$cnt) || is.null(bs$nboot) || bs$nboot <= 0) return(character(0))

  shared_freq <- bs$cnt / bs$nboot
  shared_sel <- names(shared_freq)[shared_freq >= threshold]

  tissue_sel <- character(0)
  if (!is.null(bs$tissue_cnt)) {
    tissue_freq <- bs$tissue_cnt / bs$nboot
    tissue_sel <- rownames(tissue_freq)[apply(tissue_freq >= threshold, 1, any)]
  }

  unique(c(shared_sel, tissue_sel))
}

mode_feature_string <- function(sel_features, target_gene, mode = c("rppa", "expr", "mut", "phospho")) {
  mode <- match.arg(mode)
  if (length(sel_features) == 0) return("")

  idx <- sel_features
  if (mode == "rppa") {
    idx <- idx[idx %in% feature_info$feature[feature_info$is_rppa]]
  } else if (mode == "expr") {
    idx <- idx[idx %in% feature_info$feature[feature_info$is_expr]]
  } else if (mode == "mut") {
    idx <- idx[idx %in% feature_info$feature[feature_info$is_mut]]
  } else if (mode == "phospho") {
    idx <- idx[idx %in% feature_info$feature[feature_info$is_rppa_phospho]]
  }

  if (length(idx) == 0) return("")
  idx <- idx[feature_to_gene[idx] == target_gene]
  if (length(idx) == 0) return("")
  paste(unique(idx), collapse = ";")
}

compute_target_recovery_by_modality <- function(run_id, run_subset, bootstrap_results) {
  if (is.null(bootstrap_results)) {
    return(list(per_target = tibble(), per_run_summary = tibble()))
  }

  subset_features <- feature_subsets[[run_subset]]
  subset_genes <- unique(feature_to_gene[subset_features])
  subset_genes <- subset_genes[!is.na(subset_genes)]

  per_target <- map_dfr(thresholds, function(thr) {
    map_dfr(names(bootstrap_results), function(drug_nm) {
      res <- bootstrap_results[[drug_nm]]
      if (is.null(res) || is.null(res$status) || res$status != "success") return(NULL)

      sel_features <- extract_dsen_selected(res, thr)
      sel_features <- intersect(sel_features, subset_features)

      drug_targets <- targets %>%
        filter(.data$drug_name == .env$drug_nm, target_in_features, target_gene %in% subset_genes) %>%
        distinct(target_gene)

      if (nrow(drug_targets) == 0) return(NULL)

      map_dfr(drug_targets$target_gene, function(g) {
        rppa_feat <- mode_feature_string(sel_features, g, "rppa")
        expr_feat <- mode_feature_string(sel_features, g, "expr")
        mut_feat <- mode_feature_string(sel_features, g, "mut")
        phospho_feat <- mode_feature_string(sel_features, g, "phospho")

        tibble(
          run_id = run_id,
          threshold = thr,
          drug = drug_nm,
          target_gene = g,
          in_rppa_panel = g %in% rppa_genes,
          recovered_any = (rppa_feat != "") || (expr_feat != "") || (mut_feat != "") ||
            any(feature_to_gene[sel_features] == g, na.rm = TRUE),
          recovered_via_rppa = rppa_feat != "",
          recovered_via_expr = expr_feat != "",
          recovered_via_mut = mut_feat != "",
          recovered_via_phospho_rppa = phospho_feat != "",
          rppa_feature_name = rppa_feat,
          expression_feature_name = expr_feat,
          mutation_feature_names = mut_feat,
          phospho_rppa_feature_name = phospho_feat
        )
      })
    })
  })

  per_run_summary <- per_target %>%
    group_by(run_id, threshold) %>%
    summarize(
      n_drugs = n_distinct(drug),
      n_drugs_any_recovered = n_distinct(drug[recovered_any]),
      pct_drugs_any_recovered = 100 * n_drugs_any_recovered / n_drugs,
      n_targets = n(),
      n_targets_recovered_any = sum(recovered_any),
      pct_targets_recovered_any = 100 * mean(recovered_any),
      pct_targets_recovered_via_rppa = 100 * mean(recovered_via_rppa),
      pct_targets_recovered_via_expr = 100 * mean(recovered_via_expr),
      pct_targets_recovered_via_mut = 100 * mean(recovered_via_mut),
      pct_targets_recovered_via_phospho = 100 * mean(recovered_via_phospho_rppa),
      .groups = "drop"
    )

  list(per_target = per_target, per_run_summary = per_run_summary)
}

# =============================================================================
# LOAD RUN ARTIFACTS
# =============================================================================
run_summaries <- setNames(lapply(run_ids, load_run_summary), run_ids)
run_bootstrap <- setNames(lapply(run_ids, load_run_bootstrap), run_ids)

# =============================================================================
# PHASE 1.4: TARGET RECOVERY BY MODALITY
# =============================================================================
recovery_parts <- map2(run_ids, run_to_subset[run_ids], function(run_id, subset_name) {
  compute_target_recovery_by_modality(run_id, subset_name, run_bootstrap[[run_id]])
})
names(recovery_parts) <- run_ids

recovery_per_target <- bind_rows(lapply(recovery_parts, `[[`, "per_target"))
recovery_summary <- bind_rows(lapply(recovery_parts, `[[`, "per_run_summary"))

write.csv(recovery_summary, file.path(out_dir, "target_recovery_by_modality.csv"), row.names = FALSE)
saveRDS(
  list(
    per_target = recovery_per_target,
    summary = recovery_summary,
    thresholds = thresholds
  ),
  file.path(out_dir, "target_recovery_by_modality.rds")
)

# =============================================================================
# CHECKPOINT AFTER PHASE 1
# =============================================================================
ablation_summary <- map_dfr(run_ids, function(run_id) {
  df <- run_summaries[[run_id]]
  if (is.null(df) || nrow(df) == 0) {
    return(tibble(
      run_id = run_id,
      n_features_available = length(feature_subsets[[run_to_subset[[run_id]]]]),
      n_features_after_screening = NA_real_,
      dsen_win_rate = NA_real_,
      mean_improvement = NA_real_,
      median_improvement = NA_real_,
      n_targets_recovered_at_0.5 = NA_real_
    ))
  }

  success <- df %>% filter(status == "success")
  rec_05 <- recovery_summary %>% filter(.data$run_id == .env$run_id, abs(threshold - 0.5) < 1e-12)

  tibble(
    run_id = run_id,
    n_features_available = if ("n_features_available" %in% names(success) && nrow(success) > 0) {
      round(mean(success$n_features_available, na.rm = TRUE))
    } else {
      length(feature_subsets[[run_to_subset[[run_id]]]])
    },
    n_features_after_screening = if ("mean_features_after_screening_dsen" %in% names(success)) {
      mean(success$mean_features_after_screening_dsen, na.rm = TRUE)
    } else {
      NA_real_
    },
    dsen_win_rate = 100 * mean(success$pct_improvement > 0, na.rm = TRUE),
    mean_improvement = mean(success$pct_improvement, na.rm = TRUE),
    median_improvement = median(success$pct_improvement, na.rm = TRUE),
    n_targets_recovered_at_0.5 = if (nrow(rec_05) == 0) NA_real_ else rec_05$n_drugs_any_recovered[1]
  )
})

write.csv(ablation_summary, file.path(out_dir, "ablation_summary.csv"), row.names = FALSE)

# =============================================================================
# PHASE 2: METHOD COMPARISON
# =============================================================================
dsen_all <- run_summaries[["dsen_all"]]
tglasso_summary <- {
  csv <- file.path(tglasso_dir, "tglasso_summary.csv")
  if (file.exists(csv)) read.csv(csv, stringsAsFactors = FALSE) else NULL
}

method_comparison <- NULL
if (!is.null(dsen_all) && !is.null(tglasso_summary)) {
  method_comparison <- dsen_all %>%
    select(drug, mse_en, mse_dsen) %>%
    left_join(tglasso_summary %>% select(drug, mse_tglasso), by = "drug") %>%
    mutate(
      pct_improvement_dsen = (mse_en - mse_dsen) / mse_en * 100,
      pct_improvement_tglasso = (mse_en - mse_tglasso) / mse_en * 100,
      winner = pmap_chr(
        list(mse_en, mse_dsen, mse_tglasso),
        function(en, dsen, tg) {
          vals <- c(EN = en, DSEN = dsen, `TG-LASSO` = tg)
          vals <- vals[is.finite(vals)]
          if (length(vals) == 0) return(NA_character_)
          names(vals)[which.min(vals)]
        }
      )
    )

  write.csv(method_comparison, file.path(out_dir, "method_comparison.csv"), row.names = FALSE)
}

# =============================================================================
# PHASE 5.1: MASTER TABLE
# =============================================================================
get_mse_col <- function(run_id) {
  df <- run_summaries[[run_id]]
  if (is.null(df)) return(NULL)
  df %>% select(drug, mse_dsen)
}

master <- NULL
if (!is.null(dsen_all)) {
  master <- dsen_all %>%
    select(drug, mse_en_all = mse_en, mse_dsen_all = mse_dsen)

  add_run <- function(master_df, run_id, col_name) {
    df <- get_mse_col(run_id)
    if (is.null(df)) {
      master_df[[col_name]] <- NA_real_
      return(master_df)
    }
    names(df)[2] <- col_name
    left_join(master_df, df, by = "drug")
  }

  master <- add_run(master, "dsen_no_rppa", "mse_dsen_no_rppa")
  master <- add_run(master, "dsen_expr_mut_cn", "mse_dsen_expr_mut_cn")
  master <- add_run(master, "dsen_rppa_only", "mse_dsen_rppa_only")
  master <- add_run(master, "dsen_mut_only", "mse_dsen_mut_only")

  if (!is.null(tglasso_summary)) {
    master <- left_join(master, tglasso_summary %>% select(drug, mse_tglasso), by = "drug")
  } else {
    master$mse_tglasso <- NA_real_
  }

  master <- master %>%
    mutate(
      pct_dsen_vs_en = (mse_en_all - mse_dsen_all) / mse_en_all * 100,
      pct_rppa_contribution = (mse_dsen_no_rppa - mse_dsen_all) / mse_dsen_no_rppa * 100,
      winner_overall = pmap_chr(
        list(mse_en_all, mse_dsen_all, mse_tglasso),
        function(en, dsen, tg) {
          vals <- c(EN = en, DSEN_all = dsen, `TG-LASSO` = tg)
          vals <- vals[is.finite(vals)]
          if (length(vals) == 0) return(NA_character_)
          names(vals)[which.min(vals)]
        }
      )
    )

  write.csv(master, file.path(out_dir, "master_comparison.csv"), row.names = FALSE)
}

# =============================================================================
# PHASE 5.2: FIGURE DATA TABLE
# =============================================================================
ablation_figure <- ablation_summary %>%
  mutate(
    feature_subset = run_to_subset[run_id],
    n_features = n_features_available,
    dsen_win_rate_vs_en = dsen_win_rate,
    n_targets_recovered = n_targets_recovered_at_0.5
  ) %>%
  select(
    feature_subset,
    n_features,
    dsen_win_rate_vs_en,
    mean_improvement,
    median_improvement,
    n_targets_recovered
  )

write.csv(ablation_figure, file.path(out_dir, "ablation_figure_data.csv"), row.names = FALSE)

# =============================================================================
# PHASE 5.3: KEY BIO TABLES
# =============================================================================
table_a <- recovery_per_target %>%
  filter(abs(threshold - 0.5) < 1e-12, in_rppa_panel) %>%
  transmute(
    Run = run_id,
    Drug = drug,
    `Target gene` = target_gene,
    `In RPPA panel?` = in_rppa_panel,
    `Recovered via RPPA?` = recovered_via_rppa,
    `Recovered via Expr?` = recovered_via_expr,
    `Recovered via Mut?` = recovered_via_mut,
    `RPPA feature name` = rppa_feature_name,
    `Phospho?` = recovered_via_phospho_rppa
  )

write.csv(table_a, file.path(out_dir, "tableA_rppa_vs_expression_target_recovery.csv"), row.names = FALSE)

table_b <- tibble(
  Metric = c("DSEN win rate", "Mean improvement", "Targets recovered"),
  `All RPPA` = c(
    ablation_summary$dsen_win_rate[ablation_summary$run_id == "dsen_rppa_only"],
    ablation_summary$mean_improvement[ablation_summary$run_id == "dsen_rppa_only"],
    ablation_summary$n_targets_recovered_at_0.5[ablation_summary$run_id == "dsen_rppa_only"]
  ),
  `Total only` = c(
    ablation_summary$dsen_win_rate[ablation_summary$run_id == "dsen_rppa_total"],
    ablation_summary$mean_improvement[ablation_summary$run_id == "dsen_rppa_total"],
    ablation_summary$n_targets_recovered_at_0.5[ablation_summary$run_id == "dsen_rppa_total"]
  ),
  `Phospho only` = c(
    ablation_summary$dsen_win_rate[ablation_summary$run_id == "dsen_rppa_phospho"],
    ablation_summary$mean_improvement[ablation_summary$run_id == "dsen_rppa_phospho"],
    ablation_summary$n_targets_recovered_at_0.5[ablation_summary$run_id == "dsen_rppa_phospho"]
  ),
  `All minus phospho` = c(
    ablation_summary$dsen_win_rate[ablation_summary$run_id == "dsen_all_no_phospho"],
    ablation_summary$mean_improvement[ablation_summary$run_id == "dsen_all_no_phospho"],
    ablation_summary$n_targets_recovered_at_0.5[ablation_summary$run_id == "dsen_all_no_phospho"]
  ),
  `All minus total` = c(
    ablation_summary$dsen_win_rate[ablation_summary$run_id == "dsen_all_no_total_rppa"],
    ablation_summary$mean_improvement[ablation_summary$run_id == "dsen_all_no_total_rppa"],
    ablation_summary$n_targets_recovered_at_0.5[ablation_summary$run_id == "dsen_all_no_total_rppa"]
  )
)

write.csv(table_b, file.path(out_dir, "tableB_phospho_vs_total_contributions.csv"), row.names = FALSE)

# =============================================================================
# PROGRESS LOG UPDATE
# =============================================================================
append_progress(
  phase = "Phase 5 - Aggregation",
  status = "completed",
  details = c(
    sprintf("Ablation runs detected: %d/%d", sum(!sapply(run_summaries, is.null)), length(run_summaries)),
    sprintf("Target recovery rows: %d", nrow(recovery_per_target)),
    sprintf("Method comparison available: %s", ifelse(is.null(method_comparison), "FALSE", "TRUE")),
    sprintf("Master comparison available: %s", ifelse(is.null(master), "FALSE", "TRUE"))
  )
)

cat("Saved aggregation outputs under:", out_dir, "\n")
if (nrow(ablation_summary) > 0) {
  cat("\nAblation summary:\n")
  print(ablation_summary)
}
