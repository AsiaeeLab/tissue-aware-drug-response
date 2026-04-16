#!/usr/bin/env Rscript
# Purpose: Compute threshold-50 target recovery + EN/DSEN Jaccard metrics and
#          build a side-by-side threshold_15 vs threshold_50 summary table.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})

source(here("methods", "00_paths.R"))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
threshold50_stage1_rda <- Sys.getenv(
  "THRESHOLD50_STAGE1_RDA",
  here("results", "dsen_results", "threshold_50", "dsen_production_results.Rda")
)
threshold50_bootstrap_rda <- Sys.getenv(
  "THRESHOLD50_BOOTSTRAP_RDA",
  here("results", "dsen_results", "threshold_50", "dsen_bootstrap_results.Rda")
)
threshold15_stage1_rda <- Sys.getenv(
  "THRESHOLD15_STAGE1_RDA",
  here("results", "dsen_results", "threshold_15", "dsen_production_results.Rda")
)
threshold15_bootstrap_rda <- Sys.getenv(
  "THRESHOLD15_BOOTSTRAP_RDA",
  here("results", "dsen_results", "threshold_15", "dsen_bootstrap_results.Rda")
)

out_threshold50_table_dir <- Sys.getenv(
  "OUT_THRESHOLD50_TABLE_DIR",
  here("results", "tables", "threshold_50")
)
out_comparison_csv <- Sys.getenv(
  "OUT_COMPARISON_CSV",
  here("results", "tables", "threshold_comparison_15_vs_50.csv")
)
out_threshold50_recovery_csv <- file.path(out_threshold50_table_dir, "target_recovery_extended.csv")
out_threshold50_recovery_summary_csv <- file.path(out_threshold50_table_dir, "target_recovery_extended_summary.csv")
out_threshold50_jaccard_csv <- file.path(out_threshold50_table_dir, "jaccard_en_dsen.csv")
out_threshold50_jaccard_per_drug_csv <- file.path(out_threshold50_table_dir, "jaccard_en_dsen_per_drug.csv")
out_threshold50_recovery_rda <- Sys.getenv(
  "OUT_THRESHOLD50_RECOVERY_RDA",
  here("results", "dsen_results", "threshold_50", "tissue_target_recovery.Rda")
)

jaccard_threshold <- as.numeric(Sys.getenv("JACCARD_THRESHOLD", "0.50"))
thresholds_env <- Sys.getenv("RECOVERY_THRESHOLDS", "0.005,0.05,0.10,0.25,0.50,0.75")
recovery_thresholds <- as.numeric(strsplit(thresholds_env, ",")[[1]])
recovery_thresholds <- recovery_thresholds[is.finite(recovery_thresholds)]
recovery_thresholds <- sort(unique(recovery_thresholds))
if (length(recovery_thresholds) == 0) {
  stop("RECOVERY_THRESHOLDS did not contain valid numeric values.")
}

dir.create(out_threshold50_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_comparison_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_threshold50_recovery_rda), recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  threshold50_stage1_rda,
  threshold50_bootstrap_rda,
  threshold15_stage1_rda,
  threshold15_bootstrap_rda
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(" - ", missing_files, collapse = "\n"))
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
load_stage1_df <- function(rda_file) {
  env <- new.env(parent = emptyenv())
  load(rda_file, envir = env)

  if (exists("dsen_results", envir = env, inherits = FALSE)) {
    df <- get("dsen_results", envir = env)
  } else if (exists("production_results", envir = env, inherits = FALSE)) {
    production_results <- get("production_results", envir = env)
    if (!"summary_df" %in% names(production_results)) {
      stop("Object 'production_results' in ", rda_file, " does not contain summary_df.")
    }
    df <- production_results$summary_df
  } else {
    stop("Could not find dsen_results or production_results in ", rda_file)
  }

  if (!"status" %in% names(df)) {
    df$status <- "success"
  }

  if (!"mse_en" %in% names(df) && "mse_standard" %in% names(df)) {
    df$mse_en <- df$mse_standard
  }
  if (!"pct_improvement" %in% names(df) && "mse_reduction_pct" %in% names(df)) {
    df$pct_improvement <- df$mse_reduction_pct
  }
  if (!"alpha_en" %in% names(df) && "alpha_standard" %in% names(df)) {
    df$alpha_en <- df$alpha_standard
  }

  if (!"pct_improvement" %in% names(df) && all(c("mse_en", "mse_dsen") %in% names(df))) {
    df$pct_improvement <- (df$mse_en - df$mse_dsen) / df$mse_en * 100
  }

  required <- c("drug", "status", "n_samples", "n_tissues", "mse_en", "mse_dsen", "pct_improvement")
  missing_cols <- required[!required %in% names(df)]
  if (length(missing_cols) > 0) {
    stop("Stage 1 data from ", rda_file, " missing columns: ", paste(missing_cols, collapse = ", "))
  }

  df
}

load_bootstrap_results <- function(rda_file) {
  env <- new.env(parent = emptyenv())
  load(rda_file, envir = env)

  if (exists("all_results", envir = env, inherits = FALSE)) {
    return(get("all_results", envir = env))
  }
  if (exists("bootstrap_results", envir = env, inherits = FALSE)) {
    return(get("bootstrap_results", envir = env))
  }

  stop("Could not find bootstrap results list in ", rda_file)
}

jaccard_index <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  union_n <- length(union(a, b))
  if (union_n == 0) return(NA_real_)
  length(intersect(a, b)) / union_n
}

format_pct <- function(x, digits = 1) {
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f%%"), x), NA_character_)
}

format_num <- function(x, digits = 1) {
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), NA_character_)
}

format_int <- function(x) {
  ifelse(is.finite(x), sprintf("%.0f", x), NA_character_)
}

value_at_threshold <- function(df, threshold, column) {
  row <- df %>% filter(abs(.data$threshold - .env$threshold) < 1e-12)
  if (nrow(row) == 0) return(NA_real_)
  as.numeric(row[[column]][1])
}

build_feature_gene_map <- function() {
  clean_env <- new.env(parent = emptyenv())
  load(file.path(paths$clean, "ccle.Rda"), envir = clean_env)
  feature_names <- colnames(clean_env$x)
  rm(clean_env)

  feature_info <- tibble(feature = feature_names) %>%
    mutate(
      type = case_when(
        grepl("_CN$", feature) ~ "CopyNumber",
        grepl("_mut\\.mis$", feature) ~ "Mutation_missense",
        grepl("_Expr$", feature) & !grepl("_\\d+_at_Expr$", feature) ~ "Expression",
        grepl("_mut\\.tr$", feature) ~ "Mutation_truncating",
        grepl("_mut\\.hs$", feature) ~ "Mutation_highsig",
        grepl("_mut\\.hsrec$", feature) ~ "Mutation_recurrent",
        grepl("_\\d+_at_Expr$", feature) ~ "Probe_expression",
        grepl("_Caution$", feature) ~ "QualityFlag",
        grepl("_ValidationFailed$|_ValidationUnavailable$", feature) ~ "QualityFlag",
        TRUE ~ "RPPA"
      ),
      gene = case_when(
        type %in% c(
          "CopyNumber", "Expression", "Mutation_missense",
          "Mutation_truncating", "Mutation_highsig", "Mutation_recurrent"
        ) ~ sub("_mut\\.hsrec$|_mut\\.hs$|_mut\\.tr$|_mut\\.mis$|_Expr$|_CN$", "", feature),
        TRUE ~ NA_character_
      )
    )

  rppa_mapping <- read_csv(here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)
  rppa_lookup <- setNames(rppa_mapping$hgnc_symbol, rppa_mapping$rppa_feature)
  rppa_idx <- feature_info$type == "RPPA" & feature_info$feature %in% names(rppa_lookup)
  feature_info$gene[rppa_idx] <- rppa_lookup[feature_info$feature[rppa_idx]]

  feature_info
}

standardize_targets <- function(feature_info) {
  genes_in_data <- unique(feature_info$gene[!is.na(feature_info$gene)])
  drug_targets <- read_csv(here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
  gene_mapping <- read_csv(here("data", "gene_name_mapping.csv"), show_col_types = FALSE)

  standardize_gene <- function(gene) {
    match_idx <- match(toupper(gene), toupper(gene_mapping$alias))
    if (!is.na(match_idx)) return(gene_mapping$hgnc_symbol[match_idx])
    gene
  }

  drug_targets %>%
    filter(!is.na(target_gene), source != "not_found") %>%
    mutate(
      target_gene = vapply(target_gene, standardize_gene, character(1)),
      target_in_features = target_gene %in% genes_in_data
    )
}

map_features_to_genes <- function(features, feature_info) {
  if (length(features) == 0) return(character(0))
  feature_info %>%
    filter(feature %in% features, !is.na(gene)) %>%
    pull(gene) %>%
    unique()
}

extract_gene_sets <- function(res, feature_info, threshold) {
  if (is.null(res$en_bootstrap$cnt) || is.null(res$dsen_bootstrap$cnt)) {
    return(NULL)
  }

  en_freq <- res$en_bootstrap$cnt / res$en_bootstrap$nboot
  en_selected <- names(en_freq)[en_freq >= threshold]
  en_genes <- map_features_to_genes(en_selected, feature_info)

  dsen_shared_freq <- res$dsen_bootstrap$cnt / res$dsen_bootstrap$nboot
  dsen_shared_selected <- names(dsen_shared_freq)[dsen_shared_freq >= threshold]
  dsen_shared_genes <- map_features_to_genes(dsen_shared_selected, feature_info)

  dsen_tissue_genes <- character(0)
  if ("tissue_cnt" %in% names(res$dsen_bootstrap)) {
    tissue_freq <- res$dsen_bootstrap$tissue_cnt / res$dsen_bootstrap$nboot
    selected_tissue <- rownames(tissue_freq)[apply(tissue_freq >= threshold, 1, any)]
    dsen_tissue_genes <- map_features_to_genes(selected_tissue, feature_info)
  }

  list(
    en_genes = en_genes,
    dsen_shared_genes = dsen_shared_genes,
    dsen_tissue_genes = dsen_tissue_genes,
    dsen_combined_genes = union(dsen_shared_genes, dsen_tissue_genes),
    n_en_features = length(en_selected),
    n_dsen_shared_features = length(dsen_shared_selected)
  )
}

compute_recovery <- function(bootstrap_results, targets_known, feature_info, thresholds) {
  recovery <- map_dfr(thresholds, function(threshold) {
    map_dfr(bootstrap_results, function(res) {
      if (is.null(res$drug) || is.null(res$status) || res$status != "success") return(NULL)

      drug_targets <- targets_known %>%
        filter(drug_name == res$drug, target_in_features) %>%
        pull(target_gene) %>%
        unique()

      if (length(drug_targets) == 0) {
        return(tibble(
          drug = res$drug,
          threshold = threshold,
          n_targets = 0,
          en_any = NA,
          dsen_shared_any = NA,
          dsen_tissue_any = NA,
          dsen_combined_any = NA
        ))
      }

      gene_sets <- extract_gene_sets(res, feature_info, threshold)
      if (is.null(gene_sets)) return(NULL)

      en_found <- intersect(drug_targets, gene_sets$en_genes)
      shared_found <- intersect(drug_targets, gene_sets$dsen_shared_genes)
      tissue_found <- intersect(drug_targets, gene_sets$dsen_tissue_genes)
      combined_found <- intersect(drug_targets, gene_sets$dsen_combined_genes)

      tibble(
        drug = res$drug,
        threshold = threshold,
        n_targets = length(drug_targets),
        en_any = length(en_found) > 0,
        dsen_shared_any = length(shared_found) > 0,
        dsen_tissue_any = length(tissue_found) > 0,
        dsen_combined_any = length(combined_found) > 0
      )
    })
  })

  valid <- recovery %>% filter(n_targets > 0)

  summary <- valid %>%
    group_by(threshold) %>%
    summarize(
      n_drugs = n_distinct(drug),
      en_any_pct = 100 * mean(en_any, na.rm = TRUE),
      en_any_n = sum(en_any, na.rm = TRUE),
      dsen_shared_any_pct = 100 * mean(dsen_shared_any, na.rm = TRUE),
      dsen_shared_any_n = sum(dsen_shared_any, na.rm = TRUE),
      dsen_tissue_any_pct = 100 * mean(dsen_tissue_any, na.rm = TRUE),
      dsen_tissue_any_n = sum(dsen_tissue_any, na.rm = TRUE),
      dsen_combined_any_pct = 100 * mean(dsen_combined_any, na.rm = TRUE),
      dsen_combined_any_n = sum(dsen_combined_any, na.rm = TRUE),
      .groups = "drop"
    )

  list(per_drug = valid, summary = summary)
}

compute_jaccard <- function(bootstrap_results, feature_info, threshold) {
  per_drug <- map_dfr(bootstrap_results, function(res) {
    if (is.null(res$drug) || is.null(res$status) || res$status != "success") return(NULL)

    gene_sets <- extract_gene_sets(res, feature_info, threshold)
    if (is.null(gene_sets)) return(NULL)

    en_genes <- gene_sets$en_genes
    dsen_shared_genes <- gene_sets$dsen_shared_genes
    dsen_tissue_genes <- gene_sets$dsen_tissue_genes
    dsen_combined_genes <- gene_sets$dsen_combined_genes

    neither_pct <- if (length(en_genes) == 0) {
      NA_real_
    } else {
      100 * length(setdiff(en_genes, dsen_combined_genes)) / length(en_genes)
    }

    tibble(
      drug = res$drug,
      threshold = threshold,
      n_en_genes = length(en_genes),
      n_dsen_shared_genes = length(dsen_shared_genes),
      n_dsen_tissue_genes = length(dsen_tissue_genes),
      n_dsen_combined_genes = length(dsen_combined_genes),
      en_genes_in_neither_dsen_block_pct = neither_pct,
      jaccard_en_vs_dsen_shared = jaccard_index(en_genes, dsen_shared_genes),
      jaccard_en_vs_dsen_combined = jaccard_index(en_genes, dsen_combined_genes),
      jaccard_dsen_shared_vs_tissue = jaccard_index(dsen_shared_genes, dsen_tissue_genes)
    )
  })

  summarize_jaccard <- function(df, column, label) {
    vals <- df[[column]]
    tibble(
      comparison = label,
      mean_jaccard = mean(vals, na.rm = TRUE),
      median_jaccard = median(vals, na.rm = TRUE),
      pct_zero_overlap = 100 * mean(vals == 0, na.rm = TRUE),
      n_drugs = sum(!is.na(vals)),
      threshold = threshold
    )
  }

  summary <- bind_rows(
    summarize_jaccard(per_drug, "jaccard_en_vs_dsen_shared", "en_vs_dsen_shared"),
    summarize_jaccard(per_drug, "jaccard_en_vs_dsen_combined", "en_vs_dsen_combined"),
    summarize_jaccard(per_drug, "jaccard_dsen_shared_vs_tissue", "dsen_shared_vs_tissue")
  )

  list(per_drug = per_drug, summary = summary)
}

collect_threshold_metrics <- function(stage1_df, recovery_summary, jaccard_result) {
  success <- stage1_df %>% filter(status == "success")
  if (nrow(success) == 0) stop("No successful rows in Stage 1 data.")

  max_idx <- which.max(success$pct_improvement)
  max_drug <- success$drug[max_idx]
  max_val <- success$pct_improvement[max_idx]

  jaccard_summary <- jaccard_result$summary
  jaccard_per_drug <- jaccard_result$per_drug

  list(
    n_drugs = nrow(success),
    n_tissues_mean = mean(success$n_tissues, na.rm = TRUE),
    n_samples_mean = mean(success$n_samples, na.rm = TRUE),
    dsen_win_rate = 100 * mean(success$pct_improvement > 0, na.rm = TRUE),
    dsen_mean_improvement = mean(success$pct_improvement, na.rm = TRUE),
    dsen_median_improvement = median(success$pct_improvement, na.rm = TRUE),
    dsen_max_improvement = max_val,
    dsen_max_improvement_drug = max_drug,
    en_target_recovery_50pct = value_at_threshold(recovery_summary, 0.50, "en_any_pct"),
    dsen_shared_recovery_50pct = value_at_threshold(recovery_summary, 0.50, "dsen_shared_any_pct"),
    dsen_tissue_recovery_50pct = value_at_threshold(recovery_summary, 0.50, "dsen_tissue_any_pct"),
    dsen_combined_recovery_50pct = value_at_threshold(recovery_summary, 0.50, "dsen_combined_any_pct"),
    en_target_recovery_0_5pct = value_at_threshold(recovery_summary, 0.005, "en_any_pct"),
    dsen_combined_recovery_0_5pct = value_at_threshold(recovery_summary, 0.005, "dsen_combined_any_pct"),
    jaccard_en_vs_dsen_shared = jaccard_summary %>%
      filter(comparison == "en_vs_dsen_shared") %>%
      pull(mean_jaccard) %>% first(),
    jaccard_en_vs_dsen_combined = jaccard_summary %>%
      filter(comparison == "en_vs_dsen_combined") %>%
      pull(mean_jaccard) %>% first(),
    jaccard_dsen_shared_vs_tissue = jaccard_summary %>%
      filter(comparison == "dsen_shared_vs_tissue") %>%
      pull(mean_jaccard) %>% first(),
    mean_en_genes = mean(jaccard_per_drug$n_en_genes, na.rm = TRUE),
    mean_dsen_shared_genes = mean(jaccard_per_drug$n_dsen_shared_genes, na.rm = TRUE),
    mean_dsen_tissue_genes = mean(jaccard_per_drug$n_dsen_tissue_genes, na.rm = TRUE),
    mean_dsen_combined_genes = mean(jaccard_per_drug$n_dsen_combined_genes, na.rm = TRUE),
    en_genes_in_neither_dsen_block = mean(jaccard_per_drug$en_genes_in_neither_dsen_block_pct, na.rm = TRUE)
  )
}

build_comparison_table <- function(metrics_15, metrics_50) {
  tibble(
    metric = c(
      "n_drugs",
      "n_tissues_mean",
      "n_samples_mean",
      "dsen_win_rate",
      "dsen_mean_improvement",
      "dsen_median_improvement",
      "dsen_max_improvement",
      "en_target_recovery_50pct",
      "dsen_shared_recovery_50pct",
      "dsen_tissue_recovery_50pct",
      "dsen_combined_recovery_50pct",
      "en_target_recovery_0.5pct",
      "dsen_combined_recovery_0.5pct",
      "jaccard_en_vs_dsen_shared",
      "jaccard_en_vs_dsen_combined",
      "jaccard_dsen_shared_vs_tissue",
      "mean_en_genes",
      "mean_dsen_shared_genes",
      "mean_dsen_tissue_genes",
      "mean_dsen_combined_genes",
      "en_genes_in_neither_dsen_block"
    ),
    threshold_15 = c(
      format_int(metrics_15$n_drugs),
      format_num(metrics_15$n_tissues_mean, 1),
      format_num(metrics_15$n_samples_mean, 0),
      format_pct(metrics_15$dsen_win_rate, 1),
      format_pct(metrics_15$dsen_mean_improvement, 2),
      format_pct(metrics_15$dsen_median_improvement, 2),
      sprintf("%.1f%% (%s)", metrics_15$dsen_max_improvement, metrics_15$dsen_max_improvement_drug),
      format_pct(metrics_15$en_target_recovery_50pct, 1),
      format_pct(metrics_15$dsen_shared_recovery_50pct, 1),
      format_pct(metrics_15$dsen_tissue_recovery_50pct, 1),
      format_pct(metrics_15$dsen_combined_recovery_50pct, 1),
      format_pct(metrics_15$en_target_recovery_0_5pct, 1),
      format_pct(metrics_15$dsen_combined_recovery_0_5pct, 1),
      format_num(metrics_15$jaccard_en_vs_dsen_shared, 3),
      format_num(metrics_15$jaccard_en_vs_dsen_combined, 3),
      format_num(metrics_15$jaccard_dsen_shared_vs_tissue, 3),
      format_num(metrics_15$mean_en_genes, 1),
      format_num(metrics_15$mean_dsen_shared_genes, 1),
      format_num(metrics_15$mean_dsen_tissue_genes, 1),
      format_num(metrics_15$mean_dsen_combined_genes, 1),
      format_pct(metrics_15$en_genes_in_neither_dsen_block, 1)
    ),
    threshold_50 = c(
      format_int(metrics_50$n_drugs),
      format_num(metrics_50$n_tissues_mean, 1),
      format_num(metrics_50$n_samples_mean, 0),
      format_pct(metrics_50$dsen_win_rate, 1),
      format_pct(metrics_50$dsen_mean_improvement, 2),
      format_pct(metrics_50$dsen_median_improvement, 2),
      sprintf("%.1f%% (%s)", metrics_50$dsen_max_improvement, metrics_50$dsen_max_improvement_drug),
      format_pct(metrics_50$en_target_recovery_50pct, 1),
      format_pct(metrics_50$dsen_shared_recovery_50pct, 1),
      format_pct(metrics_50$dsen_tissue_recovery_50pct, 1),
      format_pct(metrics_50$dsen_combined_recovery_50pct, 1),
      format_pct(metrics_50$en_target_recovery_0_5pct, 1),
      format_pct(metrics_50$dsen_combined_recovery_0_5pct, 1),
      format_num(metrics_50$jaccard_en_vs_dsen_shared, 3),
      format_num(metrics_50$jaccard_en_vs_dsen_combined, 3),
      format_num(metrics_50$jaccard_dsen_shared_vs_tissue, 3),
      format_num(metrics_50$mean_en_genes, 1),
      format_num(metrics_50$mean_dsen_shared_genes, 1),
      format_num(metrics_50$mean_dsen_tissue_genes, 1),
      format_num(metrics_50$mean_dsen_combined_genes, 1),
      format_pct(metrics_50$en_genes_in_neither_dsen_block, 1)
    )
  )
}

# ------------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------------
feature_info <- build_feature_gene_map()
targets_known <- standardize_targets(feature_info)

stage1_15 <- load_stage1_df(threshold15_stage1_rda)
stage1_50 <- load_stage1_df(threshold50_stage1_rda)
bootstrap_15 <- load_bootstrap_results(threshold15_bootstrap_rda)
bootstrap_50 <- load_bootstrap_results(threshold50_bootstrap_rda)

recovery_15 <- compute_recovery(bootstrap_15, targets_known, feature_info, recovery_thresholds)
recovery_50 <- compute_recovery(bootstrap_50, targets_known, feature_info, recovery_thresholds)
jaccard_15 <- compute_jaccard(bootstrap_15, feature_info, jaccard_threshold)
jaccard_50 <- compute_jaccard(bootstrap_50, feature_info, jaccard_threshold)

write_csv(recovery_50$per_drug, out_threshold50_recovery_csv)
write_csv(recovery_50$summary, out_threshold50_recovery_summary_csv)
write_csv(jaccard_50$summary, out_threshold50_jaccard_csv)
write_csv(jaccard_50$per_drug, out_threshold50_jaccard_per_drug_csv)

save(
  recovery_50,
  jaccard_50,
  feature_info,
  targets_known,
  file = out_threshold50_recovery_rda
)

metrics_15 <- collect_threshold_metrics(stage1_15, recovery_15$summary, jaccard_15)
metrics_50 <- collect_threshold_metrics(stage1_50, recovery_50$summary, jaccard_50)
comparison_tbl <- build_comparison_table(metrics_15, metrics_50)
write_csv(comparison_tbl, out_comparison_csv)

cat("Wrote threshold-50 tables:\n")
cat("  ", out_threshold50_recovery_csv, "\n", sep = "")
cat("  ", out_threshold50_recovery_summary_csv, "\n", sep = "")
cat("  ", out_threshold50_jaccard_csv, "\n", sep = "")
cat("  ", out_threshold50_jaccard_per_drug_csv, "\n", sep = "")
cat("  ", out_threshold50_recovery_rda, "\n", sep = "")
cat("Wrote comparison table:\n")
cat("  ", out_comparison_csv, "\n", sep = "")
