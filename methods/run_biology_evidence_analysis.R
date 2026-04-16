#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

source("methods/00_paths.R")

selection_threshold <- 100L
out_dir <- file.path(paths$results, "biology_evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading clean feature matrix and response metadata...")
data_env <- new.env()
load(file.path(paths$clean, "ccle.Rda"), envir = data_env)
load(file.path(paths$clean, "Sanger.Rdata"), envir = data_env)

x <- data_env$x
sanger <- data_env$Sanger

sample_idx <- match(rownames(x), rownames(sanger$dat))
sample_tissue <- sanger$tissue[sample_idx]
if (anyNA(sample_tissue)) {
  stop("Failed to align tissue labels to the clean feature matrix.")
}

message("Loading bootstrap results (single pass object load)...")
boot_env <- new.env()
load(file.path(paths$results, "dsen_results", "threshold_25_v2", "dsen_bootstrap_results.Rda"),
     envir = boot_env)
all_results <- boot_env$all_results

message("Loading lookup tables and comparison summaries...")
rppa_map <- read_csv("data/rppa_gene_mapping.csv", show_col_types = FALSE) %>%
  transmute(
    feature = rppa_feature,
    hgnc_symbol,
    modality = if_else(str_detect(feature, "_p"), "rppa_phospho", "rppa_total")
  )

gene_alias <- read_csv("data/gene_name_mapping.csv", show_col_types = FALSE) %>%
  mutate(alias = str_to_upper(alias), hgnc_symbol = str_to_upper(hgnc_symbol))

targets_curated <- read_csv("data/drug_targets_curated.csv", show_col_types = FALSE) %>%
  mutate(
    target_gene = str_to_upper(target_gene),
    target_gene = coalesce(
      gene_alias$hgnc_symbol[match(target_gene, gene_alias$alias)],
      target_gene
    )
  )

method_comparison <- read_csv("results/comparison/method_comparison.csv", show_col_types = FALSE)
lto_summary <- read_csv("results/comparison/leave_tissue_out/lto_summary.csv", show_col_types = FALSE) %>%
  mutate(
    lto_pct_improvement = 100 * (mean_mse_en - mean_mse_dsen_shared) / mean_mse_en
  )
lto_results <- read_csv("results/comparison/leave_tissue_out/lto_results.csv", show_col_types = FALSE) %>%
  mutate(
    lto_pct_improvement = 100 * (mse_en - mse_dsen_shared) / mse_en
  )
cosmic_gene_details <- read_csv("results/enrichment/cosmic_selected_gene_details.csv", show_col_types = FALSE) %>%
  mutate(gene_symbol = str_to_upper(gene_symbol))

feature_map <- tibble(feature = colnames(x)) %>%
  left_join(rppa_map, by = "feature") %>%
  mutate(
    modality = case_when(
      str_detect(feature, "_Expr$") ~ "expression",
      str_detect(feature, "_CN$") ~ "copy_number",
      str_detect(feature, "_mut\\.") ~ "mutation",
      !is.na(modality) ~ modality,
      TRUE ~ "unmapped"
    ),
    gene_symbol = case_when(
      modality == "expression" ~ str_remove(feature, "_Expr$"),
      modality == "copy_number" ~ str_remove(feature, "_CN$"),
      modality == "mutation" ~ str_remove(feature, "_mut\\..*$"),
      !is.na(hgnc_symbol) ~ hgnc_symbol,
      TRUE ~ NA_character_
    ),
    gene_symbol = str_to_upper(gene_symbol)
  ) %>%
  select(feature, modality, gene_symbol)

tau_score <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  v <- pmax(v, 0)
  vmax <- max(v)
  if (!is.finite(vmax) || vmax <= 0) return(NA_real_)
  sum(1 - v / vmax) / max(length(v) - 1, 1)
}

norm_entropy <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_)
  v <- pmax(v, 0)
  s <- sum(v)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  p <- v / s
  -sum(ifelse(p > 0, p * log(p), 0)) / log(length(p))
}

rank_tissue_name <- function(values, which_rank = 1L) {
  ord <- order(values, decreasing = TRUE, na.last = NA)
  if (length(ord) < which_rank) return(NA_character_)
  names(values)[ord[which_rank]]
}

rank_tissue_value <- function(values, which_rank = 1L) {
  ord <- order(values, decreasing = TRUE, na.last = NA)
  if (length(ord) < which_rank) return(NA_real_)
  unname(values[ord[which_rank]])
}

message("Computing lineage breadth metrics for all expression genes...")
expr_cols <- grep("_Expr$", colnames(x), value = TRUE)
expr_metrics <- map_dfr(expr_cols, function(col) {
  tissue_medians <- tapply(x[, col], sample_tissue, median, na.rm = TRUE)
  top_val <- rank_tissue_value(tissue_medians, 1L)
  median_val <- median(tissue_medians, na.rm = TRUE)
  tibble(
    gene_symbol = str_to_upper(str_remove(col, "_Expr$")),
    expr_tau = tau_score(tissue_medians),
    expr_entropy = norm_entropy(tissue_medians),
    expr_top_tissue = rank_tissue_name(tissue_medians, 1L),
    expr_second_tissue = rank_tissue_name(tissue_medians, 2L),
    expr_top_median = top_val,
    expr_second_median = rank_tissue_value(tissue_medians, 2L),
    expr_median_tissue = median_val,
    expr_top_to_median = ifelse(is.finite(top_val) && is.finite(median_val) && median_val != 0,
                                top_val / median_val, NA_real_)
  )
})

tau_q1 <- quantile(expr_metrics$expr_tau, probs = 0.25, na.rm = TRUE)
tau_q3 <- quantile(expr_metrics$expr_tau, probs = 0.75, na.rm = TRUE)
expr_metrics <- expr_metrics %>%
  mutate(
    expr_class = case_when(
      expr_tau <= tau_q1 ~ "broad",
      expr_tau >= tau_q3 ~ "restricted",
      TRUE ~ "intermediate"
    )
  )

extract_stable_rows <- function() {
  dsen_shared_rows <- vector("list", length(all_results))
  dsen_tissue_rows <- vector("list", length(all_results))
  en_rows <- vector("list", length(all_results))

  for (i in seq_along(all_results)) {
    drug <- names(all_results)[i]
    res <- all_results[[i]]
    if (!identical(res$status, "success")) next

    feature_names <- rownames(res$dsen_bootstrap$mean_coefs)
    shared_idx <- which(res$dsen_bootstrap$cnt > selection_threshold)
    if (length(shared_idx)) {
      dsen_shared_rows[[i]] <- tibble(
        drug = drug,
        method = "DSEN",
        block = "shared",
        relevant_tissue = NA_character_,
        feature = feature_names[shared_idx],
        selection_count = as.integer(res$dsen_bootstrap$cnt[shared_idx]),
        coef = as.numeric(res$dsen_bootstrap$mean_coefs[shared_idx, "Shared"])
      )
    }

    tissue_counts <- res$dsen_bootstrap$tissue_cnt
    tissue_names <- colnames(tissue_counts)
    tissue_long <- map_dfr(tissue_names, function(tis) {
      idx <- which(tissue_counts[, tis] > selection_threshold)
      if (!length(idx)) return(NULL)
      tibble(
        drug = drug,
        method = "DSEN",
        block = "tissue",
        relevant_tissue = tis,
        feature = feature_names[idx],
        selection_count = as.integer(tissue_counts[idx, tis]),
        coef = as.numeric(res$dsen_bootstrap$mean_coefs[idx, tis])
      )
    })
    dsen_tissue_rows[[i]] <- tissue_long

    en_idx <- which(res$en_bootstrap$cnt > selection_threshold)
    if (length(en_idx)) {
      en_rows[[i]] <- tibble(
        drug = drug,
        method = "EN",
        block = "en",
        relevant_tissue = NA_character_,
        feature = names(res$en_bootstrap$mean_coef)[en_idx],
        selection_count = as.integer(res$en_bootstrap$cnt[en_idx]),
        coef = as.numeric(res$en_bootstrap$mean_coef[en_idx])
      )
    }
  }

  list(
    dsen = bind_rows(dsen_shared_rows, dsen_tissue_rows),
    en = bind_rows(en_rows)
  )
}

message("Extracting stable DSEN and EN selections from bootstrap results...")
stable_rows <- extract_stable_rows()

dsen_feature_rows <- stable_rows$dsen %>%
  left_join(feature_map, by = "feature") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

en_feature_rows <- stable_rows$en %>%
  left_join(feature_map, by = "feature") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

aggregate_gene_block_rows <- function(df) {
  df %>%
    group_by(drug, method, block, gene_symbol) %>%
    summarise(
      relevant_tissues = paste(sort(unique(na.omit(relevant_tissue))), collapse = ";"),
      modalities = paste(sort(unique(modality)), collapse = ";"),
      n_supporting_features = n(),
      supporting_features = paste(sort(unique(feature)), collapse = ";"),
      max_selection_count = max(selection_count, na.rm = TRUE),
      max_abs_coef = max(abs(coef), na.rm = TRUE),
      signed_coef = coef[which.max(abs(coef))][1],
      .groups = "drop"
    ) %>%
    mutate(relevant_tissues = na_if(relevant_tissues, ""))
}

dsen_gene_block_rows <- aggregate_gene_block_rows(dsen_feature_rows)
en_gene_block_rows <- aggregate_gene_block_rows(en_feature_rows)

dsen_gene_assignment <- dsen_gene_block_rows %>%
  group_by(drug, gene_symbol) %>%
  summarise(
    block_assignment = case_when(
      all(sort(unique(block)) == "shared") ~ "shared",
      all(sort(unique(block)) == "tissue") ~ "tissue",
      TRUE ~ "both"
    ),
    shared_modalities = paste(sort(unique(modalities[block == "shared"])), collapse = ";"),
    tissue_modalities = paste(sort(unique(modalities[block == "tissue"])), collapse = ";"),
    tissue_blocks = paste(sort(unique(na.omit(relevant_tissues[block == "tissue"]))), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    shared_modalities = na_if(shared_modalities, ""),
    tissue_modalities = na_if(tissue_modalities, ""),
    tissue_blocks = na_if(tissue_blocks, "")
  )

dsen_gene_block_rows <- dsen_gene_block_rows %>%
  left_join(expr_metrics, by = "gene_symbol")

en_gene_block_rows <- en_gene_block_rows %>%
  left_join(expr_metrics, by = "gene_symbol")

target_hits <- targets_curated %>%
  filter(!is.na(target_gene), target_gene != "") %>%
  distinct(drug_name, target_gene, target_type, target_pathway)

dsen_gene_block_rows <- dsen_gene_block_rows %>%
  left_join(target_hits, by = c("drug" = "drug_name", "gene_symbol" = "target_gene")) %>%
  mutate(
    direct_target = !is.na(target_type),
    target_type = if_else(direct_target, target_type, NA_character_)
  )

en_gene_block_rows <- en_gene_block_rows %>%
  left_join(target_hits, by = c("drug" = "drug_name", "gene_symbol" = "target_gene")) %>%
  mutate(
    direct_target = !is.na(target_type),
    target_type = if_else(direct_target, target_type, NA_character_)
  )

dsen_gene_block_rows <- dsen_gene_block_rows %>%
  left_join(
    method_comparison %>%
      select(drug, pct_improvement_dsen),
    by = "drug"
  ) %>%
  left_join(
    lto_summary %>%
      select(drug, lto_pct_improvement, dsen_shared_win_rate),
    by = "drug"
  )

prompt_examples <- tribble(
  ~drug, ~gene_symbol, ~expected_block, ~example_group,
  "PLX-4720", "BRAF", "shared", "shared_control",
  "Dabrafenib", "BRAF", "shared", "shared_control",
  "Nutlin-3a (-)", "MDM2", "shared", "shared_control",
  "Cetuximab", "EGFR", "shared", "shared_control",
  "Gefitinib", "EGFR", "shared", "shared_control",
  "Afatinib", "EGFR", "shared", "shared_control",
  "Afatinib", "ERBB2", "shared", "shared_control",
  "Navitoclax", "BCL2", "shared", "shared_control",
  "Refametinib", "MAP2K1", "shared", "shared_control",
  "Refametinib", "MAP2K2", "shared", "shared_control",
  "CI-1040", "MAP2K1", "shared", "shared_control",
  "CI-1040", "MAP2K2", "shared", "shared_control",
  "Ruxolitinib", "JAK1", "tissue", "tissue_control",
  "Ruxolitinib", "JAK2", "tissue", "tissue_control",
  "PD173074", "FGFR1", "tissue", "tissue_control",
  "PD173074", "FGFR3", "tissue", "tissue_control",
  "Sunitinib", "KIT", "tissue", "tissue_control",
  "Sunitinib", "FLT3", "tissue", "tissue_control",
  "Rapamycin", "MTOR", "tissue", "tissue_control",
  "Tivozanib", "KDR", "tissue", "tissue_control",
  "Tivozanib", "KIT", "tissue", "tissue_control",
  "AKT inhibitor VIII", "AKT1", "tissue", "tissue_control"
) %>%
  mutate(gene_symbol = str_to_upper(gene_symbol))

prompt_status <- prompt_examples %>%
  left_join(
    dsen_gene_assignment,
    by = c("drug", "gene_symbol")
  ) %>%
  left_join(
    dsen_gene_block_rows %>%
      select(drug, gene_symbol, block, relevant_tissues, modalities,
             max_selection_count, signed_coef, expr_tau, expr_class,
             expr_top_tissue, expr_top_to_median, pct_improvement_dsen,
             lto_pct_improvement, dsen_shared_win_rate, direct_target, target_type) %>%
      mutate(block = paste0("dsen_", block)) %>%
      pivot_wider(
        names_from = block,
        values_from = c(relevant_tissues, modalities, max_selection_count,
                        signed_coef, direct_target, target_type),
        values_fn = \(x) paste(unique(na.omit(x)), collapse = ";")
      ),
    by = c("drug", "gene_symbol")
  ) %>%
  left_join(
    en_gene_block_rows %>%
      select(drug, gene_symbol, modalities, max_selection_count, signed_coef, direct_target, target_type) %>%
      rename(
        en_modalities = modalities,
        en_max_selection_count = max_selection_count,
        en_signed_coef = signed_coef,
        en_direct_target = direct_target,
        en_target_type = target_type
      ),
    by = c("drug", "gene_symbol")
  ) %>%
  mutate(
    observed_status = coalesce(block_assignment, "not_selected"),
    matches_expected = case_when(
      expected_block == "shared" & observed_status %in% c("shared", "both") ~ TRUE,
      expected_block == "tissue" & observed_status %in% c("tissue", "both") ~ TRUE,
      TRUE ~ FALSE
    )
  )

shared_metrics <- dsen_gene_block_rows %>%
  filter(block == "shared", !is.na(expr_tau))
tissue_metrics <- dsen_gene_block_rows %>%
  filter(block == "tissue", !is.na(expr_tau))
en_metrics <- en_gene_block_rows %>%
  filter(!is.na(expr_tau))

summary_metrics <- bind_rows(
  tibble(
    comparison = "DSEN shared",
    n_gene_blocks = nrow(shared_metrics),
    n_unique_gene_drug = n_distinct(shared_metrics$drug, shared_metrics$gene_symbol),
    n_unique_genes = n_distinct(shared_metrics$gene_symbol),
    median_expr_tau = median(shared_metrics$expr_tau, na.rm = TRUE),
    median_expr_top_to_median = median(shared_metrics$expr_top_to_median, na.rm = TRUE),
    broad_fraction = mean(shared_metrics$expr_class == "broad", na.rm = TRUE),
    restricted_fraction = mean(shared_metrics$expr_class == "restricted", na.rm = TRUE),
    direct_target_fraction = mean(shared_metrics$direct_target, na.rm = TRUE)
  ),
  tibble(
    comparison = "DSEN tissue",
    n_gene_blocks = nrow(tissue_metrics),
    n_unique_gene_drug = n_distinct(tissue_metrics$drug, tissue_metrics$gene_symbol),
    n_unique_genes = n_distinct(tissue_metrics$gene_symbol),
    median_expr_tau = median(tissue_metrics$expr_tau, na.rm = TRUE),
    median_expr_top_to_median = median(tissue_metrics$expr_top_to_median, na.rm = TRUE),
    broad_fraction = mean(tissue_metrics$expr_class == "broad", na.rm = TRUE),
    restricted_fraction = mean(tissue_metrics$expr_class == "restricted", na.rm = TRUE),
    direct_target_fraction = mean(tissue_metrics$direct_target, na.rm = TRUE)
  ),
  tibble(
    comparison = "EN",
    n_gene_blocks = nrow(en_metrics),
    n_unique_gene_drug = n_distinct(en_metrics$drug, en_metrics$gene_symbol),
    n_unique_genes = n_distinct(en_metrics$gene_symbol),
    median_expr_tau = median(en_metrics$expr_tau, na.rm = TRUE),
    median_expr_top_to_median = median(en_metrics$expr_top_to_median, na.rm = TRUE),
    broad_fraction = mean(en_metrics$expr_class == "broad", na.rm = TRUE),
    restricted_fraction = mean(en_metrics$expr_class == "restricted", na.rm = TRUE),
    direct_target_fraction = mean(en_metrics$direct_target, na.rm = TRUE)
  )
)

comparison_tests <- tibble(
  metric = c("expr_tau", "expr_top_to_median"),
  shared_vs_tissue_p = c(
    wilcox.test(shared_metrics$expr_tau, tissue_metrics$expr_tau)$p.value,
    wilcox.test(shared_metrics$expr_top_to_median, tissue_metrics$expr_top_to_median)$p.value
  ),
  shared_vs_en_p = c(
    wilcox.test(shared_metrics$expr_tau, en_metrics$expr_tau)$p.value,
    wilcox.test(shared_metrics$expr_top_to_median, en_metrics$expr_top_to_median)$p.value
  ),
  tissue_vs_en_p = c(
    wilcox.test(tissue_metrics$expr_tau, en_metrics$expr_tau)$p.value,
    wilcox.test(tissue_metrics$expr_top_to_median, en_metrics$expr_top_to_median)$p.value
  )
)

weak_examples <- dsen_gene_block_rows %>%
  filter(
    (block == "shared" & expr_class == "restricted") |
      (block == "tissue" & expr_class == "broad")
  ) %>%
  arrange(desc(max_selection_count), desc(max_abs_coef)) %>%
  select(
    drug, gene_symbol, block, relevant_tissues, modalities,
    expr_tau, expr_class, expr_top_tissue, expr_top_to_median,
    direct_target, pct_improvement_dsen, lto_pct_improvement
  )

target_summary <- dsen_gene_assignment %>%
  left_join(
    target_hits %>%
      rename(drug = drug_name, gene_symbol = target_gene) %>%
      distinct(drug, gene_symbol, target_type, target_pathway),
    by = c("drug", "gene_symbol")
  ) %>%
  filter(!is.na(target_type)) %>%
  count(block_assignment, target_type, sort = TRUE)

cosmic_pathway_cases <- cosmic_gene_details %>%
  filter(method == "DSEN", selection == "freq50") %>%
  count(drug, primary_pathway, sort = TRUE, name = "n_cosmic_genes")

message("Writing biology evidence outputs...")
write_csv(dsen_gene_block_rows, file.path(out_dir, "dsen_gene_block_metrics.csv"))
write_csv(dsen_gene_assignment, file.path(out_dir, "dsen_gene_assignments.csv"))
write_csv(en_gene_block_rows, file.path(out_dir, "en_gene_block_metrics.csv"))
write_csv(summary_metrics, file.path(out_dir, "summary_metrics.csv"))
write_csv(comparison_tests, file.path(out_dir, "shared_vs_tissue_tests.csv"))
write_csv(prompt_status, file.path(out_dir, "positive_control_status.csv"))
write_csv(weak_examples, file.path(out_dir, "weak_block_assignments.csv"))
write_csv(target_summary, file.path(out_dir, "target_recovery_block_summary.csv"))
write_csv(cosmic_pathway_cases, file.path(out_dir, "cosmic_pathway_case_counts.csv"))

top_shared <- dsen_gene_block_rows %>%
  filter(block == "shared") %>%
  arrange(desc(direct_target), expr_tau, desc(max_selection_count), desc(max_abs_coef)) %>%
  slice_head(n = 25)

top_tissue <- dsen_gene_block_rows %>%
  filter(block == "tissue") %>%
  arrange(desc(direct_target), desc(expr_tau), desc(max_selection_count), desc(max_abs_coef)) %>%
  slice_head(n = 25)

write_csv(top_shared, file.path(out_dir, "top_shared_candidates.csv"))
write_csv(top_tissue, file.path(out_dir, "top_tissue_candidates.csv"))

summary_lines <- c(
  sprintf("Stable selection threshold: > %d/200 bootstrap resamples", selection_threshold),
  sprintf("DSEN shared gene-block rows: %d", nrow(shared_metrics)),
  sprintf("DSEN tissue gene-block rows: %d", nrow(tissue_metrics)),
  sprintf("EN gene-block rows: %d", nrow(en_metrics)),
  sprintf("Median expression tau: shared %.3f, tissue %.3f, EN %.3f",
          median(shared_metrics$expr_tau, na.rm = TRUE),
          median(tissue_metrics$expr_tau, na.rm = TRUE),
          median(en_metrics$expr_tau, na.rm = TRUE)),
  sprintf("Median top/median tissue expression ratio: shared %.3f, tissue %.3f, EN %.3f",
          median(shared_metrics$expr_top_to_median, na.rm = TRUE),
          median(tissue_metrics$expr_top_to_median, na.rm = TRUE),
          median(en_metrics$expr_top_to_median, na.rm = TRUE)),
  sprintf("Shared vs tissue Wilcoxon p (tau): %.3g",
          comparison_tests$shared_vs_tissue_p[comparison_tests$metric == "expr_tau"]),
  sprintf("Shared vs tissue Wilcoxon p (top/median ratio): %.3g",
          comparison_tests$shared_vs_tissue_p[comparison_tests$metric == "expr_top_to_median"])
)

write_lines(summary_lines, file.path(out_dir, "summary_notes.txt"))

message("Biology evidence analysis complete.")
