#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "comparison_common.R"))

out_dir <- ensure_comparison_dir()
subset_file <- file.path(out_dir, "feature_subsets.rds")
correspondence_file <- file.path(out_dir, "rppa_expression_correspondence.csv")

message("Loading predictor matrix...")
clean_env <- new.env(parent = emptyenv())
load(file.path(paths$clean, "ccle.Rda"), envir = clean_env)
feature_names <- colnames(clean_env$x)
rm(clean_env)

message(sprintf("Total features detected: %d", length(feature_names)))

rppa_map <- parse_rppa_mapping()
rppa_features <- intersect(rppa_map$rppa_feature, feature_names)

expr_features <- feature_names[grepl("_Expr$", feature_names)]
cn_features <- feature_names[grepl("_CN$", feature_names)]
mut_features <- feature_names[grepl("_mut\\.", feature_names)]

rppa_map <- rppa_map %>%
  mutate(
    in_predictors = rppa_feature %in% feature_names,
    expression_feature_name = paste0(hgnc_symbol, "_Expr"),
    has_expression_feature = expression_feature_name %in% feature_names,
    mutation_feature_names = vapply(
      hgnc_symbol,
      function(g) {
        hits <- feature_names[startsWith(feature_names, paste0(g, "_mut."))]
        if (length(hits) == 0) "" else paste(hits, collapse = ";")
      },
      character(1)
    ),
    has_mutation_feature = mutation_feature_names != ""
  )

correspondence <- rppa_map %>%
  transmute(
    rppa_feature,
    hgnc_symbol,
    is_phospho,
    has_expression_feature,
    expression_feature_name = ifelse(has_expression_feature, expression_feature_name, ""),
    has_mutation_feature,
    mutation_feature_names
  )

write_csv(correspondence, correspondence_file)

rppa_phospho <- rppa_map %>% filter(is_phospho, in_predictors) %>% pull(rppa_feature)
rppa_total <- rppa_map %>% filter(!is_phospho, in_predictors) %>% pull(rppa_feature)
rppa_genes <- sort(unique(rppa_map$hgnc_symbol))

expr_of_rppa_genes <- intersect(paste0(rppa_genes, "_Expr"), feature_names)
cn_of_rppa_genes <- intersect(paste0(rppa_genes, "_CN"), feature_names)
mut_of_rppa_genes <- unique(unlist(lapply(
  rppa_genes,
  function(g) feature_names[startsWith(feature_names, paste0(g, "_mut."))]
), use.names = FALSE))

subsets <- list(
  all_features = feature_names,
  no_rppa = setdiff(feature_names, rppa_features),
  expr_mut_cn = union(union(expr_features, mut_features), cn_features),
  expr_only = expr_features,
  mut_only = mut_features,
  rppa_only = rppa_features,
  rppa_total_only = rppa_total,
  rppa_phospho_only = rppa_phospho,
  expr_of_rppa_genes = expr_of_rppa_genes,
  rppa_genes_expr = expr_of_rppa_genes,
  rppa_genes_all_modalities = sort(unique(c(expr_of_rppa_genes, cn_of_rppa_genes, mut_of_rppa_genes))),
  all_no_phospho = setdiff(feature_names, rppa_phospho),
  all_no_total_rppa = setdiff(feature_names, rppa_total)
)

subset_counts <- sapply(subsets, length)

payload <- list(
  generated_at = Sys.time(),
  n_total_features = length(feature_names),
  subset_counts = subset_counts,
  subsets = subsets,
  metadata = list(
    rppa_total_features = rppa_total,
    rppa_phospho_features = rppa_phospho,
    rppa_mapping = rppa_map
  )
)

saveRDS(payload, subset_file)

n_corr_yes <- sum(correspondence$has_expression_feature)
n_corr_no <- sum(!correspondence$has_expression_feature)

message("Feature subset creation complete")
for (nm in names(subset_counts)) {
  message(sprintf("  %-28s %6d", nm, subset_counts[[nm]]))
}
message(sprintf("RPPA with expression correspondence: %d", n_corr_yes))
message(sprintf("RPPA without expression correspondence: %d", n_corr_no))
message(sprintf("Saved subset file: %s", subset_file))
message(sprintf("Saved correspondence CSV: %s", correspondence_file))

append_progress(
  phase = "Phase 0.2 - Feature subsets",
  status = "completed",
  details = c(
    sprintf("Total features: %d", length(feature_names)),
    sprintf("RPPA in predictors: %d", length(rppa_features)),
    sprintf("RPPA phospho: %d", length(rppa_phospho)),
    sprintf("RPPA total: %d", length(rppa_total)),
    sprintf("RPPA with expression correspondence: %d", n_corr_yes),
    sprintf("RPPA without expression correspondence: %d", n_corr_no)
  )
)
