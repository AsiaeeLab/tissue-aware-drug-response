# =============================================================================
# Comparison Utilities for DSEN Experiment Suite
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

ensure_comparison_dir <- function() {
  out_dir <- here::here("results", "comparison")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir
}

append_progress <- function(phase, status, details = character(0), runtime = NA_character_) {
  out_dir <- ensure_comparison_dir()
  progress_file <- file.path(out_dir, "PROGRESS.md")

  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  header <- sprintf("## %s - %s", ts, phase)
  lines <- c(
    header,
    sprintf("- Status: %s", status),
    if (!is.na(runtime)) sprintf("- Runtime: %s", runtime) else NULL,
    if (length(details) > 0) sprintf("- %s", details) else "- Details: (none)",
    ""
  )

  if (!file.exists(progress_file)) {
    writeLines(c("# DSEN Comparison Progress", "", lines), progress_file)
  } else {
    cat(paste(lines, collapse = "\n"), file = progress_file, append = TRUE)
  }

  invisible(progress_file)
}

load_stage1_df_any <- function(rda_or_csv) {
  if (!file.exists(rda_or_csv)) {
    stop("File not found: ", rda_or_csv)
  }

  if (grepl("\\.csv$", rda_or_csv, ignore.case = TRUE)) {
    df <- read.csv(rda_or_csv, stringsAsFactors = FALSE)
  } else {
    env <- new.env(parent = emptyenv())
    load(rda_or_csv, envir = env)

    if (exists("production_results", envir = env, inherits = FALSE)) {
      df <- get("production_results", envir = env)$summary_df
    } else if (exists("dsen_results", envir = env, inherits = FALSE)) {
      df <- get("dsen_results", envir = env)
    } else if (exists("summary_df", envir = env, inherits = FALSE)) {
      df <- get("summary_df", envir = env)
    } else {
      stop("Could not find summary object in: ", rda_or_csv)
    }
  }

  if (!"status" %in% names(df)) df$status <- "success"
  if (!"mse_en" %in% names(df) && "mse_standard" %in% names(df)) df$mse_en <- df$mse_standard
  if (!"pct_improvement" %in% names(df) && all(c("mse_en", "mse_dsen") %in% names(df))) {
    df$pct_improvement <- (df$mse_en - df$mse_dsen) / df$mse_en * 100
  }
  df
}

load_bootstrap_results_any <- function(file_path) {
  if (!file.exists(file_path)) return(NULL)

  env <- new.env(parent = emptyenv())
  load(file_path, envir = env)

  if (exists("all_results", envir = env, inherits = FALSE)) return(get("all_results", envir = env))
  if (exists("bootstrap_results", envir = env, inherits = FALSE)) return(get("bootstrap_results", envir = env))

  if (exists("production_results", envir = env, inherits = FALSE)) {
    pr <- get("production_results", envir = env)
    if (all(c("summary_df", "en_bootstrap", "dsen_bootstrap") %in% names(pr))) {
      summary_df <- pr$summary_df
      drugs <- summary_df$drug
      out <- lapply(drugs, function(dr) {
        idx <- which(names(pr$en_bootstrap) == dr)
        if (length(idx) == 0) idx <- which(drugs == dr)
        i <- idx[1]
        list(
          drug = dr,
          status = summary_df$status[summary_df$drug == dr][1],
          en_bootstrap = pr$en_bootstrap[[i]],
          dsen_bootstrap = pr$dsen_bootstrap[[i]]
        )
      })
      names(out) <- drugs
      return(out)
    }
  }

  NULL
}

build_feature_gene_map <- function(paths) {
  clean_env <- new.env(parent = emptyenv())
  load(file.path(paths$clean, "ccle.Rda"), envir = clean_env)
  feature_names <- colnames(clean_env$x)
  rm(clean_env)

  feature_info <- tibble(feature = feature_names) %>%
    mutate(
      type = case_when(
        grepl("_CN$", feature) ~ "CopyNumber",
        grepl("_mut\\.", feature) ~ "Mutation",
        grepl("_Expr$", feature) ~ "Expression",
        TRUE ~ "Other"
      ),
      gene = case_when(
        type %in% c("CopyNumber", "Mutation", "Expression") ~
          sub("_mut\\.hsrec$|_mut\\.hs$|_mut\\.tr$|_mut\\.mis$|_Expr$|_CN$", "", feature),
        TRUE ~ NA_character_
      )
    )

  rppa_mapping <- read_csv(here::here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)
  rppa_lookup <- setNames(rppa_mapping$hgnc_symbol, rppa_mapping$rppa_feature)

  is_rppa <- feature_info$feature %in% names(rppa_lookup)
  feature_info$type[is_rppa] <- "RPPA"
  feature_info$gene[is_rppa] <- rppa_lookup[feature_info$feature[is_rppa]]

  feature_info
}

standardize_targets <- function(feature_info) {
  genes_in_data <- unique(feature_info$gene[!is.na(feature_info$gene)])
  drug_targets <- read_csv(here::here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
  gene_mapping <- read_csv(here::here("data", "gene_name_mapping.csv"), show_col_types = FALSE)

  standardize_gene <- function(gene) {
    idx <- match(toupper(gene), toupper(gene_mapping$alias))
    if (!is.na(idx)) return(gene_mapping$hgnc_symbol[idx])
    gene
  }

  drug_targets %>%
    filter(!is.na(target_gene), source != "not_found") %>%
    mutate(
      target_gene = vapply(target_gene, standardize_gene, character(1)),
      target_in_features = target_gene %in% genes_in_data
    )
}

parse_rppa_mapping <- function() {
  map <- read_csv(here::here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)
  # Match the production paper accounting: phospho features are those explicitly
  # encoded with "_pS/_pT/_pY" in RPPA feature names (39 features).
  map$is_phospho <- grepl("_p[STY][0-9]+", map$rppa_feature)
  map
}

safe_pct <- function(num, den) {
  ifelse(is.finite(den) & den != 0, 100 * num / den, NA_real_)
}
