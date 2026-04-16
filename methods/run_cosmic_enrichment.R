#!/usr/bin/env Rscript
# Purpose: Run COSMIC/CGC enrichment validation for DSEN bootstrap outputs.
# Outputs:
#   results/enrichment/feature_gene_mapping.csv
#   results/enrichment/cosmic_enrichment_per_drug.csv
#   results/enrichment/cosmic_enrichment_summary.csv
#   results/enrichment/cosmic_directional_enrichment.csv
#   results/enrichment/cosmic_tissue_enrichment.csv (only when Phase 2 is positive)
#   results/enrichment/figures/*.pdf

suppressPackageStartupMessages({
  library(here)
})

source(here("methods", "00_paths.R"))

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

BOOTSTRAP_FILE <- here(
  "results", "dsen_results", "threshold_25_v2", "dsen_bootstrap_results.Rda"
)
RPPA_MAP_FILE <- here("data", "rppa_gene_mapping.csv")
TARGET_FILE <- here("data", "drug_targets_curated.csv")
COSMIC_FILE <- here("data", "cosmic_cgc.csv")

OUT_DIR <- here("results", "enrichment")
FIG_DIR <- file.path(OUT_DIR, "figures")
CHECKPOINT_FILE <- file.path(OUT_DIR, "cosmic_enrichment_checkpoint.rds")

FEATURE_MAP_FILE <- file.path(OUT_DIR, "feature_gene_mapping.csv")
PER_DRUG_FILE <- file.path(OUT_DIR, "cosmic_enrichment_per_drug.csv")
SUMMARY_FILE <- file.path(OUT_DIR, "cosmic_enrichment_summary.csv")
DIRECTIONAL_FILE <- file.path(OUT_DIR, "cosmic_directional_enrichment.csv")
TISSUE_FILE <- file.path(OUT_DIR, "cosmic_tissue_enrichment.csv")
SELECTED_DETAIL_FILE <- file.path(OUT_DIR, "cosmic_selected_gene_details.csv")

FREQ_THRESHOLDS <- c(freq25 = 50L, freq50 = 100L, freq75 = 150L)
TOP_K <- 50L
CHECKPOINT_EVERY <- as.integer(Sys.getenv("CHECKPOINT_EVERY", "10"))
if (is.na(CHECKPOINT_EVERY) || CHECKPOINT_EVERY < 1L) {
  CHECKPOINT_EVERY <- 10L
}
RESUME <- Sys.getenv("RESUME", "TRUE") == "TRUE"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, paste(..., collapse = "")))
}

trim_char <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}

normalize_gene <- function(x) {
  x <- trim_char(x)
  x <- toupper(x)
  x
}

collapse_tokens <- function(x, sort_tokens = TRUE) {
  x <- trim_char(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("")
  ux <- unique(x)
  if (sort_tokens) ux <- sort(ux)
  paste(ux, collapse = ";")
}

role_flags_from_text <- function(txt) {
  txt <- tolower(paste(txt, collapse = ";"))
  is_og <- grepl("oncogene|\\bog\\b", txt)
  is_tsg <- grepl("tsg|tumou?r\\s*suppressor", txt)
  is_fusion <- grepl("fusion", txt)
  role <- c(
    if (is_og) "OG",
    if (is_tsg) "TSG",
    if (is_fusion) "fusion"
  )
  list(
    role = if (length(role) == 0L) NA_character_ else paste(role, collapse = ";"),
    is_og = is_og,
    is_tsg = is_tsg,
    is_fusion = is_fusion
  )
}

coerce_logical <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  x_chr <- trimws(tolower(as.character(x)))
  out <- x_chr %in% c("1", "true", "t", "yes", "y", "og", "tsg", "fusion")
  out[is.na(x_chr)] <- FALSE
  out
}

align_vector <- function(vec, feature_names, label) {
  if (is.null(names(vec))) {
    if (length(vec) != length(feature_names)) {
      stop(
        sprintf(
          "%s length mismatch: got %d, expected %d",
          label, length(vec), length(feature_names)
        )
      )
    }
    return(as.numeric(vec))
  }
  idx <- match(feature_names, names(vec))
  if (anyNA(idx)) {
    stop(
      sprintf(
        "%s missing %d features after name alignment",
        label, sum(is.na(idx))
      )
    )
  }
  as.numeric(vec[idx])
}

align_matrix_rows <- function(mat, feature_names, label) {
  if (is.null(rownames(mat))) {
    if (nrow(mat) != length(feature_names)) {
      stop(
        sprintf(
          "%s row mismatch: got %d, expected %d",
          label, nrow(mat), length(feature_names)
        )
      )
    }
    return(mat)
  }
  idx <- match(feature_names, rownames(mat))
  if (anyNA(idx)) {
    stop(
      sprintf(
        "%s missing %d features after rowname alignment",
        label, sum(is.na(idx))
      )
    )
  }
  mat[idx, , drop = FALSE]
}

safe_fisher <- function(selected_genes, cosmic_genes, universe_genes) {
  selected <- intersect(unique(selected_genes), universe_genes)
  cosmic <- intersect(unique(cosmic_genes), universe_genes)

  a <- length(intersect(selected, cosmic))
  b <- length(setdiff(selected, cosmic))
  c <- length(setdiff(cosmic, selected))
  d <- length(universe_genes) - a - b - c
  if (d < 0L) d <- 0L

  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  p_val <- 1
  or <- NA_real_

  if (sum(tab[1, ]) > 0L && sum(tab[, 1]) > 0L) {
    ft <- suppressWarnings(fisher.test(tab, alternative = "greater"))
    p_val <- unname(ft$p.value)
    or <- unname(ft$estimate)
  }

  or_ha <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

  list(
    a = a,
    b = b,
    c = c,
    d = d,
    p_value = p_val,
    odds_ratio = or,
    odds_ratio_ha = or_ha
  )
}

rows_to_df <- function(rows) {
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

download_ncg_proxy <- function(dest_csv) {
  url <- "http://network-cancer-genes.org/download.php"
  tmp_tsv <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp_tsv), add = TRUE)

  cmd <- sprintf(
    "curl -s -L -X POST -d %s %s -o %s",
    shQuote("downloadcancergenes=Download"),
    shQuote(url),
    shQuote(tmp_tsv)
  )
  status <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (status != 0L || !file.exists(tmp_tsv) || file.info(tmp_tsv)$size <= 0L) {
    stop("Failed to download NCG cancer driver annotation file.")
  }

  ncg <- read.delim(
    tmp_tsv,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    quote = "\""
  )
  req <- c("symbol", "cgc_annotation")
  if (!all(req %in% colnames(ncg))) {
    stop("Downloaded NCG file does not contain required columns.")
  }

  ncg$symbol <- normalize_gene(ncg$symbol)
  ncg$cgc_annotation <- trim_char(ncg$cgc_annotation)
  ncg <- ncg[!is.na(ncg$symbol) & !is.na(ncg$cgc_annotation), , drop = FALSE]
  if (nrow(ncg) == 0L) {
    stop("Downloaded NCG file has no rows with CGC annotation.")
  }

  split_gene <- split(ncg, ncg$symbol)
  cosmic_rows <- lapply(names(split_gene), function(gene) {
    gdf <- split_gene[[gene]]
    flags <- role_flags_from_text(gdf$cgc_annotation)
    tumour_sites <- trim_char(gdf$primary_site)
    tumour_sites <- tumour_sites[tolower(tumour_sites) != "multiple"]
    if (length(tumour_sites) == 0L) {
      tumour_sites <- trim_char(gdf$primary_site)
    }
    data.frame(
      gene_symbol = gene,
      role = flags$role,
      is_og = as.integer(flags$is_og),
      is_tsg = as.integer(flags$is_tsg),
      is_fusion = as.integer(flags$is_fusion),
      tumour_types_somatic = collapse_tokens(tumour_sites),
      cancer_types = collapse_tokens(gdf$cancer_type),
      source_dataset = "NCG7.2_proxy",
      source_url = url,
      stringsAsFactors = FALSE
    )
  })
  cosmic <- do.call(rbind, cosmic_rows)
  cosmic <- cosmic[order(cosmic$gene_symbol), , drop = FALSE]
  write.csv(cosmic, dest_csv, row.names = FALSE)
  cosmic
}

canonicalize_cosmic <- function(cosmic_raw) {
  cn <- colnames(cosmic_raw)
  cn_norm <- tolower(gsub("[^a-z0-9]+", "_", cn))
  colnames(cosmic_raw) <- cn_norm

  gene_col <- intersect(
    c("gene_symbol", "gene", "symbol", "hgnc_symbol", "gene_name"),
    colnames(cosmic_raw)
  )
  if (length(gene_col) == 0L) stop("No gene symbol column found in cosmic file.")
  cosmic_raw$gene_symbol <- normalize_gene(cosmic_raw[[gene_col[1]]])

  role_col <- intersect(
    c("role", "role_in_cancer", "cgc_annotation", "cancer_role"),
    colnames(cosmic_raw)
  )
  if (length(role_col) > 0L) {
    cosmic_raw$role <- trim_char(cosmic_raw[[role_col[1]]])
  } else {
    cosmic_raw$role <- NA_character_
  }

  tumour_col <- intersect(
    c(
      "tumour_types_somatic", "tumor_types_somatic",
      "primary_site", "tumour_type", "tumor_type", "cancer_type"
    ),
    colnames(cosmic_raw)
  )
  if (length(tumour_col) > 0L) {
    cosmic_raw$tumour_types_somatic <- trim_char(cosmic_raw[[tumour_col[1]]])
  } else {
    cosmic_raw$tumour_types_somatic <- NA_character_
  }

  if ("is_og" %in% colnames(cosmic_raw)) {
    cosmic_raw$is_og <- coerce_logical(cosmic_raw$is_og)
  } else {
    cosmic_raw$is_og <- FALSE
  }
  if ("is_tsg" %in% colnames(cosmic_raw)) {
    cosmic_raw$is_tsg <- coerce_logical(cosmic_raw$is_tsg)
  } else {
    cosmic_raw$is_tsg <- FALSE
  }
  if ("is_fusion" %in% colnames(cosmic_raw)) {
    cosmic_raw$is_fusion <- coerce_logical(cosmic_raw$is_fusion)
  } else {
    cosmic_raw$is_fusion <- FALSE
  }

  # Fill role-derived flags when role text exists.
  has_role <- !is.na(cosmic_raw$role)
  if (any(has_role)) {
    role_flags <- lapply(cosmic_raw$role[has_role], role_flags_from_text)
    cosmic_raw$is_og[has_role] <- cosmic_raw$is_og[has_role] |
      vapply(role_flags, function(x) x$is_og, logical(1))
    cosmic_raw$is_tsg[has_role] <- cosmic_raw$is_tsg[has_role] |
      vapply(role_flags, function(x) x$is_tsg, logical(1))
    cosmic_raw$is_fusion[has_role] <- cosmic_raw$is_fusion[has_role] |
      vapply(role_flags, function(x) x$is_fusion, logical(1))
  }

  # Collapse duplicate genes.
  cosmic_raw <- cosmic_raw[!is.na(cosmic_raw$gene_symbol), , drop = FALSE]
  split_gene <- split(cosmic_raw, cosmic_raw$gene_symbol)
  collapsed <- lapply(names(split_gene), function(g) {
    gdf <- split_gene[[g]]
    role_txt <- role_flags_from_text(c(gdf$role, ifelse(gdf$is_og, "OG", NA), ifelse(gdf$is_tsg, "TSG", NA), ifelse(gdf$is_fusion, "fusion", NA)))
    data.frame(
      gene_symbol = g,
      role = role_txt$role,
      is_og = as.integer(any(gdf$is_og)),
      is_tsg = as.integer(any(gdf$is_tsg)),
      is_fusion = as.integer(any(gdf$is_fusion)),
      tumour_types_somatic = collapse_tokens(gdf$tumour_types_somatic),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, collapsed)
  out <- out[order(out$gene_symbol), , drop = FALSE]
  out
}

prepare_cosmic <- function(cosmic_file) {
  if (!file.exists(cosmic_file)) {
    log_msg("No local data/cosmic_cgc.csv found. Downloading open proxy from NCG...")
    cosmic <- download_ncg_proxy(cosmic_file)
    return(canonicalize_cosmic(cosmic))
  }
  log_msg("Using existing cosmic file: ", cosmic_file)
  raw <- read.csv(cosmic_file, stringsAsFactors = FALSE, check.names = FALSE)
  canonicalize_cosmic(raw)
}

build_feature_mapping <- function(feature_names, rppa_map) {
  rppa_map$rppa_feature <- as.character(rppa_map$rppa_feature)
  rppa_map$hgnc_symbol <- normalize_gene(rppa_map$hgnc_symbol)

  is_mut <- grepl("_mut\\..*$", feature_names)
  is_cn <- grepl("_CN$", feature_names)
  is_rppa <- feature_names %in% rppa_map$rppa_feature

  feature_type <- ifelse(
    is_mut, "mutation",
    ifelse(is_cn, "copy_number", ifelse(is_rppa, "rppa", "expression"))
  )

  gene_symbol <- feature_names
  gene_symbol[is_mut] <- sub("_mut\\..*$", "", feature_names[is_mut])
  gene_symbol[is_cn] <- sub("_CN$", "", feature_names[is_cn])
  if (any(is_rppa)) {
    idx <- match(feature_names[is_rppa], rppa_map$rppa_feature)
    gene_symbol[is_rppa] <- rppa_map$hgnc_symbol[idx]
  }
  gene_symbol <- normalize_gene(gene_symbol)

  data.frame(
    feature_name = feature_names,
    gene_symbol = gene_symbol,
    feature_type = feature_type,
    stringsAsFactors = FALSE
  )
}

get_primary_pathway_map <- function(target_file) {
  targets <- read.csv(target_file, stringsAsFactors = FALSE)
  targets$drug_name <- as.character(targets$drug_name)
  targets$target_pathway <- as.character(targets$target_pathway)
  targets$target_type <- as.character(targets$target_type)

  primary <- targets[tolower(targets$target_type) == "primary", , drop = FALSE]
  p_map <- tapply(
    primary$target_pathway,
    primary$drug_name,
    function(x) sort(unique(trim_char(x)))[1]
  )

  all_drugs <- unique(targets$drug_name)
  missing <- setdiff(all_drugs, names(p_map))
  if (length(missing) > 0L) {
    for (drug in missing) {
      tmp <- targets$target_pathway[targets$drug_name == drug]
      tmp <- trim_char(tmp)
      tmp <- tmp[!is.na(tmp)]
      p_map[drug] <- if (length(tmp) == 0L) "Unknown" else names(sort(table(tmp), decreasing = TRUE))[1]
    }
  }
  p_map
}

# ------------------------------------------------------------------------------
# Load input data and create mapping
# ------------------------------------------------------------------------------

log_msg(strrep("=", 80))
log_msg("Starting COSMIC enrichment analysis pipeline")
log_msg(strrep("=", 80))

if (!file.exists(BOOTSTRAP_FILE)) {
  stop("Bootstrap results file not found: ", BOOTSTRAP_FILE)
}
if (!file.exists(RPPA_MAP_FILE)) {
  stop("RPPA mapping file not found: ", RPPA_MAP_FILE)
}
if (!file.exists(TARGET_FILE)) {
  stop("Drug target mapping file not found: ", TARGET_FILE)
}

log_msg("Loading bootstrap results...")
load(BOOTSTRAP_FILE) # creates all_results
if (!exists("all_results")) {
  stop("Could not find `all_results` in bootstrap file.")
}
if (length(all_results) == 0L) {
  stop("`all_results` is empty.")
}
drug_names <- names(all_results)
if (is.null(drug_names) || any(!nzchar(drug_names))) {
  drug_names <- vapply(all_results, function(x) as.character(x$drug), character(1))
  names(all_results) <- drug_names
}
log_msg("Loaded bootstrap list for ", length(all_results), " drugs.")

rppa_map <- read.csv(RPPA_MAP_FILE, stringsAsFactors = FALSE)

example_drug <- drug_names[1]
feature_names <- rownames(all_results[[example_drug]]$dsen_bootstrap$mean_coefs)
if (is.null(feature_names)) {
  stop("Feature names are missing from dsen_bootstrap$mean_coefs.")
}
log_msg("Feature matrix size: ", length(feature_names), " features.")

feature_map <- build_feature_mapping(feature_names, rppa_map)
write.csv(feature_map, FEATURE_MAP_FILE, row.names = FALSE)
log_msg("Saved feature mapping: ", FEATURE_MAP_FILE)
log_msg(
  "Feature types: ",
  paste(names(table(feature_map$feature_type)), table(feature_map$feature_type), collapse = ", ")
)
log_msg("Unique mapped genes: ", length(unique(feature_map$gene_symbol)))

cosmic_df <- prepare_cosmic(COSMIC_FILE)
write.csv(cosmic_df, COSMIC_FILE, row.names = FALSE)
log_msg("COSMIC/CGC gene list ready: ", COSMIC_FILE, " (", nrow(cosmic_df), " genes)")

primary_pathway <- get_primary_pathway_map(TARGET_FILE)

valid_idx <- which(!is.na(feature_map$gene_symbol))
gene_vec_valid <- feature_map$gene_symbol[valid_idx]
gene_factor <- factor(gene_vec_valid)
gene_levels <- levels(gene_factor)
gene_split <- split(seq_along(gene_vec_valid), gene_factor)
universe_genes <- gene_levels

cosmic_df <- cosmic_df[cosmic_df$gene_symbol %in% universe_genes, , drop = FALSE]
cosmic_set <- cosmic_df$gene_symbol
og_set <- cosmic_df$gene_symbol[cosmic_df$is_og == 1L]
tsg_set <- cosmic_df$gene_symbol[cosmic_df$is_tsg == 1L]
fusion_set <- cosmic_df$gene_symbol[cosmic_df$is_fusion == 1L]

log_msg("COSMIC genes in feature universe: ", length(cosmic_set))

cosmic_lookup <- cosmic_df
rownames(cosmic_lookup) <- cosmic_lookup$gene_symbol

summarize_gene_metrics <- function(coef_vec, count_vec) {
  coef_valid <- coef_vec[valid_idx]
  count_valid <- count_vec[valid_idx]

  max_count <- tapply(count_valid, gene_factor, max)
  max_abs_coef <- tapply(abs(coef_valid), gene_factor, max)

  signed_coef <- numeric(length(gene_split))
  names(signed_coef) <- names(gene_split)
  for (i in seq_along(gene_split)) {
    idx_local <- gene_split[[i]]
    vals <- coef_valid[idx_local]
    if (all(is.na(vals))) {
      signed_coef[i] <- NA_real_
    } else {
      signed_coef[i] <- vals[which.max(abs(vals))]
    }
  }

  list(
    max_count = max_count,
    max_abs_coef = max_abs_coef,
    signed_coef = signed_coef
  )
}

# ------------------------------------------------------------------------------
# Phase 2: Per-drug enrichment + directional analysis with checkpointing
# ------------------------------------------------------------------------------

per_rows <- list()
direction_rows <- list()
gene_rows <- list()
completed_drugs <- character(0)

if (RESUME && file.exists(CHECKPOINT_FILE)) {
  cp <- readRDS(CHECKPOINT_FILE)
  per_rows <- cp$per_rows
  direction_rows <- cp$direction_rows
  gene_rows <- cp$gene_rows
  completed_drugs <- cp$completed_drugs
  log_msg(
    "Resuming from checkpoint: ",
    length(completed_drugs), "/", length(drug_names), " drugs already processed."
  )
}

remaining_drugs <- setdiff(drug_names, completed_drugs)
log_msg("Remaining drugs to process: ", length(remaining_drugs))

write_intermediate <- function(force = FALSE) {
  n_done <- length(completed_drugs)
  if (!force && (n_done %% CHECKPOINT_EVERY != 0L) && n_done != length(drug_names)) {
    return(invisible(NULL))
  }

  per_df <- rows_to_df(per_rows)
  if (nrow(per_df) > 0L) write.csv(per_df, PER_DRUG_FILE, row.names = FALSE)

  dir_df <- rows_to_df(direction_rows)
  if (nrow(dir_df) > 0L) write.csv(dir_df, DIRECTIONAL_FILE, row.names = FALSE)

  gene_df <- rows_to_df(gene_rows)
  if (nrow(gene_df) > 0L) write.csv(gene_df, SELECTED_DETAIL_FILE, row.names = FALSE)

  saveRDS(
    list(
      per_rows = per_rows,
      direction_rows = direction_rows,
      gene_rows = gene_rows,
      completed_drugs = completed_drugs
    ),
    CHECKPOINT_FILE
  )
  log_msg("Checkpoint saved (", n_done, "/", length(drug_names), ").")
}

for (ii in seq_along(remaining_drugs)) {
  drug <- remaining_drugs[ii]
  res <- all_results[[drug]]

  if (is.null(res$dsen_bootstrap) || is.null(res$en_bootstrap)) {
    log_msg("Skipping ", drug, " because bootstrap fields are missing.")
    completed_drugs <- c(completed_drugs, drug)
    next
  }

  # DSEN shared vectors
  dsen_coefs_raw <- res$dsen_bootstrap$mean_coefs[, "Shared"]
  dsen_cnt_raw <- res$dsen_bootstrap$cnt
  dsen_coefs <- align_vector(dsen_coefs_raw, feature_names, paste0(drug, ":DSEN shared coefs"))
  dsen_cnt <- align_vector(dsen_cnt_raw, feature_names, paste0(drug, ":DSEN counts"))

  # EN vectors
  en_coefs_raw <- res$en_bootstrap$mean_coef
  en_cnt_raw <- res$en_bootstrap$cnt
  en_coefs <- align_vector(en_coefs_raw, feature_names, paste0(drug, ":EN coefs"))
  en_cnt <- align_vector(en_cnt_raw, feature_names, paste0(drug, ":EN counts"))

  dsen_gene <- summarize_gene_metrics(dsen_coefs, dsen_cnt)
  en_gene <- summarize_gene_metrics(en_coefs, en_cnt)

  method_payload <- list(
    DSEN = dsen_gene,
    EN = en_gene
  )

  for (method in names(method_payload)) {
    gsum <- method_payload[[method]]

    # Frequency-based selection settings
    for (sel in names(FREQ_THRESHOLDS)) {
      thr <- FREQ_THRESHOLDS[[sel]]
      sel_genes <- names(gsum$max_count)[gsum$max_count > thr]
      fish <- safe_fisher(sel_genes, cosmic_set, universe_genes)
      cosmic_found <- intersect(sel_genes, cosmic_set)

      per_rows[[length(per_rows) + 1L]] <- data.frame(
        drug = drug,
        method = method,
        selection = sel,
        n_selected_genes = length(sel_genes),
        n_cosmic_in_selected = fish$a,
        n_cosmic_in_universe = length(cosmic_set),
        odds_ratio = fish$odds_ratio,
        odds_ratio_ha = fish$odds_ratio_ha,
        p_value = fish$p_value,
        cosmic_genes_found = collapse_tokens(cosmic_found),
        stringsAsFactors = FALSE
      )

      if (sel == "freq50") {
        signed <- gsum$signed_coef
        signed_sel <- signed[sel_genes]
        resistance <- names(signed_sel)[!is.na(signed_sel) & signed_sel > 0]
        sensitivity <- names(signed_sel)[!is.na(signed_sel) & signed_sel < 0]

        og_res <- safe_fisher(resistance, og_set, universe_genes)
        tsg_sens <- safe_fisher(sensitivity, tsg_set, universe_genes)
        og_sens <- safe_fisher(sensitivity, og_set, universe_genes)
        tsg_res <- safe_fisher(resistance, tsg_set, universe_genes)

        direction_rows[[length(direction_rows) + 1L]] <- data.frame(
          drug = drug,
          method = method,
          selection = sel,
          n_resistance_genes = length(resistance),
          n_sensitivity_genes = length(sensitivity),
          n_og_in_resistance = og_res$a,
          n_tsg_in_sensitivity = tsg_sens$a,
          n_og_in_sensitivity = og_sens$a,
          n_tsg_in_resistance = tsg_res$a,
          odds_ratio_og_resistance = og_res$odds_ratio,
          odds_ratio_ha_og_resistance = og_res$odds_ratio_ha,
          p_og_resistance = og_res$p_value,
          odds_ratio_tsg_sensitivity = tsg_sens$odds_ratio,
          odds_ratio_ha_tsg_sensitivity = tsg_sens$odds_ratio_ha,
          p_tsg_sensitivity = tsg_sens$p_value,
          resistance_genes_found = collapse_tokens(intersect(resistance, cosmic_set)),
          sensitivity_genes_found = collapse_tokens(intersect(sensitivity, cosmic_set)),
          stringsAsFactors = FALSE
        )

        # Save per-gene directional details for summaries/figures.
        cosmic_selected <- intersect(sel_genes, cosmic_set)
        if (length(cosmic_selected) > 0L) {
          for (g in cosmic_selected) {
            cmeta <- cosmic_lookup[g, , drop = FALSE]
            gcoef <- signed[g]
            gsign <- ifelse(
              is.na(gcoef), "zero",
              ifelse(gcoef > 0, "positive", ifelse(gcoef < 0, "negative", "zero"))
            )
            gene_rows[[length(gene_rows) + 1L]] <- data.frame(
              drug = drug,
              method = method,
              selection = sel,
              gene_symbol = g,
              coef_value = gcoef,
              coef_sign = gsign,
              is_og = as.integer(cmeta$is_og),
              is_tsg = as.integer(cmeta$is_tsg),
              is_fusion = as.integer(cmeta$is_fusion),
              role = as.character(cmeta$role),
              primary_pathway = as.character(primary_pathway[drug]),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    # Top-K by gene-level maximum absolute coefficient
    ranked <- sort(gsum$max_abs_coef, decreasing = TRUE)
    ranked <- ranked[!is.na(ranked)]
    top_genes <- names(ranked)[seq_len(min(TOP_K, length(ranked)))]
    fish_top <- safe_fisher(top_genes, cosmic_set, universe_genes)
    top_found <- intersect(top_genes, cosmic_set)

    per_rows[[length(per_rows) + 1L]] <- data.frame(
      drug = drug,
      method = method,
      selection = "top50",
      n_selected_genes = length(top_genes),
      n_cosmic_in_selected = fish_top$a,
      n_cosmic_in_universe = length(cosmic_set),
      odds_ratio = fish_top$odds_ratio,
      odds_ratio_ha = fish_top$odds_ratio_ha,
      p_value = fish_top$p_value,
      cosmic_genes_found = collapse_tokens(top_found),
      stringsAsFactors = FALSE
    )
  }

  completed_drugs <- c(completed_drugs, drug)

  if (ii %% 5L == 0L || ii == length(remaining_drugs)) {
    log_msg(
      "Processed ", length(completed_drugs), "/", length(drug_names),
      " drugs (latest: ", drug, ")."
    )
  }
  write_intermediate(force = FALSE)
}

write_intermediate(force = TRUE)

per_drug <- rows_to_df(per_rows)
if (nrow(per_drug) == 0L) {
  stop("No per-drug enrichment rows were produced.")
}
per_drug$p_adj_bh <- ave(
  per_drug$p_value,
  interaction(per_drug$method, per_drug$selection, drop = TRUE),
  FUN = function(x) p.adjust(x, method = "BH")
)
write.csv(per_drug, PER_DRUG_FILE, row.names = FALSE)
log_msg("Wrote per-drug enrichment results: ", PER_DRUG_FILE)

directional <- rows_to_df(direction_rows)
if (nrow(directional) > 0L) {
  directional$p_og_resistance_bh <- NA_real_
  directional$p_tsg_sensitivity_bh <- NA_real_
  for (m in unique(directional$method)) {
    idx <- directional$method == m
    directional$p_og_resistance_bh[idx] <- p.adjust(
      directional$p_og_resistance[idx], method = "BH"
    )
    directional$p_tsg_sensitivity_bh[idx] <- p.adjust(
      directional$p_tsg_sensitivity[idx], method = "BH"
    )
  }
  write.csv(directional, DIRECTIONAL_FILE, row.names = FALSE)
  log_msg("Wrote directional enrichment results: ", DIRECTIONAL_FILE)
}

gene_detail <- rows_to_df(gene_rows)
if (nrow(gene_detail) > 0L) {
  write.csv(gene_detail, SELECTED_DETAIL_FILE, row.names = FALSE)
}

# ------------------------------------------------------------------------------
# Phase 2.2 Aggregate summary
# ------------------------------------------------------------------------------

summary_rows <- list()
add_summary <- function(...) {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)
}

sel_levels <- unique(per_drug$selection)
sel_levels <- sel_levels[order(match(sel_levels, c("freq50", "freq25", "freq75", "top50")))]

phase2_freq50_p <- NA_real_
phase2_freq50_med_dsen <- NA_real_
phase2_freq50_med_en <- NA_real_

for (sel in sel_levels) {
  tmp <- per_drug[per_drug$selection == sel, , drop = FALSE]

  for (m in c("DSEN", "EN")) {
    sub <- tmp[tmp$method == m, , drop = FALSE]
    if (nrow(sub) == 0L) next
    add_summary(
      analysis = "significance_fraction",
      selection = sel,
      method = m,
      group_name = "all_drugs",
      n = nrow(sub),
      n_sig_raw = sum(sub$p_value < 0.05, na.rm = TRUE),
      frac_sig_raw = mean(sub$p_value < 0.05, na.rm = TRUE),
      n_sig_bh = sum(sub$p_adj_bh < 0.05, na.rm = TRUE),
      frac_sig_bh = mean(sub$p_adj_bh < 0.05, na.rm = TRUE),
      median_odds_ratio_ha = median(sub$odds_ratio_ha, na.rm = TRUE),
      median_n_selected = median(sub$n_selected_genes, na.rm = TRUE),
      p_value = NA_real_
    )
  }

  dsen <- tmp[tmp$method == "DSEN", c("drug", "odds_ratio_ha", "p_value"), drop = FALSE]
  en <- tmp[tmp$method == "EN", c("drug", "odds_ratio_ha", "p_value"), drop = FALSE]
  colnames(dsen)[2:3] <- c("odds_ratio_ha_dsen", "p_dsen")
  colnames(en)[2:3] <- c("odds_ratio_ha_en", "p_en")
  wide <- merge(dsen, en, by = "drug", all = FALSE)
  if (nrow(wide) > 1L) {
    wl <- suppressWarnings(
      wilcox.test(
        log(wide$odds_ratio_ha_dsen),
        log(wide$odds_ratio_ha_en),
        paired = TRUE,
        alternative = "greater",
        exact = FALSE
      )
    )
    med_dsen <- median(wide$odds_ratio_ha_dsen, na.rm = TRUE)
    med_en <- median(wide$odds_ratio_ha_en, na.rm = TRUE)
    add_summary(
      analysis = "paired_wilcoxon",
      selection = sel,
      method = "DSEN_vs_EN",
      group_name = "all_drugs",
      n = nrow(wide),
      n_sig_raw = NA_real_,
      frac_sig_raw = NA_real_,
      n_sig_bh = NA_real_,
      frac_sig_bh = NA_real_,
      median_odds_ratio_ha = med_dsen - med_en,
      median_n_selected = NA_real_,
      p_value = wl$p.value
    )
    if (sel == "freq50") {
      phase2_freq50_p <- wl$p.value
      phase2_freq50_med_dsen <- med_dsen
      phase2_freq50_med_en <- med_en
    }
  }
}

# Pathway-stratified summary for main setting (freq50)
freq50 <- per_drug[per_drug$selection == "freq50", , drop = FALSE]
freq50$primary_pathway <- as.character(primary_pathway[freq50$drug])
freq50$primary_pathway[is.na(freq50$primary_pathway)] <- "Unknown"
pathways <- sort(unique(freq50$primary_pathway))
for (m in c("DSEN", "EN")) {
  for (pw in pathways) {
    sub <- freq50[freq50$method == m & freq50$primary_pathway == pw, , drop = FALSE]
    if (nrow(sub) == 0L) next
    add_summary(
      analysis = "pathway_enrichment",
      selection = "freq50",
      method = m,
      group_name = pw,
      n = nrow(sub),
      n_sig_raw = sum(sub$p_value < 0.05, na.rm = TRUE),
      frac_sig_raw = mean(sub$p_value < 0.05, na.rm = TRUE),
      n_sig_bh = sum(sub$p_adj_bh < 0.05, na.rm = TRUE),
      frac_sig_bh = mean(sub$p_adj_bh < 0.05, na.rm = TRUE),
      median_odds_ratio_ha = median(sub$odds_ratio_ha, na.rm = TRUE),
      median_n_selected = median(sub$n_selected_genes, na.rm = TRUE),
      p_value = NA_real_
    )
  }
}

# OG/TSG directional breakdown from selected cosmic genes (freq50)
if (nrow(gene_detail) > 0L) {
  gd <- gene_detail[gene_detail$selection == "freq50", , drop = FALSE]
  for (m in unique(gd$method)) {
    sub_all <- gd[gd$method == m, , drop = FALSE]
    if (nrow(sub_all) == 0L) next
    expected <- sum(
      (sub_all$is_og == 1L & sub_all$coef_sign == "positive") |
        (sub_all$is_tsg == 1L & sub_all$coef_sign == "negative")
    )
    mismatched <- sum(
      (sub_all$is_og == 1L & sub_all$coef_sign == "negative") |
        (sub_all$is_tsg == 1L & sub_all$coef_sign == "positive")
    )
    add_summary(
      analysis = "og_tsg_breakdown",
      selection = "freq50",
      method = m,
      group_name = "all_signs",
      n = nrow(sub_all),
      n_sig_raw = sum(sub_all$is_og == 1L),
      frac_sig_raw = mean(sub_all$is_og == 1L),
      n_sig_bh = sum(sub_all$is_tsg == 1L),
      frac_sig_bh = mean(sub_all$is_tsg == 1L),
      median_odds_ratio_ha = expected,
      median_n_selected = mismatched,
      p_value = NA_real_
    )
  }
}

phase2_positive <- is.finite(phase2_freq50_p) &&
  phase2_freq50_p < 0.05 &&
  phase2_freq50_med_dsen > phase2_freq50_med_en

add_summary(
  analysis = "phase2_decision",
  selection = "freq50",
  method = "DSEN_vs_EN",
  group_name = "decision",
  n = length(drug_names),
  n_sig_raw = ifelse(phase2_positive, 1L, 0L),
  frac_sig_raw = ifelse(phase2_positive, 1, 0),
  n_sig_bh = NA_real_,
  frac_sig_bh = NA_real_,
  median_odds_ratio_ha = phase2_freq50_med_dsen - phase2_freq50_med_en,
  median_n_selected = NA_real_,
  p_value = phase2_freq50_p
)

summary_df <- rows_to_df(summary_rows)
write.csv(summary_df, SUMMARY_FILE, row.names = FALSE)
log_msg("Wrote aggregate summary: ", SUMMARY_FILE)
log_msg(
  "Phase 2 decision (freq50): ",
  ifelse(phase2_positive, "POSITIVE", "NOT POSITIVE"),
  " (Wilcoxon p=", signif(phase2_freq50_p, 4), ")"
)

# ------------------------------------------------------------------------------
# Phase 3: Tissue-specific enrichment (only if Phase 2 positive)
# ------------------------------------------------------------------------------

if (phase2_positive) {
  tissue_alias <- list(
    LUNG = c("lung"),
    STOMACH = c("stomach"),
    CENTRAL_NERVOUS_SYSTEM = c("brain", "peripheral_nervous_system"),
    SKIN = c("skin"),
    HAEMATOPOIETIC_AND_LYMPHOID_TISSUE = c("blood"),
    OVARY = c("ovary"),
    BONE = c("bone"),
    PANCREAS = c("pancreas"),
    BREAST = c("breast"),
    UPPER_AERODIGESTIVE_TRACT = c("head_and_neck"),
    AUTONOMIC_GANGLIA = c("peripheral_nervous_system", "brain"),
    LARGE_INTESTINE = c("colorectal", "small_intestine"),
    OESOPHAGUS = c("esophagus")
  )

  tumour_tokens <- strsplit(
    tolower(ifelse(is.na(cosmic_df$tumour_types_somatic), "", cosmic_df$tumour_types_somatic)),
    ";",
    fixed = TRUE
  )
  tumour_tokens <- lapply(tumour_tokens, trimws)
  names(tumour_tokens) <- cosmic_df$gene_symbol

  tissue_cosmic <- lapply(names(tissue_alias), function(tiss) {
    pats <- tissue_alias[[tiss]]
    genes <- names(tumour_tokens)[vapply(
      tumour_tokens,
      function(tok) any(tok %in% pats),
      logical(1)
    )]
    intersect(genes, universe_genes)
  })
  names(tissue_cosmic) <- names(tissue_alias)

  tissue_rows <- list()
  for (drug in drug_names) {
    res <- all_results[[drug]]
    tcnt <- align_matrix_rows(
      res$dsen_bootstrap$tissue_cnt,
      feature_names,
      paste0(drug, ":tissue_cnt")
    )
    tissue_cols <- intersect(colnames(tcnt), names(tissue_alias))
    for (tiss in tissue_cols) {
      gene_max_cnt <- tapply(tcnt[valid_idx, tiss], gene_factor, max)
      selected <- names(gene_max_cnt)[gene_max_cnt > 100]
      tissue_set <- tissue_cosmic[[tiss]]
      fish_t <- safe_fisher(selected, tissue_set, universe_genes)

      tissue_rows[[length(tissue_rows) + 1L]] <- data.frame(
        drug = drug,
        tissue = tiss,
        mapped_tumour_sites = collapse_tokens(tissue_alias[[tiss]]),
        n_selected_genes = length(selected),
        n_cosmic_in_selected = fish_t$a,
        n_cosmic_in_universe = length(tissue_set),
        odds_ratio = fish_t$odds_ratio,
        odds_ratio_ha = fish_t$odds_ratio_ha,
        p_value = fish_t$p_value,
        cosmic_genes_found = collapse_tokens(intersect(selected, tissue_set)),
        stringsAsFactors = FALSE
      )
    }
  }

  tissue_df <- rows_to_df(tissue_rows)
  tissue_df$p_adj_bh <- ave(
    tissue_df$p_value,
    tissue_df$tissue,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(tissue_df, TISSUE_FILE, row.names = FALSE)
  log_msg("Wrote tissue-specific enrichment: ", TISSUE_FILE)
} else {
  log_msg("Skipping Phase 3 tissue analysis (Phase 2 was not positive).")
}

# ------------------------------------------------------------------------------
# Phase 4: Visualizations
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  source(here("methods", "04_plot_utils.R"))
  library(ggplot2)
  library(gridExtra)
})

# Figure 1: DSEN vs EN enrichment comparison (freq50)
panel_a <- per_drug[per_drug$selection == "freq50" & per_drug$method %in% c("DSEN", "EN"), , drop = FALSE]
if (nrow(panel_a) > 0L) {
  dsen_order <- panel_a[panel_a$method == "DSEN", c("drug", "p_value"), drop = FALSE]
  dsen_order$score <- -log10(pmax(dsen_order$p_value, 1e-300))
  drug_levels <- dsen_order$drug[order(dsen_order$score, decreasing = TRUE)]

  panel_a$drug <- factor(panel_a$drug, levels = rev(drug_levels))
  panel_a$neglog10_p <- -log10(pmax(panel_a$p_value, 1e-300))

  p_a <- ggplot(panel_a, aes(x = neglog10_p, y = drug, color = method)) +
    geom_point(size = 1.2, alpha = 0.85) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("DSEN" = COLORS$dsen, "EN" = COLORS$right)) +
    labs(
      title = "Panel A: Per-drug COSMIC Enrichment (-log10 p-value)",
      x = expression(-log[10](p)),
      y = "Drug",
      color = "Method"
    ) +
    theme_publication(base_size = 8) +
    theme(
      axis.text.y = element_text(size = 4),
      legend.position = "bottom"
    )

  dsen_wide <- panel_a[panel_a$method == "DSEN", c("drug", "odds_ratio_ha", "p_value"), drop = FALSE]
  en_wide <- panel_a[panel_a$method == "EN", c("drug", "odds_ratio_ha", "p_value"), drop = FALSE]
  colnames(dsen_wide)[2:3] <- c("odds_ratio_ha_dsen", "p_dsen")
  colnames(en_wide)[2:3] <- c("odds_ratio_ha_en", "p_en")
  wide <- merge(dsen_wide, en_wide, by = "drug", all = FALSE)
  wide$sig_group <- ifelse(
    wide$p_dsen < 0.05 | wide$p_en < 0.05,
    "Either p < 0.05", "Neither significant"
  )

  p_b <- ggplot(wide, aes(x = odds_ratio_ha_dsen, y = odds_ratio_ha_en, color = sig_group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Either p < 0.05" = COLORS$dsen, "Neither significant" = COLORS$neutral)) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "Panel B: Odds Ratio Comparison (DSEN vs EN)",
      x = "DSEN odds ratio (Haldane-adjusted, log scale)",
      y = "EN odds ratio (Haldane-adjusted, log scale)",
      color = ""
    ) +
    theme_publication(base_size = 10) +
    theme(legend.position = "bottom")

  fig1_file <- file.path(FIG_DIR, "cosmic_enrichment_comparison.pdf")
  pdf(fig1_file, width = 17, height = 12)
  gridExtra::grid.arrange(p_a, p_b, ncol = 2, widths = c(1.25, 1))
  dev.off()
  log_msg("Wrote figure: ", fig1_file)
}

# Figure 2: Directional enrichment summary (if directional rows exist)
if (nrow(gene_detail) > 0L) {
  gd <- gene_detail[
    gene_detail$selection == "freq50" &
      gene_detail$method == "DSEN",
    ,
    drop = FALSE
  ]
  if (nrow(gd) > 0L) {
    gd$category <- ifelse(
      gd$is_og == 1L & gd$coef_sign == "positive",
      "OG with positive coef",
      ifelse(
        gd$is_tsg == 1L & gd$coef_sign == "negative",
        "TSG with negative coef",
        "Other / mismatched"
      )
    )
    gd$primary_pathway[is.na(gd$primary_pathway) | !nzchar(gd$primary_pathway)] <- "Unknown"

    agg <- as.data.frame(
      table(gd$primary_pathway, gd$category),
      stringsAsFactors = FALSE
    )
    colnames(agg) <- c("primary_pathway", "category", "count")
    agg <- agg[agg$count > 0L, , drop = FALSE]

    if (nrow(agg) > 0L) {
      totals <- aggregate(count ~ primary_pathway, data = agg, FUN = sum)
      pathway_order <- totals$primary_pathway[order(totals$count, decreasing = TRUE)]
      agg$primary_pathway <- factor(agg$primary_pathway, levels = rev(pathway_order))

      p_dir <- ggplot(agg, aes(x = primary_pathway, y = count, fill = category)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c(
          "OG with positive coef" = COLORS$dsen,
          "TSG with negative coef" = COLORS$en,
          "Other / mismatched" = COLORS$neutral
        )) +
        labs(
          title = "Directional COSMIC Role Recovery by Target Pathway (DSEN, freq50)",
          x = "Primary target pathway",
          y = "Number of selected COSMIC gene instances",
          fill = ""
        ) +
        theme_publication(base_size = 10) +
        theme(legend.position = "bottom")

      fig2_file <- file.path(FIG_DIR, "cosmic_directional_enrichment.pdf")
      ggsave(fig2_file, p_dir, width = 11, height = 8)
      log_msg("Wrote figure: ", fig2_file)
    }
  }
}

log_msg(strrep("-", 80))
log_msg("COSMIC enrichment pipeline completed.")
log_msg("Outputs:")
log_msg("  ", FEATURE_MAP_FILE)
log_msg("  ", PER_DRUG_FILE)
log_msg("  ", SUMMARY_FILE)
if (file.exists(DIRECTIONAL_FILE)) log_msg("  ", DIRECTIONAL_FILE)
if (file.exists(TISSUE_FILE)) log_msg("  ", TISSUE_FILE)
log_msg("  ", FIG_DIR)
log_msg(strrep("-", 80))

