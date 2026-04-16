#!/usr/bin/env Rscript
# ==============================================================================
# Figure 4: Case studies (3 panels, 2-row layout)
#
# Top row: Panel a (overview heatmap, full width)
# Bottom row: Panels b and c side by side
#
# Output: doc/nature_methods/figures/fig4_case_studies.pdf
#         width=183mm, height=170mm
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(purrr)
  library(readr)
  library(scales)
  library(tidyr)
})

source(here("methods", "00_paths.R"))
source(here("methods", "01_data_utils.R"))
source(here("experiments", "figure_style_nature_methods.R"))

tr <- read_csv(
  here("results", "tables", "threshold_25_v2", "target_recovery_extended.csv"),
  show_col_types = FALSE
)
targets <- read_csv(here("data", "drug_targets_curated.csv"), show_col_types = FALSE)
gene_mapping <- read_csv(here("data", "gene_name_mapping.csv"), show_col_types = FALSE)
rppa_mapping <- read_csv(here("data", "rppa_gene_mapping.csv"), show_col_types = FALSE)

load(here("results", "dsen_results", "threshold_25_v2", "dsen_bootstrap_results.Rda"))
bootstrap_results <- all_results
rm(all_results)

main_data <- load_main_data("Sanger")

bootstrap_lookup <- setNames(bootstrap_results, map_chr(bootstrap_results, "drug"))

standardize_gene <- function(gene) {
  idx <- match(toupper(gene), toupper(gene_mapping$alias))
  if (!is.na(idx)) {
    return(gene_mapping$hgnc_symbol[idx])
  }
  gene
}

get_feature_modality <- function(feature_name) {
  if (grepl("_mut\\.", feature_name)) return("Mutation")
  if (grepl("_Expr$", feature_name)) return("Expression")
  if (grepl("_CN$", feature_name)) return("CN")
  if (feature_name %in% rppa_mapping$rppa_feature) return("RPPA")
  "Other"
}

get_gene_from_feature <- function(feature_name) {
  if (grepl("_mut\\.", feature_name)) return(gsub("_mut\\..*", "", feature_name))
  if (grepl("_Expr$", feature_name)) return(gsub("_Expr$", "", feature_name))
  if (grepl("_CN$", feature_name)) return(gsub("_CN$", "", feature_name))
  idx <- match(feature_name, rppa_mapping$rppa_feature)
  if (!is.na(idx)) return(rppa_mapping$hgnc_symbol[idx])
  feature_name
}

feature_display_label <- function(feature_name) {
  modality <- get_feature_modality(feature_name)
  gene <- get_gene_from_feature(feature_name)
  prefix <- c(
    "Expression" = "Expr.",
    "Mutation" = "Mut.",
    "CN" = "CN",
    "RPPA" = "RPPA",
    "Other" = "Feat."
  )[[modality]]
  paste(prefix, gene)
}

zscore_rows <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      rep(0, length(v))
    } else {
      (v - mean(v, na.rm = TRUE)) / s
    }
  }))
  z[is.na(z)] <- 0
  nm_clip(z, limit = 3)
}

lineage_palette <- c(
  "Focus tissue" = nm_colors$dsen_tissue,
  "Other tissues" = "#D9D9D9"
)

# ==============================================================================
# Panel a: Overview heatmap
# ==============================================================================

tr50 <- tr %>%
  filter(threshold == 0.5) %>%
  filter(en_any | dsen_shared_any | dsen_tissue_any) %>%
  mutate(
    drug_base = gsub("_repl_.*", "", drug),
    is_repl = grepl("_repl_", drug)
  ) %>%
  group_by(drug_base) %>%
  arrange(is_repl, desc(dsen_combined_any), desc(en_any)) %>%
  slice(1) %>%
  ungroup()

target_lookup <- targets %>%
  filter(!is.na(target_gene), target_gene != "NA") %>%
  mutate(target_gene = vapply(target_gene, standardize_gene, character(1))) %>%
  group_by(drug_name) %>%
  summarize(
    target_label = paste(unique(target_gene), collapse = ", "),
    target_genes = list(unique(target_gene)),
    .groups = "drop"
  )

tr50 <- tr50 %>%
  left_join(target_lookup, by = c("drug_base" = "drug_name")) %>%
  mutate(
    target_label = if_else(
      is.na(target_label),
      "?",
      if_else(nchar(target_label) > 30, paste0(substr(target_label, 1, 27), "..."), target_label)
    ),
    archetype = case_when(
      dsen_shared_any & dsen_tissue_any ~ "Both blocks",
      dsen_shared_any & !dsen_tissue_any ~ "Shared only",
      !dsen_shared_any & dsen_tissue_any ~ "Tissue only",
      en_any & !dsen_shared_any & !dsen_tissue_any ~ "EN only",
      TRUE ~ "Other"
    ),
    archetype = factor(archetype, levels = c("Both blocks", "Shared only", "Tissue only", "EN only"))
  ) %>%
  filter(!is.na(archetype)) %>%
  arrange(archetype, drug_base)

THRESH <- 0.5

modality_records <- map_dfr(seq_len(nrow(tr50)), function(i) {
  drug_name <- tr50$drug[i]
  drug_base <- tr50$drug_base[i]
  target_genes <- tr50$target_genes[[i]]
  boot <- bootstrap_lookup[[drug_name]]
  if (is.null(boot) || boot$status != "success" || is.null(target_genes)) {
    return(tibble())
  }

  nboot_en <- boot$en_bootstrap$nboot
  nboot_dsen <- boot$dsen_bootstrap$nboot
  en_stable <- names(boot$en_bootstrap$cnt)[boot$en_bootstrap$cnt / nboot_en >= THRESH]
  shared_stable <- names(boot$dsen_bootstrap$cnt)[boot$dsen_bootstrap$cnt / nboot_dsen >= THRESH]
  tissue_stable <- rownames(boot$dsen_bootstrap$tissue_cnt)[
    apply(boot$dsen_bootstrap$tissue_cnt / nboot_dsen >= THRESH, 1, any)
  ]

  bind_rows(
    tibble(method = "EN", feature = en_stable),
    tibble(method = "Shared", feature = shared_stable),
    tibble(method = "Tissue", feature = tissue_stable)
  ) %>%
    mutate(
      gene = vapply(feature, get_gene_from_feature, character(1)),
      modality = vapply(feature, get_feature_modality, character(1)),
      drug_base = drug_base
    ) %>%
    filter(toupper(gene) %in% toupper(target_genes))
})

dominant_modality <- modality_records %>%
  filter(modality != "Other") %>%
  count(drug_base, modality, sort = TRUE) %>%
  group_by(drug_base) %>%
  slice(1) %>%
  ungroup() %>%
  select(drug_base, modality)

tr50 <- tr50 %>%
  left_join(dominant_modality, by = "drug_base") %>%
  mutate(modality = replace_na(modality, "Other"))

overview_order <- rev(tr50$drug_base)
tr50$drug_base_f <- factor(tr50$drug_base, levels = overview_order)

tile_df <- tr50 %>%
  select(drug_base_f, drug_base, archetype, target_label, modality, en_any, dsen_shared_any, dsen_tissue_any) %>%
  pivot_longer(
    cols = c(en_any, dsen_shared_any, dsen_tissue_any),
    names_to = "method",
    values_to = "found"
  ) %>%
  mutate(
    method = recode(
      method,
      en_any = "EN",
      dsen_shared_any = "Shared",
      dsen_tissue_any = "Tissue"
    ),
    method_num = c("EN" = 1, "Shared" = 2, "Tissue" = 3)[method],
    fill_value = case_when(
      method == "EN" & found ~ "EN",
      method == "Shared" & found ~ "Shared",
      method == "Tissue" & found ~ "Tissue",
      TRUE ~ "Not found"
    )
  )

target_annot <- tr50 %>%
  transmute(drug_base_f, archetype, target_label, modality)

panel_a <- ggplot(tile_df, aes(x = method_num, y = drug_base_f, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 0.85) +
  geom_text(
    data = target_annot,
    aes(x = 4.15, y = drug_base_f, label = target_label),
    inherit.aes = FALSE,
    hjust = 0,
    family = nm_base_family,
    size = 2.25,
    color = nm_colors$dark
  ) +
  geom_point(
    data = target_annot,
    aes(x = 5.8, y = drug_base_f, color = modality),
    inherit.aes = FALSE,
    shape = 15,
    size = 2.3
  ) +
  facet_grid(archetype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(
    values = c(
      "EN" = nm_colors$en,
      "Shared" = nm_colors$dsen_shared,
      "Tissue" = nm_colors$dsen_tissue,
      "Not found" = nm_colors$not_found
    ),
    guide = "none"
  ) +
  scale_color_manual(values = nm_modality_colors, name = "Target feature class") +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("EN", "Shared", "Tissue"),
    expand = expansion(add = c(0.5, 1.8))
  ) +
  labs(x = NULL, y = NULL, tag = "a") +
  coord_cartesian(clip = "off") +
  nm_theme_compact() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 6.2),
    strip.placement = "outside",
    strip.background = element_rect(fill = alpha(nm_colors$light_gray, 0.95), color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(5, 15, 5, 5)
  )

# ==============================================================================
# Exemplar panel builder (updated palette, more breathing room)
# ==============================================================================

build_exemplar_panel <- function(drug_name, mode, title, tag, focus_tissue = NULL, highlight_gene,
                                 show_heat_legend = TRUE) {
  boot <- bootstrap_lookup[[drug_name]]
  if (is.null(boot) || boot$status != "success") {
    stop("Bootstrap result not found for ", drug_name)
  }

  drug_data <- prepare_drug_data_filtered(main_data, drug_name, min_per_tissue = 25)
  if (is.null(drug_data)) {
    stop("Modeling data not found for ", drug_name)
  }

  if (mode == "shared") {
    freq_vec <- boot$dsen_bootstrap$cnt / boot$dsen_bootstrap$nboot
    coef_vec <- boot$dsen_bootstrap$mean_coefs[, "Shared"]
    block_color <- nm_colors$dsen_combined
  } else {
    freq_vec <- boot$dsen_bootstrap$tissue_cnt[, focus_tissue] / boot$dsen_bootstrap$nboot
    coef_vec <- boot$dsen_bootstrap$mean_coefs[, focus_tissue]
    block_color <- nm_colors$dsen_tissue
  }

  feature_df <- tibble(
    feature = names(freq_vec),
    freq = unname(freq_vec),
    coef = unname(coef_vec[names(freq_vec)]),
    modality = vapply(names(freq_vec), get_feature_modality, character(1)),
    gene = vapply(names(freq_vec), get_gene_from_feature, character(1))
  ) %>%
    mutate(
      is_target = toupper(gene) == toupper(highlight_gene),
      exact_target = feature == highlight_gene
    ) %>%
    filter(freq > 0 | is_target)

  top_features <- feature_df %>%
    arrange(desc(freq), desc(abs(coef)))

  target_feature <- top_features %>%
    filter(exact_target | is_target) %>%
    arrange(desc(exact_target), desc(freq), desc(abs(coef))) %>%
    slice(1) %>%
    pull(feature)

  selected_features <- unique(top_features$feature)
  selected_features <- selected_features[seq_len(min(length(selected_features), 10))]
  if (length(target_feature) == 1 && !(target_feature %in% selected_features)) {
    selected_features[length(selected_features)] <- target_feature
  }

  selected_df <- top_features %>%
    filter(feature %in% selected_features) %>%
    arrange(desc(freq), desc(abs(coef))) %>%
    distinct(feature, .keep_all = TRUE) %>%
    slice_head(n = 10) %>%
    mutate(
      feature_label = vapply(feature, feature_display_label, character(1)),
      label_text = sprintf("%s (%.2f)", feature_label, freq),
      coef_scaled = coef * 1000
    )

  if (!highlight_gene %in% selected_df$gene) {
    warning("Highlighted gene not selected for ", drug_name, " (", highlight_gene, ")")
  }

  selected_df <- selected_df %>%
    mutate(row_id = rev(seq_len(n())))

  available_features <- selected_df$feature[selected_df$feature %in% colnames(drug_data$x)]
  selected_df <- selected_df %>%
    filter(feature %in% available_features) %>%
    mutate(row_id = rev(seq_len(n())))

  ord <- order(drug_data$y)
  x_subset <- drug_data$x[ord, selected_df$feature, drop = FALSE]
  mat <- t(as.matrix(x_subset))
  mode(mat) <- "numeric"
  z_mat <- zscore_rows(mat)
  sample_ids <- colnames(z_mat)

  heat_df <- as_tibble(z_mat, rownames = "feature") %>%
    left_join(selected_df %>% select(feature, row_id), by = "feature") %>%
    pivot_longer(cols = -c(feature, row_id), names_to = "sample_id", values_to = "z_value") %>%
    mutate(sample_rank = match(sample_id, sample_ids))

  response_df <- tibble(
    sample_rank = seq_along(ord),
    response = drug_data$y[ord],
    focus_group = if (is.null(focus_tissue)) {
      "All tissues"
    } else if_else(drug_data$tissue[ord] == focus_tissue, "Focus tissue", "Other tissues")
  )

  highlight_rows <- selected_df %>% filter(is_target)

  bar_plot <- ggplot(selected_df, aes(x = coef_scaled, y = row_id)) +
    geom_vline(xintercept = 0, color = nm_colors$neutral, linewidth = 0.4) +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id + 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id - 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    geom_col(aes(fill = coef_scaled >= 0), width = 0.7, orientation = "y", show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = nm_colors$heat_high, "FALSE" = nm_colors$heat_low)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.5, nrow(selected_df) + 0.5), expand = c(0, 0)) +
    labs(x = expression("Mean coef. (" %*% 10^{-3} * ")"), y = NULL, tag = tag) +
    nm_theme_compact(base_size = 7.8) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.title.x = element_text(size = 6.5),
      axis.text.x = element_text(size = 6),
      plot.margin = margin(4, 2, 2, 4)
    )

  heat_plot <- ggplot(heat_df, aes(x = sample_rank, y = row_id, fill = z_value)) +
    geom_raster() +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id + 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id - 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    scale_fill_gradient2(
      low = nm_colors$heat_low,
      mid = nm_colors$heat_mid,
      high = nm_colors$heat_high,
      midpoint = 0,
      limits = c(-3, 3),
      name = "Feature\nz-score",
      guide = guide_colorbar(
        direction = "vertical",
        barwidth = grid::unit(0.25, "cm"),
        barheight = grid::unit(1.8, "cm")
      )
    ) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.5, nrow(selected_df) + 0.5), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title) +
    nm_theme_compact(base_size = 7.8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
      plot.margin = margin(4, 4, 4, 4)
    )

  label_plot <- ggplot(selected_df, aes(y = row_id)) +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id + 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    geom_hline(
      data = highlight_rows,
      aes(yintercept = row_id - 0.48),
      inherit.aes = FALSE,
      color = block_color,
      linewidth = 0.35
    ) +
    geom_point(aes(x = 0.05, color = modality), shape = 15, size = 2.0, show.legend = FALSE) +
    geom_text(
      aes(
        x = 0.14,
        label = label_text,
        color = modality,
        fontface = if_else(is_target, "bold", "plain")
      ),
      hjust = 0,
      family = nm_base_family,
      size = 2.35,
      show.legend = FALSE
    ) +
    scale_color_manual(values = nm_modality_colors, guide = "none") +
    scale_x_continuous(NULL, breaks = NULL, limits = c(0, 1.1), expand = c(0, 0)) +
    scale_y_continuous(NULL, breaks = NULL, limits = c(0.5, nrow(selected_df) + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = nm_base_family) +
    theme(plot.margin = margin(4, 2, 4, 4))

  response_plot <- ggplot(response_df, aes(x = sample_rank, y = response)) +
    geom_line(color = nm_colors$en, linewidth = 0.45) +
    scale_y_continuous(
      breaks = pretty_breaks(n = 2),
      labels = number_format(accuracy = 0.1)
    ) +
    labs(x = NULL, y = "Activity\narea") +
    nm_theme_compact(base_size = 6.5) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 4, 1, 4),
      axis.title.y = element_text(size = 6)
    )

  response_bundle <- if (is.null(focus_tissue)) {
    response_plot
  } else {
    tissue_strip <- ggplot(response_df, aes(x = sample_rank, y = 1, fill = focus_group)) +
      geom_raster() +
      scale_fill_manual(values = lineage_palette, guide = "none") +
      theme_void(base_family = nm_base_family) +
      theme(plot.margin = margin(0, 4, 4, 4))

    response_plot / tissue_strip + plot_layout(heights = c(1, 0.14))
  }

  # Build a standalone z-score legend as a grob
  legend_plot <- ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_gradient2(
      low = nm_colors$heat_low,
      mid = nm_colors$heat_mid,
      high = nm_colors$heat_high,
      midpoint = 0,
      limits = c(-3, 3),
      name = "Feature\nz-score",
      guide = guide_colorbar(
        direction = "vertical",
        barwidth = grid::unit(0.25, "cm"),
        barheight = grid::unit(1.5, "cm")
      )
    ) +
    theme_void(base_family = nm_base_family) +
    theme(
      legend.position = "left",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 5.5),
      legend.margin = margin(0, 0, 0, 0)
    )
  legend_grob <- cowplot::get_legend(legend_plot)

  main_row <- (bar_plot | heat_plot | label_plot) +
    plot_layout(widths = c(2.0, 7.0, 3.0))
  response_row <- (response_bundle | plot_spacer()) +
    plot_layout(widths = c(9.0, 3.0))

  right_panel <- main_row / response_row +
    plot_layout(heights = c(5.0, 0.7))

  exemplar <- wrap_elements(legend_grob) | right_panel
  exemplar <- exemplar + plot_layout(widths = c(0.05, 1.0))

  exemplar
}

# ==============================================================================
# Build exemplar panels
# ==============================================================================

panel_b <- build_exemplar_panel(
  drug_name = "Cetuximab",
  mode = "shared",
  title = "Cetuximab shared block: EGFR",
  tag = "b",
  highlight_gene = "EGFR",
  show_heat_legend = TRUE
)

panel_c <- build_exemplar_panel(
  drug_name = "Ruxolitinib",
  mode = "tissue",
  focus_tissue = "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
  title = "Ruxolitinib tissue block: JAK2 (hematopoietic)",
  tag = "c",
  highlight_gene = "JAK2",
  show_heat_legend = FALSE
)

# ==============================================================================
# 2-row layout: panel_a on top, panels b+c side by side on bottom
# ==============================================================================

fig4 <- panel_a / panel_b / panel_c +
  plot_layout(heights = c(0.8, 0.9, 0.9))

out_path <- here("doc", "nature_methods", "figures", "fig4_case_studies.pdf")
nm_save(fig4, out_path, width_mm = 178, height_mm = 220)
file.copy(
  out_path,
  here("doc", "nature_methods", "latex", "fig4_case_studies.pdf"),
  overwrite = TRUE
)

cat("Figure 4 saved:", out_path, "\n")
cat(sprintf("  Overview drugs: %d unique recoveries at threshold 0.5\n", nrow(tr50)))
