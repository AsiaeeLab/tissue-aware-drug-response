suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

nm_base_family <- "sans"

nm_colors <- list(
  primary_blue = "#1B4F72",
  primary_orange = "#E67E22",
  light_blue = "#85C1E9",
  red_accent = "#C0392B",
  gray = "#7F8C8D",
  light_gray = "#F2F3F4",
  dark = "#2C3E50",

  # Semantic aliases

  en = "#7F8C8D",
  dsen_combined = "#1B4F72",
  dsen_shared = "#85C1E9",
  dsen_tissue = "#E67E22",
  leakage = "#C0392B",
  honest = "#1B4F72",
  win = "#1B4F72",
  loss = "#C0392B",
  neutral = "#7F8C8D",
  background = "#F2F3F4",
  not_found = "#F2F3F4",
  heat_low = "#1B4F72",
  heat_mid = "#F8F8F8",
  heat_high = "#E67E22"
)

nm_modality_colors <- c(
  "Expression" = "#85C1E9",
  "Mutation" = "#E67E22",
  "CN" = "#1B4F72",
  "RPPA" = "#C0392B",
  "Other" = "#7F8C8D"
)

nm_wrap <- function(x, width = 18) {
  vapply(
    x,
    function(one) paste(strwrap(one, width = width), collapse = "\n"),
    character(1)
  )
}

nm_theme <- function(base_size = 8, base_family = nm_base_family) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, color = nm_colors$dark),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.25),
      axis.text = element_text(size = base_size - 1, color = nm_colors$dark),
      axis.title = element_text(size = base_size, color = nm_colors$dark),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2, color = nm_colors$gray),
      plot.caption = element_text(size = base_size - 1, color = nm_colors$gray),
      plot.tag = element_text(size = base_size + 2, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width = grid::unit(0.42, "cm")
    )
}

nm_theme_compact <- function(base_size = 8, base_family = nm_base_family) {
  nm_theme(base_size = base_size, base_family = base_family) +
    theme(
      legend.spacing.x = grid::unit(0.12, "cm"),
      legend.box.spacing = grid::unit(0.08, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

nm_save <- function(plot, path, width_mm, height_mm) {
  ggsave(
    filename = path,
    plot = plot,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    device = cairo_pdf
  )
}

nm_clip <- function(x, limit = 3) {
  pmax(pmin(x, limit), -limit)
}
