# Generate revised Figure 4: Target Recovery
# Replaces the ineffective boxplot with a more informative two-panel figure

library(here)
library(tidyverse)
library(gridExtra)
library(grid)

# Read target recovery data
recovery_full <- read_csv(here("results/tables/target_recovery_full.csv"))

# --- Panel A: Summary Bar Chart ---
# Filter to drugs with non-zero recovery
recovered_drugs <- recovery_full %>%
  filter(en_recovered > 0 | dsen_recovered > 0)

# Summary statistics
total_drugs <- nrow(recovery_full)
drugs_with_targets <- sum(!is.na(recovery_full$n_targets) & recovery_full$n_targets > 0)

# Calculate recovery rates
en_any <- sum(recovery_full$en_recovered > 0, na.rm = TRUE)
dsen_any <- sum(recovery_full$dsen_recovered > 0, na.rm = TRUE)
en_full <- sum(recovery_full$en_recovered == 1, na.rm = TRUE)
dsen_full <- sum(recovery_full$dsen_recovered == 1, na.rm = TRUE)

summary_data <- data.frame(
  Method = rep(c("Elastic Net", "DSEN"), 2),
  Metric = c("Any target", "Any target", "Full recovery", "Full recovery"),
  Rate = c(en_any / drugs_with_targets * 100,
           dsen_any / drugs_with_targets * 100,
           en_full / drugs_with_targets * 100,
           dsen_full / drugs_with_targets * 100)
)
summary_data$Metric <- factor(summary_data$Metric, levels = c("Any target", "Full recovery"))
summary_data$Method <- factor(summary_data$Method, levels = c("Elastic Net", "DSEN"))

p_summary <- ggplot(summary_data, aes(x = Metric, y = Rate, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", Rate)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Elastic Net" = "#1f77b4", "DSEN" = "#ff7f0e")) +
  scale_y_continuous(limits = c(0, 8), expand = c(0, 0)) +
  labs(x = NULL, y = "Recovery rate (%)",
       title = "a") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

# --- Panel B: Case Study Table ---
# Create case study data from successful recoveries
case_studies <- recovery_full %>%
  filter(en_recovered > 0 | dsen_recovered > 0) %>%
  mutate(
    Drug = gsub("_repl_.*", "", drug),
    Targets = coalesce(dsen_targets_found, en_targets_found),
    EN = ifelse(en_recovered > 0, ifelse(en_recovered == 1, "Full", "Partial"), "-"),
    DSEN = ifelse(dsen_recovered > 0, ifelse(dsen_recovered == 1, "Full", "Partial"), "-")
  ) %>%
  select(Drug, Targets, EN, DSEN) %>%
  distinct(Drug, .keep_all = TRUE)

# Add drug mechanism annotations
mechanisms <- data.frame(
  Drug = c("PLX-4720", "SB590885", "Dabrafenib", "AZ628",
           "Nutlin-3a (-)", "Afatinib", "Gefitinib", "Cetuximab",
           "Pictilisib", "Pazopanib", "Axitinib"),
  Mechanism = c("BRAF V600E inhibitor", "BRAF inhibitor", "BRAF inhibitor", "BRAF inhibitor",
                "MDM2 antagonist", "Pan-HER TKI", "EGFR TKI", "EGFR antibody",
                "PI3K inhibitor", "Multi-TKI", "Multi-TKI")
)

case_studies <- case_studies %>%
  left_join(mechanisms, by = "Drug") %>%
  arrange(Targets) %>%
  select(Drug, Mechanism, Targets, EN, DSEN)

# Clean up for display
case_studies$Drug <- gsub(" \\(-\\)", "", case_studies$Drug)

# Create table as grob
table_theme <- ttheme_minimal(
  core = list(fg_params = list(fontsize = 10),
              bg_params = list(fill = c("grey95", "white"), col = "grey80")),
  colhead = list(fg_params = list(fontsize = 11, fontface = "bold"),
                 bg_params = list(fill = "grey85", col = "grey80"))
)

# Rename columns for display
colnames(case_studies) <- c("Drug", "Mechanism", "Target(s)", "EN", "DSEN")

table_grob <- tableGrob(case_studies, rows = NULL, theme = table_theme)

# Add title to table
title_grob <- textGrob("b", x = 0.02, hjust = 0,
                       gp = gpar(fontsize = 14, fontface = "bold"))

# Combine into single figure
pdf(here("doc/nature_methods/figures/fig4_target_recovery.pdf"), width = 10, height = 6)

grid.arrange(
  p_summary,
  arrangeGrob(title_grob, table_grob, ncol = 1, heights = c(0.1, 0.9)),
  ncol = 2,
  widths = c(1, 1.3),
  top = textGrob("Target recovery is rare but biologically coherent when it occurs",
                 gp = gpar(fontsize = 13, fontface = "italic"))
)

dev.off()

cat("Figure 4 regenerated successfully!\n")
cat(sprintf("Summary: %d/%d drugs with any target recovery (EN: %d, DSEN: %d)\n",
            max(en_any, dsen_any), drugs_with_targets, en_any, dsen_any))
