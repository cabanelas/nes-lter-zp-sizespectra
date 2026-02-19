################################################################################
## Script:  checks/check_weights.R
## Project: NES-LTER Zooplankton Size Spectra
## Author:  Alexandra Cabanelas
## Created: August 2024
## Updated: February 2026
##
## Purpose: Diagnostic plots and tables for weight QC.
##          Run after 04_weights.R and 05_weight_conversions.R.
##          Use these plots to identify taxa with unreasonable weights
##          before proceeding to ESD and biomass calculations.
################################################################################

source(here::here("R", "00_packages.R"))

zp_eqns <- readRDS(here("data", "processed", "05_zp_converted.rds"))

## ── 1) Calculated vs. published weights (1:1 comparison) ────────────────────
# Points on red dashed line = perfect agreement
# Color = weight type (C, DW, WW)
zp_eqns %>%
  filter(!is.na(weight_ug), !is.na(lit_weight_ug)) %>%
  distinct(taxa_name, stage, weight_ug, lit_weight_ug, lw_reg_weight_type) %>%
  ggplot(aes(x = lit_weight_ug, y = weight_ug,
             color = lw_reg_weight_type, label = taxa_name)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Published weight (µg)", y = "Calculated weight (µg)",
       title = "Calculated vs. published weights",
       color = "Weight type") +
  theme_bw()

## ── 2) Calculated weights ordered by mean body length ───────────────────────
# Sanity check: weights should increase generally left to right
zp_eqns %>%
  distinct(taxa_name, stage, length_mean_um, weight_ug) %>%
  filter(!is.na(weight_ug)) %>%
  arrange(length_mean_um) %>%
  mutate(taxa_stage = paste(taxa_name, stage, sep = " - "),
         taxa_stage = fct_reorder(taxa_stage, length_mean_um)) %>%
  ggplot(aes(x = weight_ug, y = taxa_stage)) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Calculated weight (µg)", y = NULL,
       title = "Taxa ordered by length (smallest to largest)",
       subtitle = "Weights should increase generally left to right") +
  theme_bw()

## ── 3) Weight distributions for the 40 most abundant taxa ───────────────────
top_taxa <- zp_eqns %>%
  filter(!is.na(weight_ug)) %>%
  group_by(taxa_name) %>%
  summarise(tot_abund = sum(abundance_10m2, na.rm = TRUE), .groups = "drop") %>%
  slice_max(tot_abund, n = 40) %>%
  pull(taxa_name)

zp_eqns %>%
  filter(taxa_name %in% top_taxa, !is.na(weight_ug)) %>%
  mutate(taxa_name = fct_reorder(taxa_name, length_mean_um,
                                 .fun = mean, na.rm = TRUE)) %>%
  ggplot(aes(x = taxa_name, y = weight_ug)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_y_log10() +
  coord_flip() +
  labs(x = NULL, y = "Predicted individual weight (µg)",
       title = "Weight distributions — top 40 most abundant taxa",
       subtitle = "Ordered by mean body length") +
  theme_bw()

## ── 4) Median log10 ratio: predicted vs. published ──────────────────────────
# 0 = perfect agreement | ±0.3 ≈ factor of 2
# Flags taxa where L-W regression diverges most from literature
comp <- zp_eqns %>%
  filter(!is.na(weight_ug), !is.na(lit_weight_ug)) %>%
  mutate(log10_ratio = log10(weight_ug / lit_weight_ug))

comp_sum <- comp %>%
  group_by(taxa_name, stage) %>%
  summarise(med_log10_ratio = median(log10_ratio, na.rm = TRUE),
            n = n(), .groups = "drop")

ggplot(comp_sum, aes(x = reorder(taxa_name, med_log10_ratio),
                     y = med_log10_ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  coord_flip() +
  labs(x = NULL,
       y = "Median log10(predicted / published)",
       title = "Divergence of predicted from published weights",
       subtitle = "0 = perfect agreement | ±0.3 ≈ factor of 2") +
  theme_bw()

comp_sum %>%
  arrange(desc(abs(med_log10_ratio))) %>%
  print(n = 30)

## ── 5) Ratio table: same weight type only (apples-to-apples) ────────────────
# Ratio = 1: perfect | >1: predicted > published | <1: predicted < published
zp_eqns %>%
  distinct(taxa_name, stage, weight_ug, lit_weight_ug,
           lw_reg_weight_type, lit_weight_type) %>%
  filter(!is.na(weight_ug), !is.na(lit_weight_ug),
         lw_reg_weight_type == lit_weight_type,
         weight_ug > 0, lit_weight_ug > 0) %>%
  mutate(ratio = weight_ug / lit_weight_ug) %>%
  arrange(desc(ratio)) %>%
  select(taxa_name, stage, lw_reg_weight_type, weight_ug, lit_weight_ug, ratio) %>%
  print(n = 40)

## ── 6) Absolute weight flags ─────────────────────────────────────────────────
# Flag biologically unlikely weights for mesozooplankton
zp_eqns %>%
  filter(!is.na(weight_ug)) %>%
  mutate(weight_flag = case_when(
    weight_ug > 10000 ~ "suspiciously heavy (>10 mg)",
    weight_ug < 0.01  ~ "suspiciously light (<0.01 µg)",
    TRUE              ~ NA_character_
  )) %>%
  filter(!is.na(weight_flag)) %>%
  distinct(taxa_name, stage, length_mean_um, weight_ug, weight_flag) %>%
  arrange(desc(weight_ug)) %>%
  print(n = Inf)

## ── 7) Weight source summary ─────────────────────────────────────────────────
zp_eqns %>%
  group_by(weight_source) %>%
  summarise(
    n_taxastage     = n_distinct(paste(taxa_name, stage)),
    n_records       = n(),
    total_abundance = sum(abundance_10m2, na.rm = TRUE)
  ) %>%
  mutate(
    pct_records   = round(100 * n_records / sum(n_records), 1),
    pct_abundance = round(100 * total_abundance / sum(total_abundance), 1)
  ) %>%
  print()
