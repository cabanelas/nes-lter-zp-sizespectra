################################################################################
## Script:  diagnostic_plots/check_weights.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Diagnostic plots and tables for weight QC.
##          Run after 03_weights.R 
##  Plots:
##    1. Calculated vs published weight (1:1 comparison, facet by weight type)
##    2. C:DW ratio per taxa/stage
##    3. DW:WW ratio per taxa/stage
##    4. Weight ~ length (log-log)
##
## Inputs:  data/processed/03_zp_weights.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp <- readRDS(here("data", "processed", "03_zp_weights.rds"))

## ------------------------------------------ ##
#      quick simple ones -----
## ------------------------------------------ ##
ggplot(zp,
       aes(x = DW_length_um, y = DW_weight_ug)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

ggplot(zp,
       aes(x = DW_length_um, y = pub_DW_ug)) +
  geom_point(alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw()

## ------------------------------------------ ##
#      1) Calculated vs published weight -----
## ------------------------------------------ ##
compare_weights <- zp %>%
  distinct(taxa_name, stage,
           C_weight_ug, DW_weight_ug, WW_weight_ug,
           pub_C_ug, pub_DW_ug, pub_WW_ug) %>%
  pivot_longer(
    cols      = c(C_weight_ug, DW_weight_ug, WW_weight_ug),
    names_to  = "weight_type",
    values_to = "calc_weight_ug"
  ) %>%
  mutate(
    weight_type   = sub("_weight_ug$", "", weight_type),
    pub_weight_ug = case_when(
      weight_type == "C"  ~ pub_C_ug,
      weight_type == "DW" ~ pub_DW_ug,
      weight_type == "WW" ~ pub_WW_ug
    )
  ) %>%
  select(taxa_name, stage, weight_type, calc_weight_ug, pub_weight_ug) %>%
  filter(!is.na(calc_weight_ug) | !is.na(pub_weight_ug))

compare_weights %>%
  filter(!is.na(calc_weight_ug), !is.na(pub_weight_ug)) %>%
  ggplot(aes(x = pub_weight_ug, y = calc_weight_ug,
             label = paste(taxa_name, stage))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(size = 2.5, max.overlaps = 15) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ weight_type) +
  labs(
    x        = "Published weight (µg)",
    y        = "Calculated weight (µg)",
    title    = "Calculated vs. published individual weights"
  ) +
  theme_bw()

## ------------------------------------------ ##
#    2) Calculated weights ordered by mean body length
## ------------------------------------------ ##
zp %>%
  distinct(taxa_name, stage,
           length_PL_um, length_TL_um, length_CL_um, length_TrL_um, length_SL_um,
           C_weight_ug, DW_weight_ug, WW_weight_ug) %>%
  mutate(length_any_um = coalesce(length_PL_um, length_TL_um, length_CL_um, 
                                  length_TrL_um, length_SL_um)) %>%
  pivot_longer(cols = ends_with("_weight_ug"),
               names_to = "weight_type", values_to = "weight_ug") %>%
  mutate(weight_type = sub("_weight_ug$", "", weight_type)) %>%
  filter(!is.na(weight_ug)) %>%
  mutate(taxa_stage = paste(taxa_name, stage, sep = " - "),
         taxa_stage = fct_reorder(taxa_stage, length_any_um)) %>%
  ggplot(aes(x = weight_ug, y = taxa_stage, color = weight_type)) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(x = "Calculated weight (µg)", y = NULL, color = "Weight type") +
  theme_bw()

## ------------------------------------------ ##
#      3) C:DW ratio per taxa/stage -----
## ------------------------------------------ ##
zp %>%
  filter(!is.na(C_weight_ug), !is.na(DW_weight_ug), DW_weight_ug > 0) %>%
  distinct(taxa_name, stage, C_weight_ug, DW_weight_ug) %>%
  mutate(
    c_dw_ratio = C_weight_ug / DW_weight_ug,
    label      = paste(taxa_name, stage),
    flag       = c_dw_ratio > 0.7 | c_dw_ratio < 0.2
  ) %>%
  ggplot(aes(x = reorder(label, c_dw_ratio), y = c_dw_ratio, color = flag)) +
  geom_point(size = 2) +
  geom_hline(yintercept = c(0.2, 0.7), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  coord_flip() +
  labs(x = NULL,
       y = "C / DW ratio") +
  theme_bw() +
  theme(legend.position = "none")

## ------------------------------------------ ##
#      4) DW:WW ratio per taxa/stage -----
## ------------------------------------------ ##
zp %>%
  filter(!is.na(DW_weight_ug), !is.na(WW_weight_ug), WW_weight_ug > 0) %>%
  distinct(taxa_name, stage, DW_weight_ug, WW_weight_ug) %>%
  mutate(
    dw_ww_ratio = DW_weight_ug / WW_weight_ug,
    label       = paste(taxa_name, stage),
    flag        = dw_ww_ratio > 0.25 | dw_ww_ratio < 0.05
  ) %>%
  ggplot(aes(x = WW_weight_ug, y = dw_ww_ratio, color = flag, label = label)) +
  geom_point(size = 2) +
  geom_text_repel(data = . %>% filter(flag), size = 2.5) +
  geom_hline(yintercept = c(0.05, 0.25), linetype = "dashed", color = "red") +
  scale_x_log10() +
  scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
  labs(
    x        = "Wet weight (µg)",
    y        = "DW / WW ratio"
  ) +
  theme_bw() +
  theme(legend.position = "none")

## ------------------------------------------ ##
#      5) Weight ~ length (log-log) -----
## ------------------------------------------ ##
zp %>%
  distinct(taxa_name, stage,
           length_PL_um, length_TL_um, length_CL_um, length_TrL_um, length_SL_um,
           C_weight_ug, DW_weight_ug, WW_weight_ug) %>%
  mutate(length_um = coalesce(length_PL_um, length_TL_um, length_CL_um,
                              length_TrL_um, length_SL_um)) %>%
  pivot_longer(
    cols      = ends_with("_weight_ug"),
    names_to  = "weight_type",
    values_to = "weight_ug"
  ) %>%
  mutate(weight_type = sub("_weight_ug$", "", weight_type)) %>%
  filter(!is.na(weight_ug), !is.na(length_um), length_um > 0) %>%
  group_by(weight_type) %>%
  mutate(
    resid = log10(weight_ug) - predict(lm(log10(weight_ug) ~ log10(length_um))),
    z     = (resid - mean(resid)) / sd(resid),
    flag  = abs(z) > 2.5,
    label = paste(taxa_name, stage)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = length_um, y = weight_ug, color = flag, label = label)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = . %>% filter(flag), size = 2.5, max.overlaps = 20) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  facet_wrap(~ weight_type, scales = "free_y") +
  labs(
    x        = "Length (µm)",
    y        = "Weight (µg)",
    title    = "Weight ~ length (log-log)") +
  theme_bw() +
  theme(legend.position = "none")
