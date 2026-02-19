################################################################################
## Script:  checks/check_esd.R
## Project: NES-LTER Zooplankton Size Spectra
## Author:  Alexandra Cabanelas
## Created: August 2024
## Updated: February 2026
##
## Purpose: Diagnostic plots and tables for ESD QC.
##          Run after 06_esd.R.
##          Expected ESD ranges (literature):
##            Nauplii:                     50 – 200 µm
##            Small copepods (Oithona):   200 – 500 µm
##            Medium copepods (Acartia):  400 – 800 µm
##            Large copepods (Calanus):   800 – 2000 µm
##            Euphausiids (adults):      2000 – 8000 µm
##            Chaetognaths:              1000 – 4000 µm
##            Appendicularia:             200 – 600 µm
##            Ostracods:                  300 – 1000 µm
##            Pteropods:                  500 – 3000 µm
################################################################################

source(here::here("R", "00_packages.R"))

zp_esd <- readRDS(here("data", "processed", "06_zp_esd.rds"))

## ── 1) Overall ESD distribution ─────────────────────────────────────────────
# Red dashed = absolute bounds (50–10000 µm)
# Orange dotted = typical mesozooplankton range (200–5000 µm)
zp_esd %>%
  filter(!is.na(ESD_um)) %>%
  ggplot(aes(x = ESD_um)) +
  geom_histogram(bins = 60, fill = "#1C5D8A", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 50,    linetype = "dashed", color = "red",    linewidth = 0.7) +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "red",    linewidth = 0.7) +
  geom_vline(xintercept = 200,   linetype = "dotted", color = "orange", linewidth = 0.7) +
  geom_vline(xintercept = 5000,  linetype = "dotted", color = "orange", linewidth = 0.7) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "ESD (µm, log scale)", y = "Count",
       title = "Distribution of ESD values",
       subtitle = "Red dashed = absolute bounds | Orange dotted = typical range (200–5000 µm)") +
  theme_bw()

## ── 2) ESD by taxa, ordered by median ───────────────────────────────────────
# Euphausiids should be largest, nauplii smallest
# Any inversions = suspect weight or ESD
zp_esd %>%
  distinct(taxa_name, stage, ESD_um) %>%
  filter(!is.na(ESD_um)) %>%
  mutate(taxa_stage = paste(taxa_name, stage, sep = " - "),
         taxa_stage = fct_reorder(taxa_stage, ESD_um, .fun = median)) %>%
  ggplot(aes(x = ESD_um, y = taxa_stage)) +
  geom_point(size = 2, color = "#0E9AA7", alpha = 0.7) +
  geom_vline(xintercept = 200,  linetype = "dashed", color = "orange", linewidth = 0.5) +
  geom_vline(xintercept = 5000, linetype = "dashed", color = "orange", linewidth = 0.5) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "ESD (µm, log scale)", y = NULL,
       title = "ESD by taxa/stage (ordered by median ESD)",
       subtitle = "Should match known body size ranking — euphausiids largest, nauplii smallest") +
  theme_bw(base_size = 9)

## ── 3) ESD vs body length ────────────────────────────────────────────────────
# ESD should be BELOW the 1:1 red line (elongated animals: ESD < prosome length)
# Points above line = ESD > body length = wrong
# Typical ratio ESD/length: ~0.1–0.5 depending on taxon
zp_esd %>%
  distinct(taxa_name, stage, length_mean_um, ESD_um) %>%
  filter(!is.na(ESD_um), !is.na(length_mean_um)) %>%
  ggplot(aes(x = length_mean_um, y = ESD_um, label = taxa_name)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", linewidth = 0.7) +
  geom_point(alpha = 0.7, size = 2, color = "#1C5D8A") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Measured body length (µm)", y = "Calculated ESD (µm)",
       title = "ESD vs measured body length (log–log)",
       subtitle = "All points should fall BELOW the 1:1 line (ESD < body length)") +
  theme_bw()

## ── 4) ESD / length ratio distribution ──────────────────────────────────────
# Typical range: 0.05–0.60 for elongated mesozooplankton
# >0.6 = ESD unrealistically close to body length (near-spherical assumption breaking down)
zp_esd %>%
  filter(!is.na(esd_length_ratio)) %>%
  ggplot(aes(x = esd_length_ratio)) +
  geom_histogram(bins = 50, fill = "#0E9AA7", color = "white", alpha = 0.8) +
  geom_vline(xintercept = 0.6, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(x = "ESD / body length ratio (log scale)", y = "Count",
       title = "ESD / length ratio distribution",
       subtitle = "Red line = 0.6 (above this = ESD suspiciously close to body length)") +
  theme_bw()

# taxa with high or low ratio — investigate these
zp_esd %>%
  filter(!is.na(esd_length_ratio)) %>%
  distinct(taxa_name, stage, length_mean_um, ESD_um, esd_length_ratio) %>%
  filter(esd_length_ratio > 0.6 | esd_length_ratio < 0.05) %>%
  arrange(desc(esd_length_ratio)) %>%
  print(n = Inf)

## ── 5) ESD vs carbon weight ──────────────────────────────────────────────────
# Smooth power-law relationship expected (tight positive trend on log-log)
# Outliers = suspect weight or conversion chain
zp_esd %>%
  distinct(taxa_name, stage, final_C_ug, ESD_um, final_weight_type) %>%
  filter(!is.na(ESD_um), !is.na(final_C_ug), final_C_ug > 0, ESD_um > 0) %>%
  ggplot(aes(x = ESD_um, y = final_C_ug, label = taxa_name)) +
  geom_point(alpha = 0.7, size = 2, color = "#0E9AA7") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "ESD (µm)", y = "Individual carbon weight (µg C)",
       title = "Individual C weight vs ESD",
       subtitle = "Tight positive relationship expected — outliers suggest conversion issues") +
  theme_bw()

## ── 6) Summary table: median ESD per taxa/stage ─────────────────────────────
zp_esd %>%
  filter(!is.na(ESD_um)) %>%
  group_by(taxa_name, stage) %>%
  summarise(
    n           = n(),
    median_ESD  = round(median(ESD_um), 1),
    min_ESD     = round(min(ESD_um), 1),
    max_ESD     = round(max(ESD_um), 1),
    median_C_ug = round(median(final_C_ug, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(median_ESD) %>%
  print(n = 50)
