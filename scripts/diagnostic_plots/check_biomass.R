################################################################################
## Script:  checks/check_biomass.R
## Project: NES-LTER Zooplankton Size Spectra
## Author:  Alexandra Cabanelas
## Created: August 2024
## Updated: February 2026
##
## Purpose: Diagnostic plots and tables for biomass QC.
##          Run after 07_biomass.R.
##          Typical total mesozooplankton biomass: ~1–100 mg C m⁻²
################################################################################

source(here::here("R", "00_packages.R"))

zp_esd <- readRDS(here("data", "processed", "07_zp_biomass.rds"))

## ── 1) Total biomass per station ─────────────────────────────────────────────
zp_esd %>%
  group_by(cruise, station) %>%
  summarise(total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE), .groups = "drop") %>%
  mutate(station = fct_reorder(as.character(station), total_mgC_m2)) %>%
  ggplot(aes(x = total_mgC_m2, y = station, fill = cruise)) +
  geom_col() +
  geom_vline(xintercept = 1,   linetype = "dashed", color = "orange") +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red") +
  labs(x = "Total biomass (mg C m⁻²)", y = NULL,
       title = "Total zooplankton biomass per station",
       subtitle = "Orange = 1 mg C m⁻² | Red = 100 mg C m⁻² (typical range bounds)",
       fill = "Cruise") +
  theme_bw()

## ── 2) Biomass distribution across all records ───────────────────────────────
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2), biomass_C_mgC_m2 > 0) %>%
  ggplot(aes(x = biomass_C_mgC_m2)) +
  geom_histogram(bins = 60, fill = "#1C5D8A", color = "white", alpha = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Biomass (mg C m⁻², log scale)", y = "Count",
       title = "Distribution of individual taxon biomass records") +
  theme_bw()

## ── 3) Top 20 taxa by total biomass ─────────────────────────────────────────
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  group_by(taxa_name) %>%
  summarise(total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total_mgC_m2, n = 20) %>%
  mutate(taxa_name = fct_reorder(taxa_name, total_mgC_m2)) %>%
  ggplot(aes(x = total_mgC_m2, y = taxa_name)) +
  geom_col(fill = "#0E9AA7") +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Total biomass (mg C m⁻², log scale)", y = NULL,
       title = "Top 20 taxa by total biomass (summed across all stations)") +
  theme_bw()

## ── 4) Biomass vs abundance ──────────────────────────────────────────────────
# Top-right = abundant AND heavy (dominant contributors)
# Top-left  = heavy but rare (large taxa, e.g. euphausiids)
# Bottom-right = abundant but light (small copepods)
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2), abundance_10m2 > 0) %>%
  group_by(taxa_name) %>%
  summarise(
    total_biomass   = sum(biomass_C_mgC_m2, na.rm = TRUE),
    total_abundance = sum(abundance_10m2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = total_abundance, y = total_biomass, label = taxa_name)) +
  geom_point(alpha = 0.7, size = 2, color = "#1C5D8A") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "Total abundance (ind / 10 m²)", y = "Total biomass (mg C m⁻²)",
       title = "Biomass vs abundance by taxon",
       subtitle = "Top-left = large rare taxa | Bottom-right = small abundant taxa") +
  theme_bw()

## ── 5) Individual C weight vs ESD ───────────────────────────────────────────
# Tight power-law relationship expected
# Outliers = suspect weight or ESD calculation
zp_esd %>%
  distinct(taxa_name, stage, final_C_ug, ESD_um) %>%
  filter(!is.na(final_C_ug), !is.na(ESD_um), final_C_ug > 0, ESD_um > 0) %>%
  ggplot(aes(x = ESD_um, y = final_C_ug, label = taxa_name)) +
  geom_point(alpha = 0.7, size = 2, color = "#0E9AA7") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "ESD (µm)", y = "Individual carbon weight (µg C)",
       title = "Individual C weight vs ESD",
       subtitle = "Tight positive power-law expected — outliers suggest conversion issues") +
  theme_bw()

## ── 6) Biomass distribution by cruise ───────────────────────────────────────
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  group_by(cruise, station) %>%
  summarise(total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = cruise, y = total_mgC_m2)) +
  geom_boxplot(fill = "#1C5D8A", color = "white", alpha = 0.7, outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, color = "#0E9AA7") +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Total biomass (mg C m⁻², log scale)",
       title = "Biomass distribution by cruise",
       subtitle = "Each point = one station") +
  theme_bw()

## ── 7) Summary table ─────────────────────────────────────────────────────────
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  group_by(taxa_name, stage) %>%
  summarise(
    n_records     = n(),
    mean_mgC_m2   = round(mean(biomass_C_mgC_m2,  na.rm = TRUE), 4),
    median_mgC_m2 = round(median(biomass_C_mgC_m2, na.rm = TRUE), 4),
    max_mgC_m2    = round(max(biomass_C_mgC_m2,    na.rm = TRUE), 4),
    mean_C_ug_ind = round(mean(final_C_ug,          na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(median_mgC_m2)) %>%
  print(n = 40)
