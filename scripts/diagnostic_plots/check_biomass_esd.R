################################################################################
## Script:  diagnostic_plots/check_biomass_esd.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Diagnostic plots for biomass and ESD QC.
##          Run after 04_biomass_esd.R.
##
##  Plots:
##    1. Biomass distribution
##    2. Top 20 taxa by total biomass
##    3. ESD vs body length 
##    4. Individual carbon weight vs ESD 
##    5. Biomass vs abundance
##    6. Biomass by cruise
##    7. ESD by length 
##
## Inputs:  data/processed/04_zp_biomass_esd.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp <- readRDS(here("data", "processed", "04_zp_biomass_esd.rds"))

## ------------------------------------------ ##
#     1) Biomass distribution -----
## ------------------------------------------ ##
zp %>%
  filter(!is.na(biomass_C_mgC_m2), biomass_C_mgC_m2 > 0) %>%
  ggplot(aes(x = biomass_C_mgC_m2)) +
  geom_histogram(bins = 60, fill = "#1C5D8A", color = "white", alpha = 0.8) +
  scale_x_log10(labels = scales::comma) +
  labs(x = "Biomass (mg C m⁻², log scale)", y = "Count") +
  theme_bw()

## ------------------------------------------ ##
#     2) Top 20 taxa by total biomass -----
## ------------------------------------------ ##
zp %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  group_by(taxa_name) %>%
  summarise(total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE),
            .groups = "drop") %>%
  slice_max(total_mgC_m2, n = 20) %>%
  mutate(taxa_name = fct_reorder(taxa_name, total_mgC_m2)) %>%
  ggplot(aes(x = total_mgC_m2, y = taxa_name)) +
  geom_col(fill = "#0E9AA7") +
  scale_x_log10(labels = scales::comma) +
  labs(
    x        = expression("Total biomass (mg C m"^{-2}~", log scale)"),
    y        = NULL,
    title    = "Top 20 taxa by total biomass",
    subtitle = "Summed across all stations and cruises"
  ) +
  theme_bw()

## ------------------------------------------ ##
#     3) ESD vs body length -----
## ------------------------------------------ ##
zp %>%
  mutate(length_um = coalesce(C_length_um, length_TL_um, length_CL_um,
                              length_PL_um, length_SL_um, length_TrL_um)) %>%
  filter(!is.na(length_um), !is.na(ESD_um)) %>%
  distinct(taxa_name, stage, length_um, ESD_um) %>%
  ggplot(aes(x = length_um, y = ESD_um)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(
    x        = "Body length (µm)",
    y        = "ESD (µm)"
  ) +
  theme_bw()

## ------------------------------------------ ##
#     4) Individual carbon weight vs ESD -----
## ------------------------------------------ ##
zp %>%
  distinct(taxa_name, stage, final_C_ug, ESD_um) %>%
  filter(!is.na(final_C_ug), !is.na(ESD_um), final_C_ug > 0, ESD_um > 0) %>%
  ggplot(aes(x = ESD_um, y = final_C_ug, label = taxa_name)) +
  geom_point(alpha = 0.7, size = 2, color = "#0E9AA7") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  geom_text(size = 2.5, vjust = -0.5, check_overlap = TRUE) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "ESD (µm)",
       y = "Individual carbon weight (µg C)") +
  theme_bw()

## ------------------------------------------ ##
#     5) Biomass vs abundance  -----
## ------------------------------------------ ##
zp %>%
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
  labs(x = "Total abundance (ind / 10 m²)", y = "Total biomass (mg C m⁻²)") +
  theme_bw()

## ------------------------------------------ ##
#     6) Biomass by cruise  -----
## ------------------------------------------ ##
zp %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  group_by(cruise, station) %>%
  summarise(total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = cruise, y = total_mgC_m2)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "#1C5D8A", color = "white", alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5, color = "#0E9AA7") +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "Total biomass (mg C m⁻², log scale)") +
  theme_bw()

## ------------------------------------------ ##
#     7) ESD by length  -----
## ------------------------------------------ ##
zp %>%
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
       title = "ESD by taxa/stage (ordered by median ESD)") +
  theme_bw(base_size = 9)
