################################################################################
## Script:  04_biomass_esd.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate carbon biomass and Equivalent Spherical Diameter (ESD)
##          for each taxon/stage/sample
##
##  Biomass:
##    biomass_C_ug_10m2  = final_C_ug  *  abundance_10m2         [ug C / 10 m2]
##    biomass_C_mgC_m2   = biomass_C_ug_10m2 / 10 * 1e-3         [mg C / m2]
##
##  ESD: each organism approximated as a sphere:
##    Conversion chain:
##      final_WW_ug (ug)  ->  WW (g)       [* 1e-6]
##      WW (g)            ->  volume (cm3) [/ 1.05 g cm-3; zooplankton density]
##      volume (cm3)      ->  ESD (cm)     [(6V / pi)^(1/3)]
##      ESD (cm)          ->  ESD (um)     [* 1e4]
##
##    Zooplankton density assumed = 1.05 g cm-3 (slightly denser than seawater)
##
## Inputs:  data/processed/03_zp_weights.rds
## Outputs: data/processed/04_zp_biomass_esd.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp_weights <- readRDS(here("data", "processed", "03_zp_weights.rds"))

## ------------------------------------------ ##
#            Biomass -----
## ------------------------------------------ ##
zp_biomass_esd <- zp_weights %>%
  mutate(
    biomass_C_ug_10m2 = abundance_10m2 * final_C_ug,    # ug C / 10 m2
    biomass_C_mgC_m2  = (biomass_C_ug_10m2 / 10) * 1e-3 # mg C / m2
  )

## ------------------------------------------ ##
#            ESD -----
## ------------------------------------------ ##
zp_biomass_esd <- zp_biomass_esd %>%
  mutate(
    ww_g       = final_WW_ug * 1e-6,          # ug -> g
    volume_cm3 = ww_g / 1.05,                 # g -> cm3 (density = 1.05 g/cm3)
    ESD_cm     = (6 * volume_cm3 / pi)^(1/3), # sphere: d = (6V/pi)^(1/3)
    ESD_um     = ESD_cm * 1e4                 # cm -> um
  )

## ------------------------------------------ ##
#            Checks -----
## ------------------------------------------ ##

# --- Biomass range ---
zp_biomass_esd %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  summarize(
    min  = min(biomass_C_mgC_m2),
    mean = mean(biomass_C_mgC_m2),
    max  = max(biomass_C_mgC_m2),
    n_missing = sum(is.na(biomass_C_mgC_m2))
  ) %>%
  print()

# --- ESD range ---
zp_biomass_esd %>%
  filter(!is.na(ESD_um)) %>%
  summarize(
    min_um    = round(min(ESD_um), 1),
    max_um    = round(max(ESD_um), 1),
    n_missing = sum(is.na(ESD_um))
  ) %>%
  print()

zp_biomass_esd %>%
  distinct(taxa_name, stage, ESD_um) %>%
  filter(!is.na(ESD_um)) %>%
  mutate(esd_flag = case_when(
    ESD_um < 50    ~ "too small (<50 um)",
    ESD_um > 10000 ~ "too large (>10000 um)",
    TRUE           ~ "ok"
  ))

## ------------------------------------------ ##
#            Save -----
## ------------------------------------------ ##
saveRDS(zp_biomass_esd, here("data", "processed", "04_zp_biomass_esd.rds"))

message("04_biomass_esd.R complete — saved to data/processed/04_zp_biomass_esd.rds")
