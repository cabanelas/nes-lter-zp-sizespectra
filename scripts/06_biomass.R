################################################################################
## Script:  06_biomass.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate carbon biomass per taxon per sample
##          Biomass = individual carbon weight (ug C) * abundance (ind / 10 m2)
##
##            biomass_C_ug_10m2  = ug C / 10m2  
##            biomass_C_mgC_m2   = mg C / m2      
##
## Inputs:  data/processed/05_zp_converted.rds
## Outputs: data/processed/06_zp_biomass.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp_eqns <- readRDS(here("data", "processed", "05_zp_converted.rds"))

## --- Calculate carbon biomass ----------------------------------------------
zp_esd <- zp_esd %>%
  mutate(
    biomass_C_ug_10m2 = abundance_10m2 * final_C_ug,    # ug C / 10m2
    biomass_C_mgC_m2  = (biomass_C_ug_10m2 / 10) * 1e-3 # mg C / m2
  )

## --- Checks ----------------------------------------------------------------
# total biomass per station
biomass_check <- zp_esd %>%
  group_by(cruise, station) %>%
  summarise(
    total_mgC_m2 = sum(biomass_C_mgC_m2, na.rm = TRUE),
    n_taxa       = n_distinct(taxa_name),
    .groups = "drop"
  ) %>%
  mutate(flag = case_when(
    total_mgC_m2 < 1   ~ "low (<1 mg C m-2)",
    total_mgC_m2 > 100 ~ "high (>100 mg C m-2)",
    TRUE               ~ "ok"
  ))

message("Biomass per station flag summary:")
print(count(biomass_check, flag))

message("\nTop 10 stations by total biomass (mg C m-2):")
biomass_check %>%
  arrange(desc(total_mgC_m2)) %>%
  print(n = 10)

message("\nBiomass range across all individual records:")
zp_esd %>%
  filter(!is.na(biomass_C_mgC_m2)) %>%
  summarise(
    min_mgC_m2  = min(biomass_C_mgC_m2),
    max_mgC_m2  = max(biomass_C_mgC_m2),
    mean_mgC_m2 = mean(biomass_C_mgC_m2)
  ) %>%
  print()

## --- Save ------------------------------------------------------------------
saveRDS(zp_esd, here("data", "processed", "06_zp_biomass.rds"))
message("06_biomass.R complete — saved to data/processed/06_zp_biomass.rds")
