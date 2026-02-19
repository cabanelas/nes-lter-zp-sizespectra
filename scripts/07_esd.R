################################################################################
## Script:  07_esd.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate Equivalent Spherical Diameter (ESD) for each taxon/stage.
##          Each organism is approximated as a sphere. 
##          ESD is derived from weight -> volume -> diameter
##
##          Conversion chain:
##            C (ug) -> DW (ug)  [÷ %C/DW or 40%]
##            DW (ug) -> WW (ug) [÷ 0.20; DW ~20% of WW]
##            WW (ug) -> WW (g)  [× 1e-6]
##            WW (g)  -> V (cm3) [÷ 1.05 g/cm3, zooplankton density]
##            V (cm3) -> ESD (cm)[= (6V/pi)^(1/3)]
##            ESD (cm)-> ESD (um)[× 1e4]
##
## Inputs:  data/processed/06_zp_biomass.rds
## Outputs: data/processed/07_zp_esd.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp_esd <- readRDS(here("data", "processed", "06_zp_biomass.rds"))

## --- Calculate ESD ---------------------------------------------------------
zp_esd <- zp_eqns %>%
  mutate(
    # step 1: get DW
    dw_for_esd_ug = case_when(
      final_weight_type == "DW"                                ~ final_weight_ug,
      final_weight_type == "C" & !is.na(c_per_dw_percent)     ~ final_weight_ug / (c_per_dw_percent / 100),
      final_weight_type == "C" &  is.na(c_per_dw_percent)     ~ final_weight_ug / 0.40,
      TRUE ~ NA_real_
    ), ####CHECK GELATINOUS

    # step 2: DW -> WW (DW ≈ 20% of WW for most mesozooplankton)
    ww_for_esd_ug = dw_for_esd_ug / 0.20,

    # step 3: WW ug -> g, then g -> volume using density = 1.05 g/cm3
    ww_for_esd_g = ww_for_esd_ug * 1e-6, # ug to g: * 10^−6
    volume_cm3   = ww_for_esd_g / 1.05,  # density of zooplankton = g/cm3

    # step 4: ESD from sphere volume formula V = (pi/6) * d3 == d = (6V/pi)^(1/3)
    ESD_cm = (6 * volume_cm3 / pi)^(1/3),

    # step 5: convert ESD from cm to um (1 cm = 1e4 um)
    ESD_um = ESD_cm * 1e4,

    esd_length_ratio = ESD_um / length_mean_um
  )

## --- Checks ----------------------------------------------------------------
# flag biologically unlikely ESD values
esd_flags <- zp_esd %>%
  distinct(taxa_name, stage, ESD_um) %>%
  filter(!is.na(ESD_um)) %>%
  mutate(esd_flag = case_when(
    ESD_um < 50    ~ "too small (<50 um)",
    ESD_um > 10000 ~ "too large (>10000 um)",
    TRUE           ~ "ok"
  ))

message("ESD flag summary:")
print(count(esd_flags, esd_flag))

flagged <- esd_flags %>% filter(esd_flag != "ok")
if (nrow(flagged) > 0) {
  message("Flagged taxa:")
  print(flagged, n = Inf)
}

message("ESD range (um): ",
        round(min(zp_esd$ESD_um, na.rm = TRUE), 1), " – ",
        round(max(zp_esd$ESD_um, na.rm = TRUE), 1))

message("Records with missing ESD: ", sum(is.na(zp_esd$ESD_um)))

## --- Save ------------------------------------------------------------------
saveRDS(zp_esd, here("data", "processed", "07_zp_esd.rds"))
message("07_esd.R complete — saved to data/processed/07_zp_esd.rds")
