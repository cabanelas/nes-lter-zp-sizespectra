################################################################################
## Script:  05_weight_conversions.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Convert all individual weights to a Carbon (ug C).
##          Step 1 (WW -> DW): A few L-W regressions return wet weight (WW).
##            Convert WW -> DW using factor of 0.25.
##          Step 2 (DW -> C): Convert DW to carbon using taxa-specific %C/DW
##            where available, otherwise 40%.
##
##          final_weight_ug   = weight after WW->DW conversion (DW or C)
##          final_weight_type = weight type of final_weight_ug
##          final_C_ug        = per-individual carbon weight (used for biomass)
##
## Inputs:  data/processed/04_zp_weights.rds
## Outputs: data/processed/05_zp_converted.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp <- readRDS(here("data", "processed", "04_zp_weights.rds"))

## --- WW -> DW conversion ---------------------------------------------------
# Most regressions return DW or C; a few WW
# WW -> DW: multiply by 0.25 (DW ~25% of WW for most mesozooplankton)
# Published/overridden weights are already in their stated type — do not convert
#need to check gelatinous
zp <- zp %>%
  mutate(
    final_weight_ug = case_when(
      weight_source %in% c("published weight (L-W override)",
                           "published weight (no equation)") ~ weight_ug,
      lw_reg_weight_type == "WW"                            ~ weight_ug * 0.25,
      TRUE                                                  ~ weight_ug
    ),
    final_weight_type = case_when(
      weight_source %in% c("published weight (L-W override)",
                           "published weight (no equation)") ~ lit_weight_type,
      lw_reg_weight_type == "WW"                            ~ "DW",
      TRUE                                                  ~ lw_reg_weight_type
    ),
    final_weight_note = case_when(
      weight_source == "published weight (L-W override)"              ~ "L-W override: published weight used",
      weight_source == "published weight (no equation)"               ~ "gap-filled: published weight used",
      weight_source == "missing (L-W overridden, no lit)"             ~ "missing: no weight available",
      weight_source == "L-W regression" & lw_reg_weight_type == "WW"  ~ "L-W regression: WW->DW conversion (x0.25)",
      weight_source == "L-W regression"                               ~ "L-W regression",
      TRUE                                                            ~ "missing"
    )
  )

## --- DW -> C conversion  ---------------------------------------------------
# Convert DW to carbon for biomass calculations
# Uses taxa-specific %C/DW from DW_to_C_convert where available
# Fallback: 40% C per DW (general mesozooplankton value)
# Note: gelatinous taxa may need a different conversion — flagged for review

zp <- zp %>%
  mutate(
    final_C_ug = case_when(
      final_weight_type == "C"                             ~ final_weight_ug,                             # already carbon
      final_weight_type == "DW" & !is.na(c_per_dw_percent) ~ final_weight_ug * (c_per_dw_percent / 100),  # taxa-specific
      final_weight_type == "DW" &  is.na(c_per_dw_percent) ~ final_weight_ug * 0.40,                      # fallback 40%
      TRUE                                                 ~ final_weight_ug
    ),
    final_C_note = case_when(
      final_weight_type == "DW" & !is.na(c_per_dw_percent) ~
        paste0(final_weight_note, "; DW->C (taxa-specific: ", c_per_dw_percent, "%)"),
      final_weight_type == "DW" &  is.na(c_per_dw_percent) ~
        paste0(final_weight_note, "; DW->C (fallback: 40%)"),
      TRUE ~ final_weight_note
    ),
    final_C_type = case_when(
      final_weight_type == "DW" ~ "C",
      TRUE                      ~ final_weight_type
    )
  )

## --- Checks ----------------------------------------------------------------
# final_C_type should be all "C" or NA/missing
message("final_C_type breakdown:")
print(count(zp, final_C_type))

message("final_weight_type breakdown:")
print(count(zp, final_weight_type))

message("Records with missing final_C_ug: ",
        sum(is.na(zp$final_C_ug)))

## --- Save ------------------------------------------------------------------
saveRDS(zp, here("data", "processed", "05_zp_converted.rds"))
message("05_weight_conversions.R complete — saved to data/processed/05_zp_converted.rds")
