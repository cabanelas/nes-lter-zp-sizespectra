################################################################################
## Script:  03_weights.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Derive individual body weights for all taxa/stage combinations.
##            final_C_ug  : per-individual carbon weight (ug C); for biomass
##            final_WW_ug : per-individual wet weight (ug WW);   for ESD
##
##          Weight derivation priority:
##
##          Carbon (final_C_ug):
##            1. L-W regression returning C directly
##            2. L-W regression returning DW  *  (taxa-specific %C/DW ÷ 100)
##            3. Published C weight
##            4. Published DW  *  (%C/DW ÷ 100)
##            5. Published WW  *  0.20  *  (%C/DW ÷ 100)
##            6. L-W WW  *  0.20  *  (%C/DW ÷ 100)
##
##          Wet weight (final_WW_ug):
##            1. L-W regression returning WW directly
##            2. L-W DW  ÷  0.20
##            3. L-W C  ÷  (%C/DW ÷ 100)  ÷  0.20
##            4. Published WW
##            5. Published DW  ÷  0.20
##            6. Published C  ÷  (%C/DW ÷ 100)  ÷  0.20
##
##          WW -> DW conversion factor: 0.20 (DW ~ 20% of WW for mesozoop)
##
## Inputs:  data/processed/01_raw_data.rds
##          data/processed/02_zp_merged.rds
## Outputs: data/processed/03_zp_weights.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

raw       <- readRDS(here("data", "processed", "01_raw_data.rds"))
zp_merged <- readRDS(here("data", "processed", "02_zp_merged.rds"))

lengths         <- raw$lengths
LW_reg          <- raw$LW_reg
pub_weights     <- raw$pub_weights
DW_to_C_convert <- raw$DW_to_C_convert %>% select(-stage_code)

## ------------------------------------------ ##
#            STEP 1 - Lengths -----
## ------------------------------------------ ##
# --- Pivot lengths to wide format so each length type has its own column ---
# Columns: length_PL_um, length_TL_um, length_CL_um, length_TrL_um, length_SL_um

lengths_wide <- lengths %>%
  select(taxa_name, stage, length_mean_um, length_mean_type) %>%
  pivot_wider(
    names_from  = length_mean_type,
    values_from = length_mean_um,
    names_glue  = "length_{length_mean_type}_um"
  ) %>%
  select(-length_NA_um)

# --- Taxa/stage with no length data ---
lengths %>%
  filter(is.na(length_mean_um)) %>%
  distinct(taxa_name, stage) %>%
  { if (nrow(.) > 0) { message("Taxa with no length data:"); print(.) }
    else message("All taxa have at least one length value") }

## ------------------------------------------ ##
#            STEP 2 - L-W regressions -----
## ------------------------------------------ ##
# some taxa/stage/weight type have more than 1 regression
# When duplicates exist, keep the equation whose length type matches
# published length type (length_mean_type)

# --- Length types available for each taxa/stage ---
zp_length_types <- lengths %>%
  distinct(taxa_name, stage, length_mean_type)

LW_reg_dedup <- LW_reg %>%
  mutate(
    length_match_priority = if_else(
      paste(taxa_name, stage, lw_reg_length_type) %in%
        paste(zp_length_types$taxa_name,
              zp_length_types$stage,
              zp_length_types$length_mean_type),
      1L, 2L
    )
  ) %>%
  group_by(taxa_name, stage, lw_reg_weight_type) %>%
  slice_min(order_by = length_match_priority, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-length_match_priority)

# --- sanity check - should be empty ---
stopifnot(
  nrow(
    LW_reg_dedup %>%
      group_by(taxa_name, stage, lw_reg_weight_type) %>%
      filter(n() > 1)
  ) == 0
)

## ------------------------------------------ ##
#    STEP 3 - Apply L-W regressions -----
## ------------------------------------------ ##
# calculate C, DW, WW weights
# Equation: log10(W) = a + b * log10(L)
# L in um, W output in ug

calc_weight <- function(wt_type) {
  # get the equations for this weight type
  eqns <- LW_reg_dedup %>% filter(lw_reg_weight_type == wt_type)
  # for each taxa/stage, join only the length that matches the eqn length type
  lengths %>%
    inner_join(eqns, by = c("taxa_name", "stage",
                            "length_mean_type" = "lw_reg_length_type")) %>%
    mutate(weight = 10^(int_a + slope_b * log10(length_mean_um))) %>%
    select(
      taxa_name, stage,
      !!paste0(wt_type, "_weight_ug")   := weight,
      !!paste0(wt_type, "_length_um")   := length_mean_um,
      !!paste0(wt_type, "_length_type") := length_mean_type,
      !!paste0(wt_type, "_eqn")         := lw_reg_eqn
    )
}

# --- join back onto zp_merged ---
zp_eqns <- zp_merged %>%
  left_join(calc_weight("C"),  by = c("taxa_name", "stage")) %>%
  left_join(calc_weight("DW"), by = c("taxa_name", "stage")) %>%
  left_join(calc_weight("WW"), by = c("taxa_name", "stage")) %>%
  left_join(lengths_wide,      by = c("taxa_name", "stage"))

# --- shows taxa missing L-W reg ---
zp_eqns %>%
  filter(is.na(C_weight_ug) & is.na(DW_weight_ug) & is.na(WW_weight_ug)) %>%
  distinct(taxa_name, stage) %>%
  { if (nrow(.) > 0) { message("Taxa with no L-W weight:"); print(.) }
    else message("All taxa have at least one L-W weight") }

## ------------------------------------------ ##
#    STEP 4 - Published weight lookup -----
## ------------------------------------------ ##
# Build a per-taxa/stage lookup with one published weight per type (C, DW, WW)
# for use in gap-filling and diagnostics

pub_weight_lookup <- zp_eqns %>%
  distinct(taxa_name, stage) %>%
  left_join(pub_weights, by = c("taxa_name", "stage")) %>%
  filter(!is.na(lit_weight_type)) %>%
  group_by(taxa_name, stage, lit_weight_type) %>%
  slice(1) %>%
  ungroup() %>%
  select(taxa_name, stage, lit_weight_type, lit_weight_ug) %>%
  pivot_wider(
    names_from  = lit_weight_type,
    values_from = lit_weight_ug,
    names_glue  = "pub_{lit_weight_type}_ug"  # pub_C_ug, pub_DW_ug, pub_WW_ug
  )

# --- Join published weights onto zp_eqns for comparison ---
zp_eqns <- zp_eqns %>%
  left_join(pub_weight_lookup, by = c("taxa_name", "stage")) %>%
  left_join(DW_to_C_convert,  by = c("taxa_name", "stage")) %>%
  mutate(c_per_dw = c_per_dw_percent / 100)

## ------------------------------------------ ##
#  STEP 5 - Get final Carbon weight for biomass  -----
## ------------------------------------------ ##
# Most L-W regressions give DW or C, but a few give WW
# Converting WW -> DW using factor of 0.20
# %C per DW values from DW_to_C_convert 

zp_weights <- zp_eqns %>%
  mutate(
    final_C_ug = case_when(
      # 1. L-W C equation
      !is.na(C_weight_ug)       ~ C_weight_ug,
      # 2. L-W DW -> C 
      !is.na(DW_weight_ug)      ~ DW_weight_ug * c_per_dw,
      # 3. Published C directly
      !is.na(pub_C_ug)          ~ pub_C_ug,
      # 4. Published DW -> C
      !is.na(pub_DW_ug)         ~ pub_DW_ug * c_per_dw,
      # 5. Published WW -> DW -> C
      !is.na(pub_WW_ug)         ~ pub_WW_ug * 0.20 * c_per_dw,
      # 6. L-W WW -> DW -> C
      !is.na(WW_weight_ug)      ~ WW_weight_ug * 0.20 * c_per_dw,
      TRUE                      ~ NA_real_
    ),

    final_C_source = case_when(
      !is.na(C_weight_ug)       ~ "L-W C direct",
      !is.na(DW_weight_ug)      ~ "L-W DW -> C",
      !is.na(pub_C_ug)          ~ "pub C direct",
      !is.na(pub_DW_ug)         ~ "pub DW -> C",
      !is.na(pub_WW_ug)         ~ "pub WW -> DW -> C",
      !is.na(WW_weight_ug)      ~ "L-W WW -> DW -> C",
      TRUE                      ~ "missing"
    )
  )

## ------------------------------------------ ##
#  STEP 6 - Get WW for ESD  -----
## ------------------------------------------ ##

zp_weights <- zp_weights %>%
  mutate(
    final_WW_ug = case_when(
      # 1. Direct L-W WW equation
      !is.na(WW_weight_ug) ~ WW_weight_ug,
      # 2. L-W DW / 0.20 -> WW
      !is.na(DW_weight_ug) ~ DW_weight_ug / 0.20,
      # 3. L-W C -> DW -> WW
      !is.na(C_weight_ug)  ~ (C_weight_ug / c_per_dw) / 0.20,
      # 4. Published WW directly
      !is.na(pub_WW_ug)    ~ pub_WW_ug,
      # 5. Published DW / 0.20
      !is.na(pub_DW_ug)    ~ pub_DW_ug / 0.20,
      # 6. Published C -> DW -> WW
      !is.na(pub_C_ug)     ~ (pub_C_ug / c_per_dw) / 0.20,
      TRUE ~ NA_real_
    ),

    final_WW_source = case_when(
      !is.na(WW_weight_ug) ~ "L-W WW direct",
      !is.na(DW_weight_ug) ~ "L-W DW -> WW",
      !is.na(C_weight_ug)  ~ "L-W C -> DW -> WW",
      !is.na(pub_WW_ug)    ~ "pub WW direct",
      !is.na(pub_DW_ug)    ~ "pub DW -> WW",
      !is.na(pub_C_ug)     ~ "pub C -> DW -> WW",
      TRUE                 ~ "missing"
    )
  )

## ------------------------------------------ ##
#  STEP 7 - Diagnostics  -----
## ------------------------------------------ ##
# --- Carbon source coverage ---
message("\n=== Carbon weight source coverage ===")
zp_weights %>%
  group_by(final_C_source) %>%
  summarize(
    n_taxastage = n_distinct(paste(taxa_name, stage)),
    n_records   = n(),
    pct_records = round(100 * n() / nrow(zp_weights), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records)) %>%
  print()

# --- WW source coverage ---
message("\n=== Wet weight source coverage ===")
zp_weights %>%
  group_by(final_WW_source) %>%
  summarize(
    n_taxastage = n_distinct(paste(taxa_name, stage)),
    n_records   = n(),
    pct_records = round(100 * n() / nrow(zp_weights), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records)) %>%
  print()

# --- Taxa missing carbon ---
missing_C <- zp_weights %>%
  filter(is.na(final_C_ug)) %>%
  distinct(taxa_name, stage)

if (nrow(missing_C) > 0) {
  message("\nTaxa missing final_C_ug:")
  print(missing_C, n = Inf)
} else {
  message("\nAll taxa have a final_C_ug value")
}

## ------------------------------------------ ##
#            Save -----
## ------------------------------------------ ##
saveRDS(zp_weights, here("data", "processed", "03_zp_weights.rds"))

message("03_weights.R complete — saved to data/processed/03_zp_weights.rds")
