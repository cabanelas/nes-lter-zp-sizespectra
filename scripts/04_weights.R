################################################################################
## Script:  04_weights.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate individual body weights using published L-W regressions.
##          - Apply log10 L-W regression: log10(W) = a + b * log10(L)
##          - Fill gaps with published mean weights where no regression exists
##          - Build published weight lookup for comparison and gap-filling
##          - Flag and override suspect L-W regression weights
##
## Inputs:  data/processed/03_zp_lengths.rds
##          data/processed/01_raw_data.rds 
## Outputs: data/processed/04_zp_weights.rds
################################################################################

## NEED TO CHECK; some weights in pub weights are in WW 

source(here::here("scripts", "00_packages.R"))

zp_lengths      <- readRDS(here("data", "processed", "03_zp_lengths.rds"))
raw             <- readRDS(here("data", "processed", "01_raw_data.rds"))

drop_cols <- c("_source_primary", "_doi", "_location", 
               "_source_secondary", "_notes")

LW_reg          <- raw$LW_reg          %>% select(-ends_with(drop_cols))
pub_weights     <- raw$pub_weights     %>% select(-ends_with(drop_cols))
DW_to_C_convert <- raw$DW_to_C_convert %>% select(-ends_with(drop_cols))

## --- 4.0 Examine L-W regression per taxa/stage -----------------------------
# how many taxa/stage combos have multiple regressions?
LW_reg %>%
  count(taxa_name, stage) %>%
  filter(n > 1) %>%
  arrange(desc(n))
# show them to me
LW_reg %>%
  group_by(taxa_name, stage, lw_reg_weight_type) %>%
  filter(n() > 1) %>%
  select(taxa_name, stage, lw_reg_weight_type, 
         lw_reg_source_primary, lw_reg_notes)

## --- 4.1 Select best L-W regression per taxa/stage -------------------------
# Priority: C > DW > WW

# some taxa have more than 1 regression, this picks based on weight type
LW_reg <- LW_reg %>%
  mutate(weight_type_priority = case_when(
    lw_reg_weight_type == "C"  ~ 1,
    lw_reg_weight_type == "DW" ~ 2,
    lw_reg_weight_type == "WW" ~ 3,
    TRUE ~ 4
  ))

LW_reg_best <- LW_reg %>%
  group_by(taxa_name, stage) %>%
  slice_min(weight_type_priority, n = 1, 
            with_ties = FALSE) %>% #if 2 regressions have same weight type, pick first one
  ungroup()

# sanity check; should be empty
stopifnot(nrow(LW_reg_best %>% 
                 count(taxa_name, stage) %>% 
                 filter(n > 1)) == 0)

## --- 4.2 Merge regressions with zooplankton data ---------------------------
zp_eqns <- zp_lengths %>%
  left_join(LW_reg_best, by = c("taxa_name", "stage")) %>%
  mutate(across(where(is.character), ~na_if(.x, "")))

# taxa missing L-W reg
zp_eqns %>%
  filter(is.na(slope_b)) %>%
  distinct(taxa_name, stage)

## --- 4.3 Apply L-W regressions ---------------------------------------------
# Equation: log10(W) = a + b * log10(L)
# Lengths in um; weights in ug

apply_lw <- function(length_um, slope_b, int_a, weight_unit) {
  weight_raw <- 10^(int_a + slope_b * log10(length_um))
  weight_ug  <- case_when(
    weight_unit == "ug" ~ weight_raw,
    weight_unit == "mg" ~ weight_raw * 1000,
    TRUE ~ NA_real_
  )
  return(weight_ug)
}

# Calculate weights
zp_eqns <- zp_eqns %>%
  mutate(weight_ug = apply_lw(length_mean_um, slope_b, int_a, lw_reg_weight_unit))

## --- 4.4 Flag length type mismatches ---------------------------------------
# note where length_mean_type (e.g. TL) != lw_reg_length_type (e.g. PL)
zp_eqns <- zp_eqns %>%
  mutate(
    LW_comment = case_when(
      is.na(length_mean_type) | is.na(lw_reg_length_type) ~ NA_character_,
      length_mean_type != lw_reg_length_type ~
        paste0("Length=", length_mean_type, "; LW expects=", lw_reg_length_type),
      TRUE ~ NA_character_
    )
  )

## ------------------------------------------ ##
#       Check weights 
## ------------------------------------------ ##
## --- 4.5 Build published weight lookup -------------------------------------
# Selects the best available published weight per taxa/stage for:
#   (a) comparison with L-W calculated weights
#   (b) gap-filling where no regression exists
# Priority: prefer weight type matching regression output, then C > DW > WW

# Fill Not_Staged gaps with Adult values 
adult_pub_weights <- pub_weights %>%
  filter(stage == "Adult")

pub_weights <- pub_weights %>%
  left_join(adult_pub_weights, by = c("taxa_name", "lit_weight_type"),
            suffix = c("", "_adult")) %>%
  mutate(
    lit_weight_ug = ifelse(
      stage == "Not_Staged" & is.na(lit_weight_ug),
      lit_weight_ug_adult,
      lit_weight_ug
    )
  ) %>%
  select(-ends_with("_adult"))

pub_weight_lookup <- zp_eqns %>%
  distinct(taxa_name, stage, lw_reg_weight_type) %>%
  left_join(pub_weights,     by = c("taxa_name", "stage")) %>%
  left_join(DW_to_C_convert, by = c("taxa_name", "stage")) %>%
  group_by(taxa_name, stage) %>%
  mutate(priority = case_when(
    !is.na(lw_reg_weight_type) & lit_weight_type == lw_reg_weight_type ~ 0,
    lit_weight_type == "C"  ~ 1,
    lit_weight_type == "DW" ~ 2,
    lit_weight_type == "WW" ~ 3,
    TRUE ~ 4
  )) %>%
  slice_min(priority, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-priority)

# join lookup back to full dataset
zp_eqns <- zp_eqns %>%
  left_join(pub_weight_lookup, by = c("taxa_name", "stage")) %>%
  select(-c(lw_reg_weight_type.y, stage_code.x.y, stage_code.y.y)) %>%
  rename(stage_code         = stage_code.x.x,
         lw_reg_weight_type = lw_reg_weight_type.x)

## --- 4.6 Gap-fill: use published weight when regression is missing ---------
zp_eqns <- zp_eqns %>%
  mutate(
    weight_ug          = coalesce(weight_ug, lit_weight_ug),
    lw_reg_weight_type = coalesce(lw_reg_weight_type, lit_weight_type),
    lw_reg_weight_unit = coalesce(lw_reg_weight_unit, "ug")
  )

## --- 4.7 Convert published weights to calculated w type for comparison -----
# some published weights may not match the weight type obtained from L-W reg...
zp_eqns <- zp_eqns %>%
  mutate(
    lit_weight_comparable_ug = case_when(
      lw_reg_weight_type == lit_weight_type                ~ lit_weight_ug,
      lit_weight_type == "DW" & lw_reg_weight_type == "C"  ~ lit_weight_ug * (c_per_dw_percent / 100),
      lit_weight_type == "C"  & lw_reg_weight_type == "DW" ~ lit_weight_ug / (c_per_dw_percent / 100),
      lit_weight_type == "WW" & lw_reg_weight_type == "DW" ~ lit_weight_ug * 0.20,
      lit_weight_type == "WW" & lw_reg_weight_type == "C"  ~ lit_weight_ug * 0.20 * (c_per_dw_percent / 100),
      TRUE ~ NA_real_
    ),
    comparison_note = case_when(
      lw_reg_weight_type == lit_weight_type ~ "same type",
      is.na(lit_weight_comparable_ug)       ~ "no convert",
      TRUE ~ paste(lit_weight_type, "->", lw_reg_weight_type)
    )
  )

## ── 4.8  Override suspect L-W weights with published values ─────────────────
# These taxa had L-W regressions giving biologically unreasonable weights
# (verified by comparison with published means in checks/check_weights.R)
bad_taxa <- c(
  "Thysanoessa gregaria",
  "Euphausia krohnii",
  "Euphausiacea",
  "Thysanoessa longicaudata",
  "Lucifer",
  "Oikopleura"
)

zp_eqns <- zp_eqns %>%
  mutate(
    # track weight source before overwriting
    weight_source = case_when(
      taxa_name %in% bad_taxa & !is.na(lit_weight_ug) ~ "published weight (L-W override)",
      taxa_name %in% bad_taxa &  is.na(lit_weight_ug) ~ "missing (L-W overridden, no lit)",
      !is.na(slope_b)                                 ~ "L-W regression",
      is.na(slope_b) & !is.na(lit_weight_ug)          ~ "published weight (no equation)",
      TRUE                                            ~ "missing"
    ),
    weight_ug = case_when(
      taxa_name %in% bad_taxa ~ lit_weight_ug,
      TRUE ~ weight_ug
    )
  )

## ── 4.9  Weight source summary ──────────────────────────────────────────────
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

## --- Save
saveRDS(zp_eqns, here("data", "processed", "04_zp_weights.rds"))
message("04_weights.R complete — saved to data/processed/04_zp_weights.rds")
