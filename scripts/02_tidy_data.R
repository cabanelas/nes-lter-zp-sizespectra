################################################################################
## Script:  02_tidy_data.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Tidy and merge zooplankton abundance data.
##          - Pivot staged abundance data from wide to long format
##          - Combine staged and unstaged data, avoiding double-counting
##          - Remove fish and unidentified plankton
##
## Inputs:  data/processed/01_raw_data.rds
## Outputs: data/processed/02_zp_merged.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

raw <- readRDS(here("data", "processed", "01_raw_data.rds"))
zp             <- raw$zp
zp_staged      <- raw$zp_staged
lengths        <- raw$lengths

# --- Stage name mapping -----------------------------------------------------
stage_mapping <- c(
  "adult_10m2"      = "Adult", # left hand side = colnames in zp staged
  "c5_10m2"         = "CV",
  "c4_10m2"         = "CIV",
  "c3_10m2"         = "CIII",
  "c2_10m2"         = "CII",
  "c1_10m2"         = "CI",
  "cryptopia_10m2"  = "Cryptopia",
  "furcilia_10m2"   = "Furcilia",
  "calyptopis_10m2" = "Calyptopis",
  "nauplius_10m2"   = "Nauplius",
  "unknown_10m2"    = "Not_Staged"
)

# --- Pivot staged data to long format ---------------------------------------
zp_long_staged <- zp_staged %>%
  pivot_longer(
    cols      = "adult_10m2":"unknown_10m2",
    names_to  = "stage_code",
    values_to = "abundance_10m2"
  ) %>%
  mutate(stage = stage_mapping[stage_code]) %>%
  filter(abundance_10m2 > 0) %>%
  select(-ends_with("_count"), -conc_10m2)

# --- Unstaged data: exclude taxa already in staged data ---------------------
# prevents double-counting taxa that appear in both datasets
staged_taxa <- zp_staged %>% distinct(taxa_name) #get names of staged taxa

zp_unstaged <- zp %>%
  rename(abundance_10m2 = conc_10m2) %>%
  anti_join(staged_taxa, by = "taxa_name") %>% #avoid double counting
  mutate(
    stage      = "Not_Staged",
    stage_code = "Not_Staged"
  ) %>%
  filter(abundance_10m2 > 0) %>%
  filter(!taxa_name %in% c("Fish", "Unidentified Plankton"))

# --- Combine staged and unstaged --------------------------------------------
zp_merged <- bind_rows(zp_unstaged, zp_long_staged)

# ----------------------------------------------------------------------------
message("Total rows in zp_merged: ", nrow(zp_merged))
message("Unique taxa:             ", n_distinct(zp_merged$taxa_name))
message("Unique stages:           ", n_distinct(zp_merged$stage))
message("Unique cruises:          ", n_distinct(zp_merged$cruise))

# --- Save -------------------------------------------------------------------
saveRDS(zp_merged, here("data", "processed", "02_zp_merged.rds"))
saveRDS(lengths,   here("data", "processed", "02_lengths_clean.rds"))

message("02_tidy_data.R complete — saved to data/processed/")
