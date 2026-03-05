################################################################################
## Script:  02_tidy_data.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Tidy and merge zooplankton abundance data.
##          - Pivot staged abundance data from wide to long format
##          - Combine staged and unstaged data, avoiding double-counting
##          - Remove fish, unidentified plankton, and Foraminifera
##          - Add net tow depth (net_max_depth_m) from bongo metadata
## Inputs:  data/processed/01_raw_data.rds
## Outputs: data/processed/02_zp_merged.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

raw <- readRDS(here("data", "processed", "01_raw_data.rds"))
zp             <- raw$zp
zp_staged      <- raw$zp_staged
bongo_metadata <- raw$bongo_metadata

## ------------------------------------------ ##
#            Stage name mapping -----
## ------------------------------------------ ##
# Maps column names in the staged data to clean stage labels
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

## ------------------------------------------ ##
#            Pivot staged data -----
## ------------------------------------------ ##
zp_long_staged <- zp_staged %>%
  pivot_longer(
    cols      = "adult_10m2":"unknown_10m2",
    names_to  = "stage_code",
    values_to = "abundance_10m2"
  ) %>%
  mutate(stage = stage_mapping[stage_code]) %>%
  filter(abundance_10m2 > 0) %>%
  select(-ends_with("_count"), -conc_10m2)

## ------------------------------------------ ##
#            Combine staged + unstaged -----
## ------------------------------------------ ##
# Exclude from unstaged any taxa that appear in the staged dataset
# to prevent double-counting

staged_taxa <- zp_staged %>% distinct(taxa_name) #get names of staged taxa

zp_unstaged <- zp %>%
  rename(abundance_10m2 = conc_10m2) %>%
  anti_join(staged_taxa, by = "taxa_name") %>% #avoid double counting
  mutate(
    stage      = "Not_Staged",
    stage_code = "Not_Staged"
  ) %>%
  filter(abundance_10m2 > 0) %>%
  filter(!taxa_name %in% c("Fish", "Unidentified Plankton", "Foraminifera"))

# --- Combine staged and unstaged 
zp_merged <- bind_rows(zp_unstaged, zp_long_staged)

## ------------------------------------------ ##
#            Tow metadata -----
## ------------------------------------------ ##
# Add net_max_depth_m 

# sample_name == 1 tow == cruise_station_net+cast
# Ring net = "R"; Bongo = "B" 
# NOTE: this analysis (and abundance data) only include samples collected with 
# 335 um mesh Bongo net 

ring_cruises <- c("AR28B", "AR31A", "AR34B", "AR39B", "AR61B", "AR66B")

# --- add sample name column
zp_merged <- zp_merged %>%
  mutate(
    net         = if_else(cruise %in% ring_cruises, "R", "B"),
    sample_name = paste(cruise, station, paste0(net, cast), sep = "_")
  )

bongo_metadata <- bongo_metadata %>%
  mutate(
    net         = if_else(cruise %in% ring_cruises, "R", "B"),
    sample_name = paste(cruise, station, paste0(net, cast), sep = "_")
  )

zp_merged <- zp_merged %>%
  left_join(
    bongo_metadata %>% select(sample_name, net_max_depth_m),
    by = "sample_name"
  ) %>%
  relocate(sample_name, .after = cast) %>%
  relocate(net,          .after = sample_name) %>%
  relocate(net_max_depth_m, .after = day_night) %>%
  select(-net)

# ----------------------------------------------------------------------------
message("Total rows in zp_merged: ", nrow(zp_merged))
message("Unique taxa:             ", n_distinct(zp_merged$taxa_name))
message("Unique stages:           ", n_distinct(zp_merged$stage))
message("Unique cruises:          ", n_distinct(zp_merged$cruise))

## ------------------------------------------ ##
#            Save -----
## ------------------------------------------ ##
saveRDS(zp_merged, here("data", "processed", "02_zp_merged.rds"))

message("02_tidy_data.R complete — saved to data/processed/02_zp_merged.rds")
