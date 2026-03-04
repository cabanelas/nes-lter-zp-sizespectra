################################################################################
## Script:  01_load_data.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Load all raw data files required.  
##          Outputs are saved as .rds for use in downstream scripts.
##
## Inputs (data/raw/):
##   - nes-lter-zp-abundance-335um-unstaged10m2.csv  (EDI: knb-lter-nes.25.2)
##   - nes-lter-zp-abundance-335um-staged10m2.csv    (EDI: knb-lter-nes.25.2)
##   - nes-lter-zooplankton-tow-metadata-v2.csv      (EDI:
##   - mean_lengths.csv
##   - length_weight_regressions.csv
##   - mean_weights.csv
##   - dw_to_c_conversions.csv
##
## Outputs (data/processed/):
##   - 01_raw_data.rds   (named list of all loaded data frames)
################################################################################

source(here::here("scripts", "00_packages.R"))

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

# --- Zooplankton abundance data (EDI: knb-lter-nes.25.2) ---
# 335 um mesh bongo net; abundances in individuals per 10 m2
zp <- read_csv(here("data", "raw", #v2 abundance package
                    "nes-lter-zp-abundance-335um-unstaged10m2.csv"))

zp_staged <- read_csv(here("data", "raw",
                           "nes-lter-zp-abundance-335um-staged10m2.csv"))

# --- Net tow metadata (zooplankton inventory package) ---
bongo_metadata <- read_csv(here("data", "raw",
                                "nes-lter-zooplankton-tow-metadata-v2.csv"))

# --- Zooplankton morphometrics 
drop_cols <- c("_source_primary", "_doi", "_location", 
               "_source_secondary", "_notes")

# --- Published mean body lengths (um) ---
lengths <- read_csv(here("data", "raw", "mean_lengths.csv")) %>%
  select(-ends_with(drop_cols))

# --- Published length-weight regressions ---
# log10(W) = a + b * log10(L)
LW_reg <- read_csv(here("data", "raw", "length_weight_regressions.csv")) %>%
  select(-ends_with(drop_cols))

# --- Published mean individual weights (ug) ---
pub_weights <- read_csv(here("data", "raw", "mean_weights.csv")) %>%
  select(-ends_with(drop_cols))

# --- Dry weight to carbon conversion factors (% C per DW) ---
DW_to_C_convert <- read_csv(here("data", "raw", "dw_to_c_conversions.csv")) %>%
  select(-ends_with(drop_cols))

## ------------------------------------------ ##
#           Quick checks -----
## ------------------------------------------ ##
message("zp rows:           ", nrow(zp))                #23460
message("zp_staged rows:    ", nrow(zp_staged))         #2550
message("lengths rows:      ", nrow(lengths))           #215
message("LW_reg rows:       ", nrow(LW_reg))            #253
message("pub_weights rows:  ", nrow(pub_weights))       #351
message("DW_to_C rows:      ", nrow(DW_to_C_convert))   #178

## ------------------------------------------ ##
#           Save -----
## ------------------------------------------ ##
saveRDS(
  list(
    zp             = zp,
    zp_staged      = zp_staged,
    bongo_metadata  = bongo_metadata,
    lengths        = lengths,
    LW_reg         = LW_reg,
    pub_weights    = pub_weights,
    DW_to_C_convert = DW_to_C_convert
  ),
  here("data", "processed", "01_raw_data.rds")
)

message("01_load_data.R complete — saved to data/processed/01_raw_data.rds")
