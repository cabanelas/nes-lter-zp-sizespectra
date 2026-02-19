################################################################################
## Script:  03_lengths.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Add published mean body lengths (um) to zooplankton data.
##          - Join lengths by taxa_name + stage
##          - For Not_Staged taxa missing a length, use Adult length
##
## Inputs:  data/processed/02_zp_merged.rds
##          data/processed/02_lengths_clean.rds
## Outputs: data/processed/03_zp_lengths.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp_merged <- readRDS(here("data", "processed", "02_zp_merged.rds"))
lengths   <- readRDS(here("data", "processed", "02_lengths_clean.rds"))

# --- Join lengths to zooplankton data ---------------------------------------
zp_lengths <- zp_merged %>%
  left_join(lengths, by = c("taxa_name", "stage")) %>%
  select(-stage_code.y) %>%
  rename(stage_code = stage_code.x)

# --- if Not_Staged length is missing, use Adult length ----------------------
adult_lengths <- lengths %>%
  filter(stage == "Adult")

zp_lengths <- zp_lengths %>%
  left_join(adult_lengths, by = "taxa_name", 
            suffix = c("", "_adult")) %>%
  mutate(
    length_mean_um = ifelse(
      stage == "Not_Staged" & is.na(length_mean_um),
      length_mean_um_adult,
      length_mean_um
    )
  ) %>%
  select(-ends_with("_adult"))

# --- Checks -----------------------------------------------------------------
# taxa I was not able to find length data for
zp_lengths %>%
  filter(is.na(length_mean_um)) %>%
  distinct(taxa_name, stage) %>% 
  print(n=500)

# full length table 
zp_lengths %>%
  distinct(taxa_name, stage, length_mean_um) %>%
  filter(!is.na(length_mean_um)) %>%
  arrange(length_mean_um) %>%
  print(n = Inf)

# --- Save -------------------------------------------------------------------
saveRDS(zp_lengths, here("data", "processed", "03_zp_lengths.rds"))
message("03_lengths.R complete — saved to data/processed/03_zp_lengths.rds")
