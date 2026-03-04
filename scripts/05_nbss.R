################################################################################
## Script:  05_nbss.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate Normalized Biomass Size Spectra (NBSS).
##          - log2 (octave) ESD bins
##          - Sum biomass per station * bin; normalize by bin width
##          - Fit log10-log10 linear regression per station (NBSS slope)
##          - Plot spectra per station, mean ± SD per cruise, slope distributions
##
##  NBSS normalization:
##    B_norm = biomass_bin / bin_width_um    [mg C m-2 um-1]
##
##  Regression per station:
##    log10(B_norm) ~ log10(ESD_midpoint)
##    slope = NBSS slope (expected range ~ -1 to -2; typically -1.2 to -1.5)
##
## Note:  Net tows use 335 um mesh; organisms < ~200 um are undersampled.
##        Drop-off at the small end of the size spectrum is expected.
##
## Inputs:  data/processed/04_zp_biomass_esd.rds
## Outputs: data/processed/05_nbss.rds
##          data/processed/05_nbss_slopes.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp  <- readRDS(here("data", "processed", "04_zp_biomass_esd.rds"))

## ------------------------------------------ ##
#    Define log2 (octave) ESD bins -----
## ------------------------------------------ ##
# octave bins spanning the mesozooplankton range (300-10000 um)
log2_breaks <- seq(log2(300), log2(10000), by = 1) #by=0.5 for half-octave bins
breaks_um   <- 2^log2_breaks

## ------------------------------------------ ##
#    Assign records to bins -----
## ------------------------------------------ ##

zp_binned <- zp %>%
  filter(
    !is.na(ESD_um), !is.na(biomass_C_mgC_m2),
    ESD_um > 0, biomass_C_mgC_m2 > 0
  ) %>%
  mutate(
    esd_bin      = cut(ESD_um, breaks = breaks_um,
                       include.lowest = TRUE, right = FALSE),
    bin_min_um   = breaks_um[as.integer(esd_bin)],
    bin_max_um   = breaks_um[as.integer(esd_bin) + 1],
    bin_mid_um   = sqrt(bin_min_um * bin_max_um),  # geometric midpoint
    bin_width_um = bin_max_um - bin_min_um
  )

# bin coverage
print(count(zp_binned, esd_bin) %>% arrange(esd_bin))

## ------------------------------------------ ##
#  Sum biomass per station x bin; normalize -----
## ------------------------------------------ ##

nbss <- zp_binned %>%
  group_by(cruise, station, cast, sample_name, net_max_depth_m,
           esd_bin, bin_min_um, bin_max_um, bin_mid_um, bin_width_um) %>%
  summarise(
    biomass_C_mgC_m2_bin = sum(biomass_C_mgC_m2, na.rm = TRUE),
    n_taxa               = n_distinct(taxa_name),
    .groups = "drop"
  ) %>%
  mutate(
    B_norm_mgC_m2_um = biomass_C_mgC_m2_bin / bin_width_um, # normalize by bin width
    log10_mid        = log10(bin_mid_um),
    log10_Bnorm      = log10(B_norm_mgC_m2_um)
  ) %>%
  filter(is.finite(log10_Bnorm))  # drop empty/zero bins

# bins per station
nbss %>%
  group_by(cruise, station) %>%
  summarise(n_bins = n(), .groups = "drop") %>%
  pull(n_bins) %>%
  summary() %>%
  print()

## ------------------------------------------ ##
#  Fit NBSS slope per station -----
## ------------------------------------------ ##

nbss_slopes <- nbss %>%
  group_by(cruise, station) %>%
  filter(n() >= 3) %>%
  summarise(
    tidy(lm(log10_Bnorm ~ log10_mid)),
    n_bins = n(),
    .groups = "drop"
  ) %>%
  filter(term == "log10_mid") %>%
  rename(
    slope    = estimate,
    slope_se = std.error,
    slope_p  = p.value
  ) %>%
  select(cruise, station, slope, slope_se, slope_p, n_bins)

summary(nbss_slopes$slope)

## ------------------------------------------ ##
#  Plots -----
## ------------------------------------------ ##

# --- 1) NBSS per station, color by cruise ---
ggplot(
  nbss,
  aes(x = bin_mid_um, y = B_norm_mgC_m2_um,
      group = interaction(cruise, station),
      color = cruise)) +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "ESD (µm)",
       y = expression("Normalized biomass (mg C m"^{-2}~"µm"^{-1}~")"),
       color = "Cruise") +
  theme_bw()

# --- 2) Mean NBSS ± 1 SD per cruise ---
nbss %>%
  group_by(cruise, bin_mid_um) %>%
  summarise(
    mean_Bnorm = mean(B_norm_mgC_m2_um, na.rm = TRUE),
    sd_Bnorm   = sd(B_norm_mgC_m2_um,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = bin_mid_um, y = mean_Bnorm, color = cruise, fill = cruise)) +
  geom_ribbon(aes(ymin = mean_Bnorm - sd_Bnorm,
                  ymax = mean_Bnorm + sd_Bnorm),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(x = "ESD (µm)",
       y = expression("Normalized biomass (mg C m"^{-2}~"µm"^{-1}~")"),
       title = "Mean NBSS by cruise (\u00b1 1 SD)",
       color = "Cruise", fill  = "Cruise") +
  theme_bw()

# --- 3) NBSS slope distribution per cruise ---
ggplot(nbss_slopes, aes(x = slope, fill = cruise)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = -1.2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "darkred") +
  labs(x = "NBSS slope",
       y = "Count") +
  theme_bw()

## ------------------------------------------ ##
#            Save -----
## ------------------------------------------ ##
saveRDS(nbss,        here("data", "processed", "05_nbss.rds"))
saveRDS(nbss_slopes, here("data", "processed", "05_nbss_slopes.rds"))

message("05_nbss.R complete — saved to data/processed/")
