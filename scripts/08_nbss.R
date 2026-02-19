################################################################################
## Script:  08_nbss.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Calculate and plot the Normalized Biomass Size Spectrum (NBSS).
##          - Bin organisms into log2 (octave) ESD bins
##          - Sum biomass per station × bin, normalize by bin width
##          - Fit log10-log10 linear regression per station to get NBSS slope
##          - Plot spectra per station, mean +- SD per cruise, slope distributions
##
##          NBSS normalization: B_norm = biomass_bin / bin_width_um
##          Units: mg C m-2 um-1
##
##          Expected NBSS slope: approximately -1 to -2 (typically -1.2 to -1.5)
##
## Note:    Net tows w 335 um mesh; undersample organisms < ~200 um.
##          Drop-off at small end of spectrum is expected
##
## Inputs:  data/processed/06_zp_biomass.rds
## Outputs: data/processed/08_nbss.rds
##          data/processed/08_nbss_slopes.rds
################################################################################

source(here::here("scripts", "00_packages.R"))
library(broom)

zp_esd <- readRDS(here("data", "processed", "07_zp_biomass.rds"))

## --- 8.1  Define log2 (octave) ESD bins ------------------------------------
# 1-octave bins spanning full mesozooplankton range (~50–10000 um)
# Increase by = 0.5 for finer resolution if needed
log2_breaks <- seq(log2(50), log2(10000), by = 1)
breaks_um   <- 2^log2_breaks

## --- 8.2  Assign each record to a bin --------------------------------------
zp_binned <- zp_esd %>%
  filter(!is.na(ESD_um), !is.na(biomass_C_mgC_m2),
         ESD_um > 0, biomass_C_mgC_m2 > 0) %>%
  mutate(
    esd_bin      = cut(ESD_um, breaks = breaks_um, include.lowest = TRUE, right = FALSE),
    bin_min_um   = breaks_um[as.integer(esd_bin)],
    bin_max_um   = breaks_um[as.integer(esd_bin) + 1],
    bin_mid_um   = sqrt(bin_min_um * bin_max_um),   # geometric midpoint
    bin_width_um = bin_max_um - bin_min_um
  )

message("Bin coverage check:")
print(count(zp_binned, esd_bin) %>% arrange(esd_bin))

## --- 8.3  Sum biomass per station * bin, then normalize --------------------
nbss <- zp_binned %>%
  group_by(cruise, station, cast, esd_bin,
           bin_min_um, bin_max_um, bin_mid_um, bin_width_um) %>%
  summarise(
    biomass_C_mgC_m2_bin = sum(biomass_C_mgC_m2, na.rm = TRUE),
    n_taxa               = n_distinct(taxa_name),
    .groups = "drop"
  ) %>%
  mutate(
    B_norm_mgC_m2_um = biomass_C_mgC_m2_bin / bin_width_um,  # normalize by bin width
    log10_mid        = log10(bin_mid_um),
    log10_Bnorm      = log10(B_norm_mgC_m2_um)
  ) %>%
  filter(is.finite(log10_Bnorm))  # remove empty bins

message("Bins per station summary:")
nbss %>%
  group_by(cruise, station) %>%
  summarise(n_bins = n(), .groups = "drop") %>%
  pull(n_bins) %>%
  summary() %>%
  print()

## --- 8.4  Fit NBSS slope per station ---------------------------------------
# log10(B_norm) ~ log10(ESD_midpoint)
# slope = NBSS slope; intercept = log10(B at 1 um)
nbss_slopes <- nbss %>%
  group_by(cruise, station) %>%
  filter(n() >= 3) %>%   
  summarise(
    tidy(lm(log10_Bnorm ~ log10_mid)),
    n_bins = n(),
    .groups = "drop"
  ) %>%
  filter(term == "log10_mid") %>%
  rename(slope    = estimate,
         slope_se = std.error,
         slope_p  = p.value) %>%
  select(cruise, station, slope, slope_se, slope_p, n_bins)

message("\nNBSS slope summary:")
summary(nbss_slopes$slope)

# flag suspect slopes
nbss_slopes %>%
  mutate(flag = case_when(
    slope > -0.5  ~ "too flat (suspect)",
    slope < -2.5  ~ "too steep (suspect)",
    TRUE          ~ "ok"
  )) %>%
  count(flag) %>%
  print()

nbss_slopes %>%
  mutate(flag = case_when(
    slope > -0.5 ~ "too flat",
    slope < -2.5 ~ "too steep",
    TRUE         ~ "ok"
  )) %>%
  filter(flag != "ok") %>%
  print(n = Inf)

## --- 8.5 Plot: all station spectra -----------------------------------------
ggplot(nbss, aes(x = bin_mid_um, y = B_norm_mgC_m2_um,
                 group = interaction(cruise, station),
                 color = cruise)) +
  geom_line(alpha = 0.4, linewidth = 0.6) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::comma) +
  labs(
    x     = "ESD (um)",
    y     = "Normalized biomass (mg C m-2 um-1)",
    title = "Normalized Biomass Size Spectrum — all stations",
    color = "Cruise"
  ) +
  theme_bw()

## --- 8.6  Plot: mean NBSS +- SD per cruise ---------------------------------
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
  labs(
    x        = "ESD (µm)",
    y        = "Normalized biomass (mg C m-2 um-1)",
    title    = "Mean NBSS by cruise (+- 1 SD)",
    color    = "Cruise",
    fill     = "Cruise"
  ) +
  theme_bw()

## --- 8.7 Plot: NBSS slope distribution per cruise --------------------------
ggplot(nbss_slopes, aes(x = slope, fill = cruise)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = -1.2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "darkred") +
  labs(
    x        = "NBSS slope",
    y        = "Count",
    title    = "Distribution of NBSS slopes by cruise",
    subtitle = "Dashed lines = typical expected range (-1.2 to -1.5)"
  ) +
  theme_bw()

## --- Save ------------------------------------------------------------------
saveRDS(nbss,        here("data", "processed", "08_nbss.rds"))
saveRDS(nbss_slopes, here("data", "processed", "08_nbss_slopes.rds"))

message("08_nbss.R complete — saved to data/processed/")
