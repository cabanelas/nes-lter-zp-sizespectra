################################################################################
## Project: LTER Pelagic Synthesis Working Group
##          Zooplankton Size Spectra
## LTER Site: NES-LTER
##
## Author:  Alexandra C. Cabanelas
## Affiliation: Woods Hole Oceanographic Institution
## 
## Script:  00_packages.R
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Install and load all packages required for the full pipeline.
##          Run this script first before any other script.
##
## Notes:   Uncomment install.packages() lines on first run.
################################################################################

# --- Install packages (first run only) ---
# install.packages("here")
# install.packages("tidyverse")
# install.packages("broom")
# install.packages("scales")
# install.packages("ggrepel") # needed for plots only

# --- Load packages ---
library(here)       # reproducible file paths
library(tidyverse)  # dplyr, ggplot2, tidyr, readr, forcats, stringr, purrr
library(broom)      # tidy model output (NBSS slope fitting)
library(scales)     # axis formatting in ggplot2
library(ggrepel)    # optional: diagnostic plots only

# --- Session info ---
# This project uses renv for package version management.
# To restore the exact package environment used here, run:
#   install.packages("renv")
#   renv::restore()
# Package versions are recorded in renv.lock (tracked by git).
sessionInfo()
