# NES-LTER Zooplankton Size Spectra

**Author:** Alexandra C. Cabanelas  
**Project:** LTER Synthesis Working Group: Pelagic Community Structure
**LTER Site:** Northeast U.S. Shelf Long-Term Ecological Research (NES-LTER)  
**Created:** August 2024 | **Last updated:** February 2026

---

## Working Group Context

This repository is part of the a cross-site **LTER Pelagic Synthesis Working Group**: *Interannual variability and long-term change in pelagic community structure across a latitudinal gradient*. The working group brings together four marine LTER sites spanning a wide latitudinal gradient to compare how pelagic communities respond to environmental variability at annual and longer time scales:

| Site | Description |
|------|-------------|
| **NES** - Northeast U.S. Shelf | Rapidly warming temperate shelf in the northwest Atlantic Ocean. Supports productive fisheries and dynamic planktonic food webs shaped by strong seasonal forcing. |
| **CCE** - California Current Ecosystem | Eastern boundary upwelling system off the U.S. West Coast. High spatial variability in ocean conditions supports diverse assemblages across trophic levels. |
| **NGA** - Northern Gulf of Alaska | Subarctic coastal system with strong freshwater influence from glacial runoff. Characterized by a large spring bloom of large-celled phytoplankton and diapausing zooplankton, transitioning to a smaller-bodied summer community. |
| **PAL** - Palmer, Antarctica | Rapidly warming polar system west of the Antarctic Peninsula. Sea ice seasonality is a dominant driver of food web structure across trophic levels from phytoplankton to marine mammals. |

The working group is pursuing three interconnected projects:

1. **Normalized Biomass Size Spectra (NBSS)** = *this repository*
2. **Trophic Amplification**: testing whether climate-driven biomass declines are amplified at higher trophic levels and lower latitudes
3. **Double Integration Hypothesis**: examining how cumulative atmospheric forcing drives long-term variability in pelagic populations

This repository contains the **NES-LTER pipeline** for Project 1. Equivalent pipelines for the other three sites are maintained in separate repositories.

---
## Overview

This repository calculates zooplankton **Normalized Biomass Size Spectra (NBSS)** from NES-LTER net tow abundance data. It converts zooplankton abundance into individual carbon weights and equivalent spherical diameters (ESD), then bins organisms into log<sub>2</sub> size classes to produce NBSS slopes per station and cruise.

Size spectra analysis provides a taxonomically aggregated view of community structure: the slope of the NBSS reflects how biomass is distributed across body sizes, and shifts in that slope can signal changes in resource availability, trophic efficiency, and ecosystem state. Comparing NBSS slopes across the four LTER sites is a core goal of the synthesis working group.

---

## Data Sources

| Dataset | Source |
|---|---|
| Zooplankton abundance (335 &micro;m bongo net) | [EDI: knb-lter-nes.25.2](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nes.25.2) (v2) |
| Published mean body lengths | Compiled from literature (see `data/raw/mean_lengths.csv`) |
| Length-weight regressions | Compiled from literature (see `data/raw/length_weight_regressions.csv`) |
| Published mean weights | Compiled from literature (see `data/raw/mean_weights.csv`) |
| DW-to-C conversion factors | Compiled from literature (see `data/raw/dw_to_c_conversions.csv`) |

---

## Repository Structure

```
nes_zp_size_spectra/
├── scripts/
│   ├── 00_packages.R              # load all required packages
│   ├── 01_load_data.R             # read raw CSVs
│   ├── 02_tidy_data.R             # merge staged/unstaged abundance data
│   ├── 03_lengths.R               # add published body lengths
│   ├── 04_weights.R               # L-W regressions + weight QC + gap-filling
│   ├── 05_weight_conversions.R    # WW->DW->C conversions
│   ├── 06_biomass.R               # calculate carbon biomass
│   ├── 07_esd.R                   # calculate ESD
│   ├── 08_nbss.R                  # bin + normalize + fit NBSS slopes + plots
│   └── diagnostic_plots/
│       ├── check_weights.R        # diagnostic plots for weight QC
│       ├── check_esd.R            # diagnostic plots for ESD QC
│       └── check_biomass.R        # diagnostic plots for biomass QC
├── data/
│   ├── raw/                       # input CSVs (not tracked by git)
│   └── processed/                 # intermediate .rds files (not tracked by git)
├── outputs/                       # figures and final tables
├── .gitignore
└── README.md
```

---

## How to Run

1. Clone this repository
2. Place all raw CSV files in `data/raw/` (see Data Sources above)
3. Open R and set your working directory to the project root, or open the `.Rproj` file
4. Run scripts in numbered order:

```r
source("R/00_packages.R")
source("R/01_load_data.R")
source("R/02_tidy_data.R")
source("R/03_lengths.R")
source("R/04_weights.R")
source("R/05_weight_conversions.R")
source("R/06_biomass.R")
source("R/07_esd.R")
source("R/08_nbss.R")
```

Each script saves an `.rds` checkpoint to `data/processed/` so individual steps can be re-run independently without rerunning the full pipeline.

Diagnostic check scripts can be run at any point after their respective pipeline step:

```r
source("R/checks/check_weights.R")   # after step 05
source("R/checks/check_biomass.R")   # after step 06
source("R/checks/check_esd.R")       # after step 07
```

---

## Methods Summary

### Weight calculation
Individual body weights are estimated using published length-weight (L-W) regressions of the form:

```
log10(W) = a + b * log10(L)
```

where *L* is mean length (um) and *W* is individual weight (ug). Regressions return: carbon weight (C), dry weight (DW), or wet weight (WW) depending what was available in the literature. Where no regression exists, published mean weights are used as gap-fillers.

### Weight conversions
- **WW -> DW:** multiply by 0.25 (DW = 25% of WW)
- **DW -> C:** multiply by taxa-specific %C/DW where available; fallback = 40%

### ESD calculation
Each organism is approximated as a sphere. ESD is calculated from weight via:

```
DW ->  WW (/ 0.20)  ->  volume (/ 1.05 g/cm<sup>3)  ->  ESD = (6V/pi)^(1/3)
```

### NBSS
Organisms are binned into log<sub>2</sub> (octave) ESD bins. Biomass per bin is normalized by bin width (mg C m<sup>-2</sup> &micro;m<sup>-1</sup>). The NBSS slope is estimated by fitting a log10-log10 linear regression per station:

```
log10(B_norm) ~ log10(ESD_midpoint)
```

---

## Notes

- Net tows use 335 um mesh. Organisms < ~200 um are undersampled. Drop-off at the small end of the size spectrum is expected.

---

## Dependencies

```r
install.packages(c("here", "tidyverse", "broom", "scales", "forcats"))
```

R version used: 
<!-- Add R ver here -->

---

## Citation

If you use this code, please cite:
<!-- Add citation here -->
> Cabanelas, A. (*year*). NES-LTER Zooplankton Size Spectra pipeline. GitHub. [URL]

---

## License

<!-- Add license here -->
