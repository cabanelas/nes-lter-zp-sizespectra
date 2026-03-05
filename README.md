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
| Net tow metadata | [EDI: knb-lter-nes.24.2](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nes.24.2) (v2) |
| Published mean body lengths | Compiled from literature (see `data/raw/mean_lengths.csv`) |
| Length-weight regressions | Compiled from literature (see `data/raw/length_weight_regressions.csv`) |
| Published mean weights | Compiled from literature (see `data/raw/mean_weights.csv`) |
| DW-to-C conversion factors | Compiled from literature (see `data/raw/dw_to_c_conversions.csv`) |  

See `data/raw/data_dictionary.csv` for column descriptions of all raw input files.  
See `data/raw/references.bib` for literature sources for L-W regressions, mean weights, mean lengths, and DW-to-C conversions (Zotero export).

---

## Repository Structure

```
nes_zp_size_spectra/
├── scripts/
│   ├── 00_packages.R               # load all required packages
│   ├── 01_load_data.R              # read raw CSVs
│   ├── 02_tidy_data.R              # merge staged/unstaged abundance data
│   ├── 03_weights.R                # lengths + L-W regressions + weight conversions
│   ├── 04_biomass_esd.R            # carbon biomass + ESD calculation
│   ├── 05_nbss.R                   # bin + normalize + fit slopes + plots
│   └── diagnostic_plots/
│       ├── check_weights.R         # weight QC plots 
│       └── check_biomass_esd.R     # biomass and ESD QC plots
│       └── check_weight_coverage.R # shows which taxa + stage have which weights
├── data/
│   ├── raw/                        # input CSVs
│   └── processed/                  # intermediate .rds checkpoints (not tracked by git)
├── outputs/                        # figures and final tables
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
source("scripts/00_packages.R")
source("scripts/01_load_data.R")
source("scripts/02_tidy_data.R")
source("scripts/03_weights.R")
source("scripts/04_biomass_esd.R")
source("scripts/05_nbss.R")
```

Each script saves an `.rds` checkpoint to `data/processed/` so individual steps can be re-run independently without rerunning the full pipeline.

Diagnostic check scripts can be run at any point after their respective pipeline step:

```r
source("scripts/diagnostic_plots/check_weights.R")            # after 03_weights.R
source("scripts/diagnostic_plots/check_weight_coverage.R")    # after 03_weights.R
source("scripts/diagnostic_plots/check_biomass_esd.R")        # after 04_biomass_esd.R
```

---

## Methods Summary

### Abundance data

Zooplankton abundance data are from 335 &micro;m mesh bongo net tows, in units of individuals per 10 m<sup>2</sup>. Staged taxa (primarily *Calanus* copepodite stages) are handled separately to preserve stage-specific length and weight information, then merged with unstaged taxa for the full analysis. Fish, unidentified plankton, and Foraminifera are excluded.

### Weight calculation
Individual body weights are estimated using published length-weight (L-W) regressions of the form:

```
log10(W) = a + b * log10(L)
```

where *L* is mean length (&micro;m) and *W* is individual weight (&micro;g). Regressions return: carbon weight (C), dry weight (DW), or wet weight (WW) depending what was available in the literature. Where no regression exists, published mean weights are used.

The resulting weights are used to derive two final weight columns:

**Carbon weight (`final_C_ug`)** for biomass:

| Priority | Source |
|----------|--------|
| 1 | L-W regression returning C directly |
| 2 | L-W DW × (taxa-specific %C/DW) |
| 3 | Published mean C weight |
| 4 | Published DW × (%C/DW) |
| 5 | Published WW × 0.20 × (%C/DW) |
| 6 | L-W WW × 0.20 × (%C/DW) |

**Wet weight (`final_WW_ug`)** for ESD:

| Priority | Source |
|----------|--------|
| 1 | L-W regression returning WW directly |
| 2 | L-W DW ÷ 0.20 |
| 3 | L-W C ÷ (%C/DW) ÷ 0.20 |
| 4 | Published mean WW |
| 5 | Published DW ÷ 0.20 |
| 6 | Published C ÷ (%C/DW) ÷ 0.20 |

### Weight conversions

- **WW -> DW:** multiply by 0.20 (DW = 20% of WW)
- **DW -> C:** multiply by taxa-specific %C/DW where available; fallback = 40%

### Equivalent Spherical Diameter (ESD)

Each organism is approximated as a sphere. ESD is calculated from weight via:

```
WW  ->  volume  ->  ESD 
WW (µg) -> WW (g) -> volume (cm³) [÷ 1.05 g cm⁻³] -> ESD (cm) [(6V/π)^(1/3)] -> ESD (µm)
```
Zooplankton density is assumed to be 1.05 g cm⁻³ (slightly denser than seawater).


### NBSS
Organisms are binned into log<sub>2</sub> (octave) ESD bins. Biomass per bin is normalized by bin width (mg C m<sup>-2</sup> &micro;m<sup>-1</sup>). The NBSS slope is estimated by fitting a log10-log10 linear regression per station:

```
log10(B_norm) ~ log10(ESD_midpoint)
```

---

## Notes

- Net tows use 335 &micro;m mesh. Organisms < ~200 &micro;m are undersampled. Drop-off at the small end of the size spectrum is expected.

---

## Dependencies

```r
install.packages(c("here", "tidyverse", "broom", "scales", "ggrepel"))
```

R version used: 
<!-- Add R ver here -->

This project uses `renv` for package version management. To restore the exact environment:

```r
install.packages("renv")
renv::restore()
```

---

## Citation

If you use this code, please cite:
<!-- Add citation here -->
> Cabanelas, A. (*year*). NES-LTER Zooplankton Size Spectra pipeline. GitHub. [URL]

---

## License

<!-- Add license here -->
