################################################################################
## Script:  diagnostic_plots/check_weight_coverage.R
## Project: NES-LTER Zooplankton Size Spectra (Pelagic Synthesis WG)
## Author:  Alexandra C. Cabanelas
## Created: August 2024  |  Updated: February 2026
##
## Purpose: Plot showing which weights are available per stage + taxa
##          Run after 03_weights.R.
## Inputs:  data/processed/03_zp_weights.rds
################################################################################

source(here::here("scripts", "00_packages.R"))

zp <- readRDS(here("data", "processed", "03_zp_weights.rds"))

## ------------------------------------------ ##
#      get coverage by weight type -----
## ------------------------------------------ ##
coverage <- zp %>%
  distinct(taxa_name, stage,
           C_weight_ug, DW_weight_ug, WW_weight_ug,
           pub_C_ug,    pub_DW_ug,    pub_WW_ug) %>%
  pivot_longer(
    cols = c(C_weight_ug, DW_weight_ug, WW_weight_ug,
             pub_C_ug,    pub_DW_ug,    pub_WW_ug),
    names_to  = "column",
    values_to = "value"
  ) %>%
  mutate(
    has_value = !is.na(value),
    taxa_stage = paste(taxa_name, stage, sep = " | "),
    column = factor(column, 
                    levels = c("C_weight_ug", "DW_weight_ug", "WW_weight_ug",
                               "pub_C_ug",    "pub_DW_ug",    "pub_WW_ug"))
  )

# --- function to plot staged and unstaged separately so its not too crowded ---
plot_coverage <- function(data, title = NULL) {
  taxa_order_local <- data %>%
    group_by(taxa_stage) %>%
    summarise(n_filled = sum(has_value), .groups = "drop") %>%
    arrange(desc(n_filled)) %>%
    pull(taxa_stage)
  
  data %>%
    mutate(
      taxa_stage = factor(taxa_stage, levels = rev(taxa_order_local)),
      source = case_when(
        column %in% c("C_weight_ug", "DW_weight_ug", "WW_weight_ug") ~ "L-W Regression",
        TRUE ~ "Published"
      ),
      weight_type = case_when(
        column %in% c("C_weight_ug",  "pub_C_ug")  ~ "C",
        column %in% c("DW_weight_ug", "pub_DW_ug") ~ "DW",
        column %in% c("WW_weight_ug", "pub_WW_ug") ~ "WW"
      ),
      source = factor(source, levels = c("L-W Regression", "Published"))
    ) %>%
    ggplot(aes(x = weight_type, y = taxa_stage, fill = has_value)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(has_value, "✓", "")),
              size = 3, color = "white") +
    scale_fill_manual(values = c("FALSE" = "#f0f0f0", "TRUE" = "#2e9e5b")) +
    scale_x_discrete(expand = c(0, 0)) +   
    scale_y_discrete(expand = c(0, 0)) +   
    facet_grid(. ~ source, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, fill = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.y      = element_text(size = 7),
      axis.text.x      = element_text(size = 9, face = "bold"),
      strip.text       = element_text(size = 9, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      legend.position  = "none",
      panel.grid       = element_blank()
    )
}

## ------------------------------------------ ##
#      Plots -----
## ------------------------------------------ ##
# --- 1) unstaged
coverage %>%
  filter(stage == "Not_Staged") %>%
  mutate(taxa_stage = taxa_name) %>%  
  plot_coverage()

# --- 2) staged
coverage %>%
  filter(stage != "Not_Staged") %>%
  mutate(taxa_stage = paste(taxa_name, stage, sep = " | ")) %>%
  arrange(stage, taxa_name) %>% 
  mutate(
    taxa_stage = factor(taxa_stage, levels = rev(unique(taxa_stage)))  
  ) %>%
  plot_coverage()

# --- 3) staged; cleaner with facets
coverage %>%
  filter(!stage %in% c("Not_Staged", "Cryptopia")) %>%
  mutate(
    source = case_when(
      column %in% c("C_weight_ug", "DW_weight_ug", "WW_weight_ug") ~ "L-W Regression",
      TRUE ~ "Published"
    ),
    weight_type = case_when(
      column %in% c("C_weight_ug",  "pub_C_ug")  ~ "C",
      column %in% c("DW_weight_ug", "pub_DW_ug") ~ "DW",
      column %in% c("WW_weight_ug", "pub_WW_ug") ~ "WW"
    ),
    source      = factor(source, levels = c("L-W Regression", "Published")),
    weight_type = factor(weight_type, levels = c("C", "DW", "WW"))
  ) %>%
  ggplot(aes(x = weight_type, y = taxa_name, fill = has_value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(has_value, "✓", "")), size = 2.5, color = "white") +
  scale_fill_manual(values = c("FALSE" = "#f0f0f0", "TRUE" = "#2e9e5b")) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(stage ~ source, scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, fill = NULL,
       title = "Weight data coverage by taxa and stage") +
  theme_bw() +
  theme(
    axis.text.y      = element_text(size = 7),
    axis.text.x      = element_text(size = 9, face = "bold"),
    strip.text       = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    legend.position  = "none",
    panel.grid       = element_blank()
  )
