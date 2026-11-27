# paths (data_in, raw_data, annot_dir); run from the project root
source("script/Proteome/00_path_definitions.R")

#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(patchwork)


## Load data ---------------------------
input <- file.path(data_in,"05_0_diff_data_limma_dapar_impute_complete_group_omnibus.csv")
data <- read_csv(input)
input2 <- file.path(data_in,"01_genes.csv")
genes <- read_csv(input2)
if (interactive()) View(data)
if (interactive()) View(genes)

# ## output data  ---------------------------
output <- file.path(data_in,"/05_2_diff_mobility_limma_interaction.csv")

# define colors for the annotations
annotation_colors= c("Frac1" = "#E7872B", "Frac2" = "#F2B341", 
                     "Frac3" = "#F0E442", "Frac4" = "#009E73", 
                     "Frac5" = "#56B4E9", "Frac6" = "#0072B2")

# add time column to data
data <- data %>%
  dplyr::mutate(
    time = stringr::str_extract(comparison, "(?<=_in_).+$")
  )

# filter for interactions
data_movement <- data %>% filter(Contrast_Type == "Differential_Movement") %>%
  # Apply the confirmed logic for movement direction
  dplyr::mutate(
    movement = dplyr::if_else(
      # If diff (logFC) is POSITIVE, movement is TOWARDS Fraction A
      diff > 0,
      paste0(Fraction_B, "->", Fraction_A),
      
      # If diff (logFC) is NEGATIVE, movement is TOWARDS Fraction B
      paste0(Fraction_A, "->", Fraction_B),
      
      # Handle cases where diff is exactly 0 (no significant movement)
      "No Movement"
    ))

# add classifier for the movement
data_movement <- data_movement %>%
  mutate(
    start_frac = as.numeric(str_extract(movement, "(?<=Frac)\\d(?=->)")),
    end_frac = as.numeric(str_extract(movement, "(?<=->Frac)\\d"))
  )
if (interactive()) View(data_movement)

# Add the new movement classification column with no-movement categories
data_movement <- data_movement %>%
  mutate(
    classified_movement = case_when(
      # Inter-fractional movements
      (start_frac %in% c(1, 2) & end_frac %in% c(3, 4)) ~ "cytosol->membrane",
      (start_frac %in% c(3, 4) & end_frac %in% c(1, 2)) ~ "membrane->cytosol",
      (start_frac %in% c(1, 2) & end_frac %in% c(5, 6)) ~ "cytosol->nuclear",
      (start_frac %in% c(5, 6) & end_frac %in% c(1, 2)) ~ "nuclear->cytosol",
      (start_frac %in% c(3, 4) & end_frac %in% c(5, 6)) ~ "membrane->nuclear",
      (start_frac %in% c(5, 6) & end_frac %in% c(3, 4)) ~ "nuclear->membrane",
      
      # Intra-fractional "no-movement" classifications
      (start_frac %in% c(1, 2) & end_frac %in% c(1, 2)) ~ "cytosol",
      (start_frac %in% c(3, 4) & end_frac %in% c(3, 4)) ~ "membrane",
      (start_frac %in% c(5, 6) & end_frac %in% c(5, 6)) ~ "nuclear",
      
      TRUE ~ NA_character_ # Fills in NA for any other cases
    )
  )
if (interactive()) View(data_movement)

# Filter out general abundance changes
#Isolate the relevant Omnibus results
omnibus_results <- data %>%
  # Filter only for the Omnibus contrasts
  filter(Contrast_Type == "Omnibus_Treatment_Effect") %>%
  # Select only the relevant columns for joining and filtering
  dplyr::select(
    pg_protein_groups, 
    Treatment,
    omnibus_adj_pval = adj_pval,  # Rename for clarity
    omnibus_diff = diff
  )

# Join the movement data with the Omnibus results
data_movement_omnibus <- data_movement %>%
  # Join by the unique protein identifier and the treatment name
  left_join(
    omnibus_results, 
    by = c("pg_protein_groups", "Treatment")
  )
if (interactive()) View(data_movement_omnibus)

## Save data ---------------------------
write_csv(data_movement_omnibus, output)
