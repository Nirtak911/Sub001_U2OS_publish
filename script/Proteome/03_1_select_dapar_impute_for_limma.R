#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)

## Load data ---------------------------------------------------------------
input_norm <- file.path(data_in, "02_loess_data_log2.csv")
meta <- file.path(data_in, "01_meta_data.csv")
missing <- file.path(data_in, "/03_0_missingness_protti_data_log2.csv")
input_impute <- file.path(data_in, "/03_imputeFRsep_data_log2.csv")

data_norm_input <- read_csv(input_norm, col_types = cols(...1 = col_skip()))
meta_input <- read_csv(meta)
missing_input <- read_csv(missing, col_types = cols(...1 = col_skip()))
data_impute_input <- read_csv(input_impute, col_types = cols(...1 = col_skip()))


View(data_norm_input)
View(meta_input)
View(missing_input)
View(data_impute_input)


## output data -------------------------------------------------------------
limma_output <- file.path(data_in, "/03_1_dapar_limma_input_data_log2.csv")
limma_output_matrix <- file.path(data_in, "/03_1_dapar_limma_input_data_log2_matrix.csv")

# Step 1: Summarise missingness per (pg_protein_groups, sample) pair
missing_collected <- missing_input %>%
  group_by(pg_protein_groups, sample) %>%
  summarise(
    # Check if all values are identical
    are_all_same = n_distinct(missingness) == 1,
    
    # Check for presence of 'MNAR' and 'MAR' in the group
    has_complete = any(missingness == "complete"),
    has_mnar = any(missingness == "MNAR"),
    has_mar = any(missingness == "MAR"),
    
    # Apply the complex conditional rule
    final_missingness = case_when(
      # 1. SIMPLE CASE: If all values are identical, use that value.
      are_all_same ~ missingness[1],
      
      # 2. COMPLEX CASE (Different values): Check for MNAR/MAR priority.
      has_complete ~ "complete",
      #    Rule A: If different AND includes 'MNAR', assign 'MNAR'.
      has_mnar ~ "MNAR",
      
      #    Rule B: If different AND does NOT have 'MNAR' but DOES have 'MAR', assign 'MAR'.
      has_mar ~ "MAR",
      
      # 3. DEFAULT CASE: If none of the above are true, assign NA.
      TRUE ~ NA_character_
    ),
    .groups = 'drop'
  )
View(missing_collected)

# data_impute_protti_condensed <- missing_input %>%
#   # 1. Keep only the rows that have a valid imputed intensity (removes the NA row)
#   filter(!is.na(imputed_intensity_log2)) %>%
#   distinct(pg_protein_groups, sample, imputed_intensity_log2, .keep_all = TRUE) 
# View(data_impute_protti_condensed)

# join imputed data to the normalised data
data <- data_norm_input %>%
  left_join(data_impute_input %>% dplyr::select(pg_protein_groups, sample, imputed_intensity_log2),
            by = c("pg_protein_groups", "sample")) %>%
  left_join(missing_collected %>% dplyr::select(pg_protein_groups, sample, final_missingness)) %>% 
  mutate( # Create the final intensity column based on the missingness mechanism
    final_intensity_log2 = case_when(
      # 1. If final_missingness is "complete" OR "MNAR", use the imputed value
      final_missingness %in% c("MNAR") ~ imputed_intensity_log2,
      
      # 2. If final_missingness is "MAR" OR NA, use the original normalised intensity
      #    (The MAR/NA cases are often handled by a different imputation method,
      #     or kept as-is if no imputation was done for them yet.)
      final_missingness %in% c("complete","MAR", NA_character_) ~ normalised_intensity_log2,
      
      # 3. Default (If any other unexpected value exists, default to normalised_intensity)
      TRUE ~ normalised_intensity_log2
    )
  )
View(data)

# create matrix
data_matrix <- data %>%
  dplyr::select(pg_protein_groups, sample, final_intensity_log2) %>%
  pivot_wider(names_from = sample, values_from = final_intensity_log2) %>%
  arrange(pg_protein_groups) %>%
  column_to_rownames(var = "pg_protein_groups")
View(data_matrix)

# save as csv-------------------------------------------------------------
write.csv(data, limma_output)
write.csv(data_matrix, limma_output_matrix)

# End of script -----------------------------------------------------------