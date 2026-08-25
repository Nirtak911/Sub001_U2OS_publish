# -------------------------------------------------------------------------
# SETUP: Load Libraries & Define Paths
# -------------------------------------------------------------------------

# Load required R packages for data manipulation, visualization, and analysis.
library(tidyverse)      # Core package for data manipulation and visualization (includes dplyr, ggplot2, etc.)
library(qvalue)

# --- Define File Paths ---

# Input data paths
diff_input <- file.path(data_in, "05_2_diff_mobility_limma_interaction.csv")
norm_input <- file.path(data_in, "04_scale_not_abs_loess_data_log2.csv")
genes_input <- file.path(data_in, "01_genes.csv")

# Output data paths
combined_output <- file.path(data_in, "06_0_diff_and_scaled_mobility_per_rep_final.csv")
median_output <- file.path(data_in, "06_0_diff_and_scaled_mobility_top_movements_median.csv")

# -------------------------------------------------------------------------
# DATA LOADING & PREPARATION
# -------------------------------------------------------------------------

# Load the datasets from the specified CSV files.
diff_data <- read_csv(diff_input)
norm_data_full <- read_csv(norm_input)
genes <- read_csv(genes_input)
View(norm_data_full)
View(diff_data)

# calculate avarage scaled intensity and pct_changing
norm_data <- norm_data_full %>% 
  filter(treatment != "CTRL")

# Assuming 'norm_data_avg' is your original data frame
movement_fraction_data <- norm_data %>%
  
  # 1. Group the data by both grouping variables
  group_by(pg_protein_groups, treatment, r_replicate) %>%
  
  # 2. Find the fraction ID AND the corresponding value for the extreme changes
  summarise(
    # --- Value Calculations ---
    # Calculate the actual highest positive value
    highest_positive_pct_change = max(pct_changing * (pct_changing > 0), na.rm = TRUE),
    
    # Calculate the actual lowest negative value
    # We use min(..., na.rm=TRUE) on the original vector, filtered for < 0
    lowest_negative_pct_change = min(pct_changing[pct_changing < 0], na.rm = TRUE),
    
    # --- Fraction ID Calculations ---
    # Find the fraction ID where avg_pct_changing is the highest positive
    highest_positive_fraction = fraction[which.max(pct_changing * (pct_changing > 0))],
    
    # Find the fraction ID where avg_pct_changing is the lowest negative
    lowest_negative_fraction = fraction[which.min(pct_changing * (pct_changing < 0))],
    
    # Calculate the total movement for filtering/context
    total_movement = sum(abs(pct_changing), na.rm = TRUE),
    
    # count the number of fractions 
    count_frac_changing = n_distinct(
      fraction[
        # Check if intensity is NOT zero AND is NOT missing
        abs(scaled_intensity) != 0 & !is.na(scaled_intensity)
      ]
    ),
    mobility_score = sum(abs(pct_changing), na.rm = TRUE)/count_frac_changing,
    
    .groups = 'drop' # Ungroup the data immediately
  ) %>%
  
  # 3. Create the 'movement' string
  mutate(
    dominant_movement = paste0(lowest_negative_fraction, "->", highest_positive_fraction)
  ) 

# View the final summary table
View(movement_fraction_data)

# -------------------------------------------------------------------------
# 1. Start with the Classified Data
classified_data <- norm_data %>%
  filter(abs(pct_changing) > 0) %>%
  mutate(
    change_type = case_when(
      pct_changing < 0 ~ "Loss",
      pct_changing > 0 ~ "Gain",
      TRUE ~ "None"
    ),
    pct_changing = pct_changing
  )

# 2. Split the Data into Loss and Gain Tables
loss_data <- classified_data %>%
  filter(change_type == "Loss") %>%
  dplyr::select(pg_protein_groups, treatment, r_replicate,
                loss_fraction = fraction,
                loss_pct_change = pct_changing)

gain_data <- classified_data %>%
  filter(change_type == "Gain") %>%
  dplyr::select(pg_protein_groups, treatment, r_replicate,
                gain_fraction = fraction,
                gain_pct_change = pct_changing)

# 3. Perform the Cartesian Join (Cross-Join) by Group
# This is an EQUI-JOIN on the grouping keys, which forces the Cartesian product
# to happen ONLY within the same protein/treatment/time group, but it does it
# in a single, highly optimized operation, replacing the slow 'do()' loop.

movement_combinations_data <- loss_data %>%
  inner_join(
    gain_data,
    by = c("pg_protein_groups", "treatment", "r_replicate"),
    relationship = "many-to-many" # Explicitly allows the Cartesian product nature
  ) %>%
  
  # 4. Calculate Final Metrics (as before)
  mutate(
    movement = paste0(loss_fraction, "->", gain_fraction),
    pair_movement_score = (abs(loss_pct_change) + abs(gain_pct_change)) / 2
  ) %>%
  
  # 5. Add Contextual Metadata (optional)
  left_join(
    movement_fraction_data %>% dplyr::select(pg_protein_groups, treatment, r_replicate, 
                                             dominant_movement, highest_positive_fraction, lowest_negative_fraction, 
                                             total_movement, count_frac_changing, mobility_score),
    by = c("pg_protein_groups", "treatment", "r_replicate")
  )

View(movement_combinations_data)

# -------------------------------------------------------------------------
# add movement_fraction_data to diff_data
combined_data <- movement_combinations_data %>%
  left_join(diff_data %>% 
              dplyr::select(pg_protein_groups, Treatment, comparison,
                            diff, pval, adj_pval, movement, classified_movement) %>%
              rename("treatment" = "Treatment"),
            by = c("pg_protein_groups", "treatment", "movement")) 
View(combined_data)    

combined_data_filtered <- combined_data %>%
  filter(!is.na(comparison))
View(combined_data_filtered)

# select movements that are the same in three reps
movements_in_all_reps <- combined_data_filtered %>%
  distinct(pg_protein_groups, treatment, r_replicate, movement) %>% 
  group_by(pg_protein_groups, treatment, movement) %>%
  summarise(n_reps_with_movement = n(), .groups = "drop") %>%
  group_by(pg_protein_groups, treatment) %>%
  summarise(
    movements_present_in_all = list(movement[n_reps_with_movement >= 3]),
    .groups = "drop"
  ) %>%
  filter(lengths(movements_present_in_all) > 0)
View(movements_in_all_reps)

movements_in_all_reps_expanded <- movements_in_all_reps %>%
  tidyr::unnest(movements_present_in_all)
View(movements_in_all_reps_expanded)

# filter combined_data to only include movements present in all reps
combined_data_final <- combined_data_filtered %>%
  inner_join(movements_in_all_reps_expanded %>%
               rename(movement = movements_present_in_all),
             by = c("pg_protein_groups", "treatment", "movement"))
View(combined_data_final)

# Count how many paths the script gets to choose from
movement_counts <- combined_data_final %>%
  group_by(pg_protein_groups, treatment) %>%
  summarise(k_possible_movements = n_distinct(movement), .groups = "drop")

# median values per top movements
combined_data_median <- combined_data_final %>%
  group_by(pg_protein_groups, treatment, movement, classified_movement) %>%
  summarise(
    diff = median(diff, na.rm = TRUE),
    pval = median(pval, na.rm = TRUE),
    adj_pval = median(adj_pval, na.rm = TRUE),
    median_total_movement = median(total_movement, na.rm = TRUE),
    median_mobility_score = median(mobility_score, na.rm = TRUE),
    median_count_frac_changing = median(count_frac_changing, na.rm = TRUE),
    .groups = "drop"
  )
View(combined_data_median)
 
# -------------------------------------------------------------------------
# SAVE RESULTS
# -------------------------------------------------------------------------

# Save the protein-to-cluster assignments to a CSV file.
write.csv(combined_data_final, combined_output, row.names = FALSE)
write.csv(combined_data_median, median_output, row.names = FALSE)
