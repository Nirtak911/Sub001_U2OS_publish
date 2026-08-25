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
View(data)
View(genes)

# ## output data  ---------------------------
output <- file.path(data_in,"/05_2_diff_mobility_limma_interaction.csv")
output_translocation <- file.path(data_in,"/05_2_translocation_summary_pval.csv")
facet_output <- file.path(data_in,"/Plots/05_2_diff_mobility_limma_interaction_facet_plot.pdf")

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
View(data_movement)

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
View(data_movement)

# plot pval distribution
ggplot(data_movement, aes(x = pval)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
  theme_minimal() +
  labs(
    title = "Distribution of P-Values",
    x = "P-Value",
    y = "Frequency"
  )

# Calculate the standard deviation of all pairwise differences
sd_diffs <- sd(data_movement$diff, na.rm = TRUE)
mean_diffs <- mean(data_movement$diff, na.rm = TRUE)

# Calculate the cutoffs
two_sd_cutoff <- 2 * sd_diffs
three_sd_cutoff <- 3 * sd_diffs

# Print the results
print(paste("Standard Deviation:", sd_diffs))
print(paste("2xSD Cutoff:", two_sd_cutoff))
print(paste("3xSD Cutoff:", three_sd_cutoff))

# Create the histogram and add the cutoff lines
ggplot(data_movement, aes(x = diff)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  
  # --- Add 2 SD Cutoffs (Mean +/- 2*SD) ---
  # Negative side (Mean - 2*SD)
  geom_vline(xintercept = mean_diffs - two_sd_cutoff, 
             linetype = "solid", 
             color = "red", 
             size = 0.8) +
  # Positive side (Mean + 2*SD)
  geom_vline(xintercept = mean_diffs + two_sd_cutoff, 
             linetype = "solid", 
             color = "red", 
             size = 0.8) +
  
  # --- Add 3 SD Cutoffs (Mean +/- 3*SD) ---
  # Negative side (Mean - 3*SD)
  geom_vline(xintercept = mean_diffs - three_sd_cutoff, 
             linetype = "dashed", 
             color = "blue", 
             size = 0.8) +
  # Positive side (Mean + 3*SD)
  geom_vline(xintercept = mean_diffs + three_sd_cutoff, 
             linetype = "dashed", 
             color = "blue", 
             size = 0.8) +
  
  theme_minimal() +
  labs(
    title = "Distribution of Final Scores",
    x = "Diff of logFC",
    y = "Frequency"
  )


# Label significant------------------------------------------------------------
# Create a new column for significant
data_translocation_classified <- data_movement %>%
  mutate(
    # Create the significance category using ordered tiers (case_when)
    significance = case_when(
      
      # TIER I: Highest Confidence (FDR-controlled + High Effect)
      # adj_pval <= 0.05 & abs(diff) >= three_sd_cutoff
      adj_pval <= 0.05 & abs(diff) >= three_sd_cutoff ~ "Class_I_FDR_Sig_High_FC",
      
      # TIER II (NEW): FDR-controlled Significant (FDR-controlled + Medium Effect)
      # This addresses your request, assuming medium effect (two_sd_cutoff)
      # adj_pval <= 0.05 & abs(diff) >= two_sd_cutoff 
      adj_pval <= 0.1 & abs(diff) >= two_sd_cutoff ~ "Class_II_FDR_Sig_Medium_FC",
      
      # TIER III (Was II): Highly Significant (Strict P-value + High Effect)
      # pval <= 0.01 & abs(diff) >= three_sd_cutoff
      pval <= 0.01 & abs(diff) >= three_sd_cutoff ~ "Class_III_Highly_Sig_High_FC",
      
      # TIER IV (Was III): Highly Significant (Strict P-value + Medium Effect)
      # pval <= 0.01 & abs(diff) >= two_sd_cutoff
      pval <= 0.01 & abs(diff) >= two_sd_cutoff ~ "Class_IV_Highly_Sig_Medium_FC",
      
      # # TIER V (Was IV): Significant (Moderate P-value + Medium Effect)
      # # pval <= 0.05 & abs(diff) >= two_sd_cutoff
      # pval <= 0.05 & abs(diff) >= two_sd_cutoff ~ "Class_V_Significant_Medium_FC",
      
      # Default/Catch-all: If none of the above thresholds are met
      TRUE ~ "Not_Significant"
    )
  )

data_translocation_classified <- data_translocation_classified %>%
  mutate(
    # 1. Create the simple binary classification column
    significance_plot_group = ifelse(
      # If the tier is NOT 'Not_Significant', classify it as 'Significant'
      significance != "Not_Significant",
      "Significant",
      "Not_Significant"
    )
  ) %>%
  mutate(
    # 2. Update the color group column to use the correct tier column
    plot_color_group = ifelse(
      # Only assign movement color if it passed any significance tier
      significance != "Not_Significant",
      as.character(movement),
      "Not_Significant"
    )
  )

View(data_translocation_classified)

# Count the number of translocations in each significance tier
all_class_counts <- data_translocation_classified %>%
  # Group by the newly created tier column
  group_by(significance) %>%
  # Count the number of rows (translocations) in each group
  summarise(Count = n())

# Print the result
print("Counts by Significance Tier:")
print(all_class_counts)

unique_pairs_per_tier <- data_translocation_classified %>%
  # 1. Group the data by the significance tier
  group_by(significance) %>%
  
  # 2. Summarize the count of unique pairs within each tier
  # We use n_distinct() to count the unique combinations of the three columns
  summarise(
    Num_Unique_Pairs = n_distinct(gene_PTM_key, Treatment, classified_movement),
    .groups = 'drop' # Ensure ungrouping after summary
  ) %>%
  # Optional: Order by the number of unique pairs
  arrange(desc(Num_Unique_Pairs))

# Print the result
print("Number of unique (Gene, Treatment, Movement) pairs per Significance Tier:")
print(unique_pairs_per_tier)

# Filter out general abundance changes
#Isolate the relevant Omnibus results
omnibus_results <- data %>%
  # Filter only for the Omnibus contrasts
  filter(Contrast_Type == "Omnibus_Treatment_Effect") %>%
  # Select only the relevant columns for joining and filtering
  dplyr::select(
    ptm_collapse_key, 
    Treatment,
    omnibus_adj_pval = adj_pval,  # Rename for clarity
    omnibus_diff = diff
  )

# Join the movement data with the Omnibus results
data_movement_omnibus <- data_translocation_classified %>%
  # Join by the unique protein identifier and the treatment name
  left_join(
    omnibus_results, 
    by = c("ptm_collapse_key", "Treatment")
  )
View(data_movement_omnibus)

## Save data ---------------------------
write_csv(data_movement_omnibus, output)
