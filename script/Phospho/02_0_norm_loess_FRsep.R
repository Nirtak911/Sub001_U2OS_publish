#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)
library(limma)

## Load data ---------------------------------------------------------------
input <- file.path(data_in, "01_data_log2.csv")
meta <- file.path(data_in, "01_meta_data.csv")
genes <- file.path(data_in, "01_genes.csv")
data <- read_csv(input)
meta <- read_csv(meta)
genes <- read_csv(genes)
View(data)
View(meta)
View(genes)

## output data -------------------------------------------------------------
norm_output <- file.path(data_in, "/02_loess_data_log2.csv")
norm_matrix_output <- file.path(data_in, "/02_loess_data_log2_matrix.csv")
plot_norm_output <- file.path(data_in, "/Plots/02_loess_data_log2_plot.pdf")

## Valid value filtering----------------------------------------------------
#convert intensities <1 to NA
df <- data %>%
  mutate(log2_ptm_quantity = ifelse(log2_ptm_quantity < 0, NA, log2_ptm_quantity))

## Split data frame by Fractions and lysis------------------------------
FR_list <- df %>% group_split(fraction)
print(FR_list)

# Function at least identified in 3 out of 4 reps per condition
msk_rows_function <- function(df) {
  msk_rows <- df %>%
    group_by(r_condition, ptm_collapse_key) %>%
    summarise(na_sum = sum(!is.na(ptm_quantity)) >= 3) %>%
    group_by(ptm_collapse_key) %>%
    summarise(cond_sum = sum(na_sum)) %>%
    filter(cond_sum >= 1) %>%
    ungroup()
  
  # Filter df based on the pg_genes in msk_rows
  df_vv <- df %>%
    filter(ptm_collapse_key %in% msk_rows$ptm_collapse_key)
  
  # Return the filtered data frame
  return(df_vv)
}

FR_list_vv <- lapply(FR_list, msk_rows_function)
print(FR_list_vv)

# combine lists
data_vv <- bind_rows(FR_list_vv)

# data_filter matrix
data_vv_matrix <- data_vv %>% 
  dplyr::select(sample, ptm_collapse_key, log2_ptm_quantity) %>%
  pivot_wider(names_from = sample, 
              values_from = log2_ptm_quantity)
View(data_vv_matrix)

# split data_vv by fraction
FR_list_vv <- data_vv %>% group_split(fraction)
print(FR_list_vv)

## Normalisation -------------------------------------------------------------
#normalize data
# Create an empty list to store the results
FR_list_norm <- list()

# Loop through each data frame in the list
for (df in FR_list_vv) {
  # df to wide format
  df_wide <- df %>%
    pivot_wider(names_from = sample, values_from = log2_ptm_quantity, id_cols = ptm_collapse_key) %>%
    column_to_rownames(var = "ptm_collapse_key")
  
  norm_df <- normalizeCyclicLoess(df_wide)
  
  # Convert the normalised data frame back to long format
  df_long <- as.data.frame(norm_df) %>%
    rownames_to_column("ptm_collapse_key") %>%
    pivot_longer(cols = -ptm_collapse_key, names_to = "sample", values_to = "normalised_intensity_log2")
  
  # add back all info from df
  df_long <- df %>%
    left_join(df_long, by = c("ptm_collapse_key", "sample"))
  

  # Store the result in the result list
  FR_list_norm[[length(FR_list_norm) + 1]] <- df_long
}

# run intensities distribution normalised
# Create an empty list to store plots
FR_list_norm_plots <- list()

for (df in FR_list_norm) {
p1_Intensity <- qc_intensity_distribution(
                      data = df,
                      sample = sample,
                      grouping = ptm_collapse_key,
                      intensity_log2 = log2_ptm_quantity,
                      plot_style = "boxplot"
                    )
                    p1_Intensity

p1_Intensity_norm <- qc_intensity_distribution(
                        data = df,
                        sample = sample,
                        grouping = ptm_collapse_key,
                        intensity_log2 = normalised_intensity_log2,
                        plot_style = "boxplot"
                      )
# Append the plots to the list
FR_list_norm_plots[[length(FR_list_norm_plots) + 1]] <- p1_Intensity
FR_list_norm_plots[[length(FR_list_norm_plots) + 1]] <- p1_Intensity_norm
}
FR_list_norm_plots

print(FR_list_norm)

## Merge data frames into one--------------
norm_df <- bind_rows(FR_list_norm)

# create matrix
norm_df_matrix <- norm_df %>% 
  dplyr::select(sample, ptm_collapse_key, normalised_intensity_log2) %>%
  pivot_wider(names_from = sample, 
              values_from = normalised_intensity_log2)
View(norm_df_matrix)

## Save data ---------------------------
write_csv(norm_df, norm_output)
write_csv(norm_df_matrix, norm_matrix_output)


## Save plots --------------------------
pdf(file = plot_norm_output, width = 18, height = 10)
# Loop through each plot in the list
for(i in seq_along(FR_list_norm_plots)) {
  # Print the plot to the PDF
  print(FR_list_norm_plots[[i]])
}
dev.off()

