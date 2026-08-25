#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)

## Load data ---------------------------------------------------------------
input <- file.path(data_in, "02_loess_data_log2.csv")
meta <- file.path(data_in, "01_meta_data.csv")

data_input <- read_csv(input, col_types = cols(...1 = col_skip()))
meta_input <- read_csv(meta)

View(data_input)
View(meta_input)

## output data -------------------------------------------------------------
impute_output <- file.path(data_in, "/03_0_missingness_protti_data_log2.csv")
impute_matrix_output <- file.path(data_in, "/03_0_missingness_protti_data_log2_matrix.csv")
plot_impute_output <- file.path(data_in, "/Plots/03_0_missingness_protti_data_log2_plot.pdf")

## Split data frame by Fractions and cell line------------------------------
FR_list_norm <- data_input %>% group_split(cell_line, fraction)

## Imputation -------------------------------------------------------------
#identify missing values
# Create an empty list to store the results
FR_list_missing <- list()

# Loop through each data frame in the list
for (df in FR_list_norm) {
  df <- df %>% mutate(ref_condition = ifelse(str_detect(r_condition, "CTRL"), r_condition, NA))
  
  # Extract the first reference condition that contains "CTRL"
  ref_condition <- df$ref_condition[!is.na(df$ref_condition)][1]
  
  data_missing <- df %>%
    assign_missingness(sample = sample,
                       condition = r_condition,
                       grouping = pg_protein_groups,
                       intensity = normalised_intensity_log2,
                       ref_condition = ref_condition,
                       completeness_MAR = 0.8, # means that if 1 out of 3 is missing, it will be imputed
                       completeness_MNAR = 0.3, # if there is 1 values in one cond and the other cond is complete, it will be imputed
                       retain_columns = c(pg_genes, 
                                          pg_quantity,
                                          log2_pg_quantity,
                                          r_condition,
                                          r_replicate,
                                          r_file_name,
                                          cell_line,
                                          treatment,
                                          fraction))
  
  # Store the result in the result list
  FR_list_missing[[length(FR_list_missing) + 1]] <- data_missing
}
View(FR_list_missing[[1]])

# Create an empty list to store the results
FR_list_impute <- list()

# Loop through each data frame in the list
for (df in FR_list_missing) {
  data_impute <- df %>% 
    protti::impute(sample = sample,
                   condition = r_condition,
                   grouping = pg_protein_groups,
                   intensity = normalised_intensity_log2,
                   comparison = comparison,
                   missingness = missingness,
                   method = "ludovic",
                   skip_log2_transform_error = TRUE,
                   retain_columns = c(
                                      pg_genes, 
                                      pg_quantity,
                                      log2_pg_quantity,
                                      r_condition,
                                      r_replicate,
                                      r_file_name,
                                      cell_line,
                                      treatment,
                                      fraction))
  
  # Store the result in the result list
  FR_list_impute[[length(FR_list_impute) + 1]] <- data_impute
}
print(FR_list_impute)
View(FR_list_impute[[1]])

# display intensity distribution after imputation
# Create an empty list to store plots
FR_list_impute_plots <- list()

for (df in FR_list_impute) {
  p1_imputed_distribution <- qc_intensity_distribution(
    data = df,
    sample = imputed,
    grouping = pg_protein_groups,
    intensity_log2 = imputed_intensity,
    plot_style = "histogram"
  )
  
  p2_imputed_distribution <- qc_intensity_distribution(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity_log2 = imputed_intensity,
    plot_style = "histogram"
  )
  
  # Append the plots to the list
  FR_list_impute_plots[[length(FR_list_impute_plots) + 1]] <- p1_imputed_distribution
  FR_list_impute_plots[[length(FR_list_impute_plots) + 1]] <- p2_imputed_distribution
}
FR_list_impute_plots

## Merge data frames into one--------------
# combine the imputed data list
data_impute <- bind_rows(FR_list_impute) 

data_impute$imputed_intensity_log2 <- data_impute$imputed_intensity
data_impute$imputed_intensity <- 2^(data_impute$imputed_intensity_log2)
View(data_impute) 

missingness_counts <- data_impute %>%
  # Keep only one row for each unique combination of the case identifiers
  # The 'missingness' value for this unique combination will be preserved
  distinct(pg_protein_groups, comparison, missingness) %>%
  
  # Count the frequency of each 'missingness' value
  count(missingness, name = "count")

# Print the resulting counts
print(missingness_counts)

# create matrix
data_impute_matrix <- data_impute %>% 
  dplyr::select(sample, pg_protein_groups, imputed_intensity_log2) %>%
  pivot_wider(names_from = sample, 
              values_from = imputed_intensity_log2,
              values_fn = function(x) mean(x, na.rm = TRUE))

View(data_impute_matrix)

p3_imputed_distribution <- qc_intensity_distribution(
  data = data_impute,
  sample = sample,
  grouping = pg_protein_groups,
  intensity_log2 = imputed_intensity_log2,
  plot_style = "histogram")
p3_imputed_distribution

## Save data ---------------------------
write.csv(data_impute, impute_output)
write.csv(data_impute_matrix, impute_matrix_output)

## Save plots --------------------------
pdf(file = plot_impute_output, width = 18, height = 10)
# Loop through each plot in the list
for(i in seq_along(FR_list_impute_plots)) {
  # Print the plot to the PDF
  print(FR_list_impute_plots[[i]])
}
p3_imputed_distribution
dev.off()


