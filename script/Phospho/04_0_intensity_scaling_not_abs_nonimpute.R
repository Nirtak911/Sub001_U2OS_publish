#load libraries
library(tidyverse)
library(protti)

## Load data ---------------------------
input <- file.path(data_in, "/02_loess_data_log2.csv")
meta <- file.path(data_in, "01_meta_data.csv")
genes <- file.path(data_in, "01_genes.csv")
data <- read_csv(input)
meta <- read_csv(meta)
genes <- read_csv(genes)
View(data)
## output data  ---------------------------
scale_impute_output <- file.path(data_in, "/04_0_scale_not_abs_data_log2.csv")
scale_impute_matrix_output <- file.path(data_in, "/04_0_scale_not_abs_data_log2_matrix.csv")


## scale and count function ------------------------------------------
scale_and_count <- function(data) {
  
  #scale intensity over fractions
  data_scale <- data %>%
    group_by(cell_line, treatment, r_replicate, ptm_collapse_key) %>%
    mutate(total_intensity = sum(original_intensity, na.rm = TRUE)) %>% 
    mutate(scaled_intensity = ifelse(total_intensity == 0, 0, original_intensity/total_intensity)) %>%
    ungroup()
  View(data_scale)
  
  # count in how many fractions a protein was identified
  fraction_count <- data_scale %>%
    group_by(cell_line, treatment, r_replicate, ptm_collapse_key) %>%
    summarise(fraction_count = n_distinct(fraction)) %>%
    ungroup()  
  
  # calculate percentage of proteins identified in all fractions/ at least 2 fractions  
  count_identified <- fraction_count  %>% 
    group_by(cell_line, treatment, r_replicate) %>%
    summarise(count_identified_all = sum(fraction_count >= 6), 
              count_identified_2 = sum(fraction_count >= 2),
              count_identified_total = sum(fraction_count >= 1)) %>%
    ungroup() %>% 
    mutate(pct_all = count_identified_all/count_identified_total,
           pct_2 = count_identified_2/count_identified_total)
  View(count_identified)
  
  return(list(data_scale = data_scale, 
              count_identified = count_identified))
}

## scale imputed intensities ------------------------------------------
#convert convert log2 imputed intensity back to original intensity
data_impute <- data %>%
  mutate(original_intensity = 2^(normalised_intensity_log2))

#scale data and count
data_impute <- data_impute %>% scale_and_count()
data_impute_scale <- data_impute$data_scale
count_identified_impute <- data_impute$count_identified

# Calulate the percentage of changing proteins
data_impute_scale <- data_impute_scale %>% 
  group_by(cell_line, r_replicate, fraction, ptm_collapse_key) %>% 
  mutate(
    CTRL_intensity = scaled_intensity[treatment == "CTRL"],
    pct_changing = scaled_intensity - CTRL_intensity) %>%
  ungroup()
View(data_impute_scale)

## bring data to matrix format---------------------------------------
# create matrix
data_scale_matrix <- data_impute_scale %>% 
  dplyr::select(sample, ptm_collapse_key, scaled_intensity) %>%
  pivot_wider(names_from = sample, 
              values_from = scaled_intensity)
View(data_scale_matrix)  

# convert NA to 0
data_scale_matrix_NA <- data_scale_matrix %>% 
  replace(is.na(.), 0)
View(data_scale_matrix_NA)

data_scale_NA <- data_scale_matrix_NA %>%
  pivot_longer(cols = -ptm_collapse_key, 
               names_to = "sample", 
               values_to = "scaled_intensity") %>%
  left_join(meta, by = "sample") %>%
  group_by(cell_line, r_replicate, fraction, ptm_collapse_key) %>% 
  mutate(
    CTRL_intensity = scaled_intensity[treatment == "CTRL"],
    pct_changing = scaled_intensity - CTRL_intensity) %>%
  ungroup() %>%
  left_join(genes, by = "ptm_collapse_key")
View(data_scale_NA)


## Save data ---------------------------
write_csv(data_scale_NA, scale_impute_output)
write_csv(data_scale_matrix_NA, scale_impute_matrix_output)

