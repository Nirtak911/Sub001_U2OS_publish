#load libraries
library(tidyverse)
library(DAPAR)
library(MSnbase)
library(imp4p)

## Load data ---------------------------------------------------------------
input <- file.path(data_in, "/02_loess_data_log2.csv")
meta <- file.path(data_in, "01_meta_data.csv")

data_input <- read_csv(input)
meta_input <- read_csv(meta)

View(data_input)
View(meta_input)

## output data -------------------------------------------------------------
impute_output <- file.path(data_in, "/03_imputeFRsep_data_log2.csv")
impute_matrix_output <- file.path(data_in, "/03_imputeFRsep_data_log2_matrix.csv")
plot_impute_output <- file.path(data_in, "/Plots/03_imputeFRsep_data_log2_plot.pdf")

## Split data frame by Fractions and cell line------------------------------
FR_list_norm <- data_input %>% group_split(fraction, cell_line)
print(FR_list_norm)

## Imputation -------------------------------------------------------------
# Create an empty list to store the results
FR_list_impute <- list()

# Loop through each data frame in the list
for (df in FR_list_norm) {
  
  ## Process data -----------------------------------------------------------
  # create matrix
  norm_df_matrix <- df %>% 
    dplyr::select(sample, ptm_collapse_key, normalised_intensity_log2) %>%
    pivot_wider(names_from = sample, 
                values_from = normalised_intensity_log2)
  print(norm_df_matrix)
  
  # extract samples from data
  samples <- norm_df_matrix %>%
    pivot_longer(cols = -ptm_collapse_key, names_to = "sample", values_to = "value") %>%
    dplyr::select(sample) %>% 
    distinct()
  print(samples)
  
  data <- norm_df_matrix %>% 
    column_to_rownames("ptm_collapse_key")
  
  meta <- samples %>%
    left_join(meta_input, by = "sample") %>%
    dplyr::select(sample, r_condition, r_replicate) %>%
    dplyr::rename(Condition = r_condition, Replicate = r_replicate)
  
  # Create MSnSet
  obj <- MSnSet(exprs = as.matrix(data),
                pData = meta %>%
                  mutate(Sample.name = sample) %>% 
                  column_to_rownames("sample"),
                fData = data)
  
  # Copied from the createMSnset function to have the correct format
  pep_prot_data <- "protein"
  obj@experimentData@other$typeOfData <- pep_prot_data
  obj@experimentData@other$Prostar_Version <- NA
  obj@experimentData@other$proteinId <- NULL
  obj@experimentData@other$keyId <- "PG.UniProt"
  obj@experimentData@other$RawPValues <- FALSE
  
  metacell <- NULL
  metacell <- BuildMetaCell(from = "maxquant", level = "protein", 
                            qdata = Biobase::exprs(obj), conds = Biobase::pData(obj)$Condition, 
                            df = metacell)
  Biobase::fData(obj) <- cbind(Biobase::fData(obj), metacell, 
                               deparse.level = 0)
  obj@experimentData@other$names_metacell <- colnames(metacell)
  
  
  # Imputation
  prot_imp_KNN<-wrapper.impute.slsa(obj) # imputes only partially observed missing values
  prot_imp_2 <- wrapper.impute.detQuant(prot_imp_KNN, na.type = c("Missing MEC", "Missing POV"), qval = 0.025, factor = 1) #missing in entire condition
  
  data_impute <- prot_imp_2@assayData[["exprs"]] %>%
    as.data.frame() 
  print(data_impute)
  
  data_impute_long <- data_impute %>%
    rownames_to_column("ptm_collapse_key") %>%
    pivot_longer(cols = -ptm_collapse_key, names_to = "sample", values_to = "imputed_intensity_log2") %>%
    mutate(imputed_intensity = 2^imputed_intensity_log2) %>%
    left_join(meta_input, by = "sample")
  print(data_impute_long)
  
    # Store the result in the result list
  FR_list_impute[[length(FR_list_impute) + 1]] <- data_impute_long
}
View(FR_list_impute[[1]])

# display intensity distribution after imputation
# Create an empty list to store plots
FR_list_impute_plots <- list()

for (df in FR_list_impute) {
  p1_imputed_distribution <- qc_intensity_distribution(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity_log2 = imputed_intensity_log2,
    plot_style = "histogram"
  )
  # Append the plots to the list
  FR_list_impute_plots[[length(FR_list_impute_plots) + 1]] <- p1_imputed_distribution
}
FR_list_impute_plots


# combine the imputed data list
data_impute <- bind_rows(FR_list_impute)
View(data_impute) 
  
# create matrix
data_impute_matrix <- data_impute %>% 
  dplyr::select(sample, ptm_collapse_key, imputed_intensity_log2) %>%
  pivot_wider(names_from = sample, 
              values_from = imputed_intensity_log2) 

View(data_impute_matrix)

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
dev.off()




