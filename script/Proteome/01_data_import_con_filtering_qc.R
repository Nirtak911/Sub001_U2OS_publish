#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(limma)
library(htmlwidgets)
library(stringi)

## Load data---------------------
input <- file.path(raw_data, "20250221_110126_ProteinGroups_KS_normal (Normal).tsv")
data <- read_protti(input)
View(data)

## output data------------------
filter_output <- file.path(data_in, "/01_data_log2.csv")
raw_matrix <- file.path(data_in, "/01_data_log2_matrix.csv")
plot_output <-  file.path(data_in, "/Plots/01_QC.pdf")
interactive_plot_output <-  file.path(data_in, "/Plots/Interactive/01_QC_interactive")
count_output <-  file.path(data_in, "/01_counts.tsv")
meta_output <- file.path(data_in, "/01_meta_data.csv")
genes_output <- file.path(data_in, "/01_genes.csv")

## Data filtering and QC---------------------
# add sample column
data$r_file_name <- str_replace_all(data$r_file_name, "UTOS", "U2OS")
data$r_condition <- str_replace_all(data$r_condition, "UTOS", "U2OS") 
data$sample <- paste(data$r_condition, data$r_replicate, sep = "_")

# add cell line column
data$cell_line <- map_chr(str_split(data$sample, "_"), 1)

# add treatment column
data$treatment <- map_chr(str_split(data$sample, "_"), 2)

# add fraction column
data$fraction <- map_chr(str_split(data$sample, "_"), 3)

# filter data
data_filter <- data %>% 
  filter(!str_detect(pg_protein_groups, "CON_")) #filter contaminants

# remove extra ;
data_filter$pg_genes <- str_replace_all(data_filter$pg_genes, "^;", "")
data_filter$pg_protein_groups <- str_replace_all(data_filter$pg_protein_groups, "^;", "")

# only keep the first gene name
data_filter$pg_genes <- str_replace_all(data_filter$pg_genes, ";.*", "")
data_filter$pg_protein_groups <- str_replace_all(data_filter$pg_protein_groups, ";.*", "")

#log2 transform PTM Quantity column
data_filter$log2_pg_quantity <- log2(data_filter$pg_quantity)
View(data_filter)

# data_filter matrix
data_filter_matrix <- data_filter %>% 
  dplyr::select(sample, pg_protein_groups, log2_pg_quantity) %>%
  pivot_wider(names_from = sample, 
              values_from = log2_pg_quantity,
              values_fn = list(log2_pg_quantity = ~ mean(.x, na.rm = TRUE)))
View(data_filter_matrix)

# get metadata
meta_data <- data_filter %>% 
  dplyr::select(sample, cell_line, treatment, fraction, r_replicate, r_condition, r_file_name) %>%
  distinct()
View(meta_data)

# get pg_genes and corresponding pg_protein_groups
genes <- data_filter %>% 
  dplyr::select(pg_protein_groups, pg_genes) %>%
  distinct()

# data_output
data_filter_output <- data_filter_matrix %>%
                      pivot_longer(cols = -pg_protein_groups, 
                                   names_to = "sample", 
                                   values_to = "log2_pg_quantity") %>%
                      left_join(meta_data, by = "sample") %>%
                      mutate(pg_quantity = 2^log2_pg_quantity) %>%
                      left_join(genes, by = "pg_protein_groups")
View(data_filter_output)

##calculate and plot CVs
#CV table
t1_CV <- qc_cvs(
  data = data_filter_output,
  grouping = pg_protein_groups,
  condition = r_condition,
  intensity = pg_quantity,
  plot = FALSE
)
View(t1_CV)

##plot identifications per sample
#IDs per sample table
t1_protein_IDs <- qc_ids(
  data = data_filter_output,
  sample = sample,
  grouping = pg_protein_groups,
  intensity = pg_quantity,
  condition = r_condition,
  title = "Protein identifications per sample",
  plot = FALSE
)
View(t1_protein_IDs)


## Split data frame by cell line------------------------------
Cell_list <- data_filter_output %>% group_split(cell_line)
print(Cell_list)

# Create an empty list to store plots
cell_list_plots <- list()

for (df in Cell_list) {
  
  ##calculate and plot CVs-----------------------
  #CV plot
  p1_CV <- qc_cvs(
    data = df,
    grouping = pg_protein_groups,
    condition = r_condition,
    intensity = pg_quantity,
    plot = TRUE
  )
  p1_CV
  
  #CV violin plot
  p2_CV <- qc_cvs(
    data = df,
    grouping = pg_protein_groups,
    condition = r_condition,
    intensity = pg_quantity,
    plot = TRUE,
    plot_style = "violin"
  )
  p2_CV
  
  
  ##plot identifications per sample
  #IDs per sample plot
  p1_protein_IDs <- qc_ids(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    condition = r_condition,
    title = "Protein identifications per sample",
    plot = TRUE
  )
  p1_protein_IDs
  
  
  # run intensities distribution
  p1_Intensity <- qc_intensity_distribution(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity_log2 = log2_pg_quantity,
    plot_style = "boxplot"
  )
  p1_Intensity
  
  
  P2_Intensity <- qc_median_intensities(
    data =  df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity
  )
  P2_Intensity
  
  
  # data completeness
  p1_data_completeness <- qc_data_completeness(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    plot = TRUE
  )
  p1_data_completeness
  
  
  # display intensity distribution
  p1_inensity_distribution <- qc_intensity_distribution(
    data = df,
    grouping = pg_protein_groups,
    intensity_log2 = log2_pg_quantity,
    plot_style = "histogram"
  )
  p1_inensity_distribution
  
  # display PCA
  p1_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    condition = r_condition,
    digestion = NULL,
    plot_style = "scree"
  )
  p1_pca
  
  p2_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    condition = r_condition,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p2_pca
  
  p3_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    condition = treatment,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p3_pca
  
  p4_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity = log2_pg_quantity,
    condition = fraction,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p4_pca
  
  # Plot ranked peptide intensities
  p1_ranked <- qc_ranked_intensities(
    data = df,
    sample = sample,
    grouping = pg_genes,
    intensity_log2 = log2_pg_quantity,
    plot = TRUE
  )
  p1_ranked
  
  # Append the plots to the list
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_CV
  cell_list_plots[[length(cell_list_plots) + 1]] <- p2_CV
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_protein_IDs
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_Intensity
  cell_list_plots[[length(cell_list_plots) + 1]] <- P2_Intensity
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_data_completeness
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_inensity_distribution
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_pca
  cell_list_plots[[length(cell_list_plots) + 1]] <- p2_pca
  cell_list_plots[[length(cell_list_plots) + 1]] <- p3_pca
  cell_list_plots[[length(cell_list_plots) + 1]] <- p4_pca
  cell_list_plots[[length(cell_list_plots) + 1]] <- p1_ranked
}
cell_list_plots


# Create an empty list to store interactive plots
cell_list_interactive <- list()

for (df in Cell_list) {
  
  # display correlation
  p1_sample_correlation <- qc_sample_correlation(
    data = df,
    sample = sample,
    grouping = pg_protein_groups,
    intensity_log2 = log2_pg_quantity,
    condition = fraction,
    interactive = TRUE
  )
  p1_sample_correlation
  
  # Append the plots to the list
  cell_list_interactive[[length(cell_list_interactive) + 1]] <- p1_sample_correlation
}
cell_list_interactive


## Save data ---------------------------
write_csv(data_filter_output, filter_output)
write_csv(data_filter_matrix, raw_matrix)
write_csv(meta_data, meta_output)
write_csv(genes, genes_output)
write_tsv(t1_protein_IDs, count_output)

## Save plots --------------------------
pdf(file = plot_output, width = 18, height = 10)
# Loop through each plot in the list
for(i in seq_along(cell_list_plots)) {
  # Print the plot to the PDF
  print(cell_list_plots[[i]])
}
dev.off()

# Save the correlation to an HTML file
for(i in seq_along(cell_list_interactive)) {
  # Save the plot to an HTML file
  htmlwidgets::saveWidget(cell_list_interactive[[i]], file = paste0(interactive_plot_output, i, ".html"))
}
