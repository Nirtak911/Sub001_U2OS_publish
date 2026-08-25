#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(limma)
library(htmlwidgets)
library(stringi)
library(ggprism)

## Load data---------------------
input <- file.path(data_in, "/02_loess_data_log2.csv")

data <- read_csv(input)
View(data)
## output data------------------
count_output <-  file.path(data_in, "/02_counts.tsv")
plot_output<- file.path(data_in, "/Plots/02_1_QC_loess.pdf")
interactive_plot_output<- file.path(data_in, "/Plots/Interactive/02_1_QC_interactive_loess")

# add normalised intensity not log2 transformed
data$normalised_intensity <- 2^(data$normalised_intensity_log2)

##plot identifications per sample
#IDs per sample table
t1_protein_IDs <- qc_ids(
  data = data,
  sample = sample,
  grouping = ptm_collapse_key,
  intensity = normalised_intensity,
  condition = r_condition,
  title = "Protein identifications per sample",
  plot = FALSE
)
View(t1_protein_IDs)

# add cell line column
t1_protein_IDs$cell_line <- map_chr(str_split(t1_protein_IDs$r_condition, "_"), 1)

# add treatment column
t1_protein_IDs$treatment <- map_chr(str_split(t1_protein_IDs$r_condition, "_"), 2)

# add fraction column
t1_protein_IDs$fraction <- map_chr(str_split(t1_protein_IDs$r_condition, "_"), 3)

# Create the plot
# Ensure r_condition is ordered by fraction
t1_protein_IDs <- t1_protein_IDs %>%
  arrange(fraction, r_condition) %>%
  mutate(r_condition = factor(r_condition, levels = unique(r_condition)))

# split t1_protein_IDs by time
t1_protein_IDs_list <- t1_protein_IDs %>% group_split(cell_line)

ID_plot_list <- list()

for (df in t1_protein_IDs_list) {
  ID_plot <- 
    ggplot(df, aes(x = r_condition, y = count, fill = treatment)) + 
    stat_summary(geom = "col", fun = mean, colour = "black", linewidth = 0.9, alpha = 0.6, width = 0.75) + 
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1), 
                 width = 0.2, colour = "black", linewidth = 0.9) + 
    geom_jitter(aes(colour = treatment), width = 0.2, size = 4, alpha = 0.6) + 
    facet_wrap(~ fraction, scales = "free_x", nrow = 1, switch = "x") +
    scale_x_discrete(labels = function(x) df$fraction[match(x, df$r_condition)]) +
    scale_y_continuous(limits = c(0, 20000), expand = c(0, 0)) +
    theme_prism(base_size = 14) +
    scale_colour_prism(palette = "viridis") + 
    scale_fill_prism(palette = "viridis") +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x=element_blank()) +
    ylab("Localised Phosphosites")+
    theme(strip.text = element_text(size = 14), strip.placement = "outside", strip.background = element_blank())+
    ggtitle(paste(df$time[1]))
  
  ID_plot_list[[length(ID_plot_list)+1]] <- ID_plot
}
ID_plot_list

## Split data frame by lysis-----------------------------
Cell_list <- data %>% group_split(cell_line)
print(Cell_list)

# Create an empty list to store plots
cell_list_plots <- list()

for (df in Cell_list) {
  
  ##calculate and plot CVs-----------------------
  #CV plot
  p1_CV <- qc_cvs(
    data = df,
    grouping = ptm_collapse_key,
    condition = r_condition,
    intensity = normalised_intensity,
    plot = TRUE
  )
  p1_CV
  
  #CV violin plot
  p2_CV <- qc_cvs(
    data = df,
    grouping = ptm_collapse_key,
    condition = r_condition,
    intensity = normalised_intensity,
    plot = TRUE,
    plot_style = "violin"
  )
  p2_CV
  
  
  ##plot identifications per sample
  #IDs per sample plot
  p1_protein_IDs <- qc_ids(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    condition = r_condition,
    title = "Protein identifications per sample",
    plot = TRUE
  )
  p1_protein_IDs
  
  
  # run intensities distribution
  p1_Intensity <- qc_intensity_distribution(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "boxplot"
  )
  p1_Intensity
  
  
  P2_Intensity <- qc_median_intensities(
    data =  df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2
  )
  P2_Intensity
  
  
  # data completeness
  p1_data_completeness <- qc_data_completeness(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    plot = TRUE
  )
  p1_data_completeness
  
  
  # display intensity distribution
  p1_inensity_distribution <- qc_intensity_distribution(
    data = df,
    grouping = ptm_collapse_key,
    intensity_log2 = normalised_intensity_log2,
    plot_style = "histogram"
  )
  p1_inensity_distribution
  
  # display PCA
  p1_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    condition = r_condition,
    digestion = NULL,
    plot_style = "scree"
  )
  p1_pca
  
  p2_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    condition = r_condition,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p2_pca
  
  p3_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    condition = treatment,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p3_pca
  
  p4_pca <- qc_pca(
    data = df,
    sample = sample,
    grouping = ptm_collapse_key,
    intensity = normalised_intensity_log2,
    condition = fraction,
    components = c("PC1", "PC2"),
    plot_style = "pca"
  )
  p4_pca
  
  # Plot ranked peptide intensities
  p1_ranked <- qc_ranked_intensities(
    data = df,
    sample = sample,
    grouping = gene_PTM_key,
    intensity_log2 = normalised_intensity_log2,
    plot = TRUE,
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
  grouping = ptm_collapse_key,
  intensity_log2 = normalised_intensity_log2,
  condition = fraction,
  interactive = TRUE
)
p1_sample_correlation

# Append the plots to the list
cell_list_interactive[[length(cell_list_interactive) + 1]] <- p1_sample_correlation
}
cell_list_interactive

# Save data
write_tsv(t1_protein_IDs, count_output)

## Save plots --------------------------
pdf(file = plot_output, width = 18, height = 10)
# Loop through each plot in the list
for(i in seq_along(ID_plot_list)) {
  # Print the plot to the PDF
  print(ID_plot_list[[i]])
}
# Loop through each plot in the list
for(i in seq_along(cell_list_plots)) {
  # Print the plot to the PDF
  print(cell_list_plots[[i]])
}
dev.off()
# 
# # Save the correlation to an HTML file
# for(i in seq_along(cell_list_interactive)) {
#   # Save the plot to an HTML file
#   htmlwidgets::saveWidget(cell_list_interactive[[i]], file = paste0(interactive_plot_output, i, ".html"))
# }
# 

