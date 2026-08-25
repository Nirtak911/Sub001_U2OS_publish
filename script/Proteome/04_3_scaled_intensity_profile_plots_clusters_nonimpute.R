#load libraries
library(tidyverse)
library(protti)
library(plotly)
library(viridis)
library(ggprism)
library(gridExtra)
library(patchwork)

## Load data---------------------
input <- file.path(data_in,"/4_1_Cluster_pheatmap_scale_data_log2.csv.csv")
data <- read_csv(input)
View(data)
## output data------------------
plot_output <- file.path(data_in,"/Plots/04_3_Profile_plots_cluster_nonimpute.pdf")

View(selected_markers)
## Definitions------------------
# Define custom color palette
custom_colors <- c("Frac1" = "#E7872B", "Frac2" = "#F2B341", 
                   "Frac3" = "#F0E442", "Frac4" = "#009E73", 
                   "Frac5" = "#56B4E9", "Frac6" = "#0072B2",
                   "NA"  = "lightgrey",                    
                   "Cytoplasm" = "#33A02C",            
                   "Cell membrane"  = "#FB9A99",       
                   "Mitochondrion" = "#E31A1C",         
                   "Endoplasmic reticulum" = "#FF7F00",
                   "Golgi apparatus" = "#FDBF6F",
                   "Peroxisome"   = "#CAB2D6",
                   "Lysosome/Vacuole"  = "#6A3D9A", 
                   "Nucleus" = "#1F78B4",   
                   "Extracellular" = "#A6CEE3")

## Add localisation data---------------------
# only keep the first pg_protein_groups
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, "^;", "")
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, ";.*", "")

# remove "Frac" from fraction
data$fraction <- str_replace_all(data$fraction, "Frac", "")

## Profile plots of scaled intensity data---------------------
# split dataset by unique entries in deeploc_single
cluster_list <- data %>% group_split(row_cluster)

# Define the function to process each dataframe
profile_function <- function(df) {
  # Calculate centroid intensity
  centroid_intensity <- df %>%
    group_by(fraction) %>%
    summarize(centroid_intensity = mean(scaled_intensity, na.rm = TRUE))
  
  # Merge centroid intensity with the original dataframe
  df <- df %>%
    left_join(centroid_intensity, by = "fraction")

  # Filter out groups with zero standard deviation
  df_filtered <- df %>%
    group_by(pg_protein_groups, treatment, r_replicate) %>%
    filter(sd(scaled_intensity, na.rm = TRUE) != 0)

  # Calculate Pearson correlation for each replicate
  pearson_correlation <- df_filtered %>%
    group_by(pg_protein_groups, treatment, r_replicate) %>%
    summarize(pearson_correlation = cor(scaled_intensity, centroid_intensity, use = "complete.obs"))

  # Merge Pearson correlation back to the original dataframe
  df_filtered <- df_filtered %>%
    left_join(pearson_correlation, by = c("pg_protein_groups", "treatment", "r_replicate"))
  
  return(df_filtered)
}

#Apply the function to each dataframe in the list deeploc
cluster_list <- map(cluster_list, profile_function)
View(cluster_list[[1]])

## Create plots for the Deeploc data-----------------------------------
# Create an empty list to store the results
cluster_list_plots <- list()

# Create profile plots for each entry in loc_list
for (df in cluster_list) {
  # Get the main location for the title
  main_location_title <- unique(df$row_cluster)
  
  # Calculate the number of unique pg_protein_groups
  num_pg_protein_groups <- df %>% 
    pull(pg_protein_groups) %>% 
    unique() %>% 
    length()
  
  # Create plot
  p1 <- df %>%
    mutate(fraction = as.numeric(fraction)) %>%
    ggplot(aes(x = fraction, y = scaled_intensity, color = pearson_correlation)) +
    geom_line(aes(group = interaction(pg_protein_groups, r_replicate, treatment)), alpha = 0.01) +
    scale_color_gradient2(low = "#00CD00", mid = "yellow", high = "violetred3", limits = c(-1, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = main_location_title,
         x = "Fraction",
         y = "Scaled intensity",
         color = "Pearson Corr.") +
    scale_x_continuous(breaks = seq(1, 6, 1), limits = c(1, 6), expand = c(0, 0))
    # geom_line(aes(x = fraction, y = centroid_intensity, group = 1), color = "#FFBF00", linewidth = 2, linetype = "solid", inherit.aes = FALSE) +
    # annotate("text", x = 5.7, y = 0.97, label = paste("n =", num_pg_protein_groups), size = 2)
  p1
  
  # Append the plots to the list
  cluster_list_plots[[length(cluster_list_plots) + 1]] <- p1
}
cluster_list_plots

# Combine all plots into one grid with a shared legend
combined_cluster_plot <- wrap_plots(cluster_list_plots) + plot_layout(guides = "collect")
combined_cluster_plot


# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 8, height = 3)
combined_cluster_plot
dev.off()