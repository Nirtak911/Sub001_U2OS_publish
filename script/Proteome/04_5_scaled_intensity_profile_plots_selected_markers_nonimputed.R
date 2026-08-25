#load libraries
library(tidyverse)
library(protti)
library(plotly)
library(viridis)
library(ggprism)
library(gridExtra)
library(patchwork)

## Load data---------------------
input <- file.path(data_in,"/04_scale_not_abs_loess_data_log2.csv")
data <- read_csv(input)
input2 <- "raw/Subcell_markers_selected.csv"
selected_markers <- read_protti(input2)

## output data------------------
plot_output <- file.path(data_in,"/Plots/04_5_Profile_plots_selected_markers_CTRLalone_nonimpute.pdf")
plot_output_separate <- file.path(data_in,"/Plots/04_5_Profile_plots_selected_markers_CTRLalone_nonimpute_separate.pdf")

View(data)
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

# filter for CTRL only
data <- data %>% filter(treatment == "CTRL")
# Rename columns if necessary
colnames(selected_markers)[colnames(selected_markers) == "uniprot_id"] <- "pg_protein_groups"
# Add deeploc data
data_select <- data %>% left_join(selected_markers, by = "pg_protein_groups")
View(data_select)

# Filter for each deeploc localisation the top100
data_select <- data_select %>% filter(!is.na(hsap_location_marker))
View(data_select)

## Profile plots of scaled intensity data---------------------
# split dataset by unique entries in deeploc_single
deeploc_list <- data_select %>% group_split(hsap_location_marker)

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
deeploc_list <- map(deeploc_list, profile_function)
View(deeploc_list[[1]])

## Create plots for the Deeploc data-----------------------------------
# Create an empty list to store the results
deeploc_list_plots <- list()

# Create profile plots for each entry in loc_list
for (df in deeploc_list) {
  # Get the main location for the title
  main_location_title <- unique(df$hsap_location_marker)
  
  # Calculate the number of unique pg_protein_groups
  num_pg_protein_groups <- df %>% 
    pull(pg_protein_groups) %>% 
    unique() %>% 
    length()
  
  # Create plot
  p1 <- df %>%
    mutate(fraction = as.numeric(fraction)) %>%
    ggplot(aes(x = fraction, y = scaled_intensity, color = pearson_correlation)) +
    geom_line(aes(group = interaction(pg_protein_groups, r_replicate, treatment)), alpha = 0.5) +
    scale_color_gradient(low = "white", high = "darkblue", limits = c(-1, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(title = main_location_title,
         x = "Fraction",
         y = "Scaled intensity",
         color = "Pearson Corr.") +
    geom_line(aes(x = fraction, y = centroid_intensity, group = 1), color = "#FFBF00", linewidth = 1, linetype = "solid", inherit.aes = FALSE) +
    scale_x_continuous(breaks = seq(1, 6, 1), limits = c(1, 6), expand = c(0, 0)) +
    annotate("text", x = 5.4, y = 0.97, label = paste("n =", num_pg_protein_groups), size = 2)
  p1
  
  # Append the plots to the list
  deeploc_list_plots[[length(deeploc_list_plots) + 1]] <- p1
}
deeploc_list_plots

# Combine all plots into one grid with a shared legend
combined_deeploc_plot <- wrap_plots(deeploc_list_plots) + plot_layout(guides = "collect")
combined_deeploc_plot


# Create an empty list to store the results
separate_list_plots <- list()

for (df in deeploc_list) {
  
  main_location_title <- unique(df$hsap_location_marker)
  
  # Create plot
  p_individual <- df %>%
    mutate(fraction = as.numeric(fraction)) %>%
    ggplot(aes(x = fraction, y = scaled_intensity)) +
    # 1. The Centroid Line (Background reference for the class)
    geom_line(aes(y = centroid_intensity, group = 1), 
              color = "orange", linewidth = 1, linetype = "dashed", alpha = 0.7) +
    # 2. The Individual Protein Lines (Colored by replicate)
    geom_line(aes(group = interaction(r_replicate, treatment), color = as.factor(r_replicate)), 
              linewidth = 0.8) +
    # 3. Facet by protein ID
    facet_wrap(~gene_name, scales = "fixed") + 
    theme_bw() +
    theme(strip.text = element_text(size = 6), # Font size for protein IDs
          axis.text = element_text(size = 6),
          panel.background = element_rect(fill = "white")) +
    scale_color_brewer(palette = "Set1") +
    scale_x_continuous(breaks = seq(1, 6, 1)) +
    labs(title = paste("Location:", main_location_title),
         subtitle = "Dashed line = Class Centroid",
         x = "Fraction",
         y = "Scaled Intensity",
         color = "Replicate")
  
  # Append to list
  separate_list_plots[[length(separate_list_plots) + 1]] <- p_individual
}
separate_list_plots 

# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 10, height = 8)
combined_deeploc_plot
dev.off()

pdf(file = plot_output_separate, width = 20, height = 16)
for(i in seq_along(separate_list_plots )) {
  # Print the plot to the PDF
  print(separate_list_plots [[i]])
}
dev.off()
