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
plot_output <- file.path(data_in,"/Plots/04_9_Profile_plots_selected_markers_nonimpute_Treatments_centroids.pdf")


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


profile_function <- function(df) {
  
  # 1. Master Centroid (all replicates combined)
  master_centroid <- df %>%
    group_by(fraction) %>%
    summarize(master_intensity = mean(scaled_intensity, na.rm = TRUE), .groups = "drop")
  
  # 2. Individual Replicate Centroids
  treatment_centroids <- df %>%
    group_by(fraction, treatment) %>%
    summarize(treat_centroid_intensity = mean(scaled_intensity, na.rm = TRUE), .groups = "drop")
  
  # Merge both back to the main dataframe
  df <- df %>%
    left_join(master_centroid, by = "fraction") %>%
    left_join(treatment_centroids, by = c("fraction", "treatment"))
  
  # Filter out groups with zero standard deviation
  df_filtered <- df %>%
    group_by(pg_protein_groups, treatment, r_replicate) %>%
    filter(sd(scaled_intensity, na.rm = TRUE) != 0)
  
  # Calculate Pearson correlation using the replicate-specific centroid
  pearson_correlation <- df_filtered %>%
    group_by(pg_protein_groups, treatment, r_replicate) %>%
    summarize(
      pearson_correlation = cor(scaled_intensity, treat_centroid_intensity, use = "complete.obs"),
      .groups = "drop"
    )
  
  # Merge Pearson correlation back
  df_filtered <- df_filtered %>%
    left_join(pearson_correlation, by = c("pg_protein_groups", "treatment", "r_replicate"))
  
  return(df_filtered)
}

# Apply the updated function
deeploc_list <- map(deeploc_list, profile_function)

# Define a function to generate the profile plot for each localization group
plot_profile_group <- function(df) {
  
  # Extract the localization name for titles and color mapping
  loc_name <- unique(df$hsap_location_marker)[1]
  line_color <- custom_colors[loc_name]
  if (is.na(line_color)) line_color <- "black" # Fallback if not in custom_colors
  
  n_proteins <- n_distinct(df$pg_protein_groups)
  
  # Ensure fraction is a factor or numeric for correct line drawing ordering
  df <- df %>% mutate(fraction = as.numeric(fraction))
  
  ggplot(df, aes(x = fraction)) +
    # 1. Background individual protein tracks
    geom_line(
      aes(y = scaled_intensity, group = interaction(pg_protein_groups, treatment)), 
      color = "grey80", alpha = 0.15, linewidth = 0.4
    ) +
    
    # 2. Individual Replicate Centroids (FIXED: forced to factor)
    geom_line(
      aes(
        y = treat_centroid_intensity, 
        group = treatment, 
        color = factor(treatment)),
      linewidth = 1.0,
      alpha = 0.9
    ) +
    
    # # 3. The Master Centroid of everything combined
    # geom_line(
    #   aes(y = master_intensity, group = 1),
    #   color = "#FFBF00", 
    #   linewidth = 1.5,
    #   linetype = "solid"
    # ) +
    # 
    # Adding aesthetics and formatting match standard publication styles
    theme_prism(base_size = 8) +
    labs(
      title = paste0(loc_name, " (n = ", n_proteins, ")"),
      x = "Fraction",
      y = "Scaled Intensity",
      linetype = "Treatment"  # Cleans up the legend title
    ) +
    scale_x_continuous(breaks = unique(df$fraction))
}

# Generate all plots again
plot_list <- map(deeploc_list, plot_profile_group)
plot_list

# Combine all plots into one grid with a shared legend
combined_plot <- wrap_plots(plot_list) + plot_layout(guides = "collect")
combined_plot


# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 10, height = 8)
combined_plot
dev.off()
