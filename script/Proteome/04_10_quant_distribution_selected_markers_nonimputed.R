#load libraries
library(tidyverse)
library(protti)
library(plotly)
library(viridis)
library(ggprism)
library(gridExtra)
library(patchwork)

## Load data---------------------
input <- file.path(data_in,"/02_loess_data_log2.csv")
data <- read_csv(input)
input2 <- "raw/Subcell_markers_selected.csv"
selected_markers <- read_protti(input2)
input3 <- file.path(data_in,"/01_genes.csv")
genes <- read_csv(input3)

# ## output data------------------
plot_output <- file.path(data_in,"/Plots/04_10_quant_distribution_selected_markers_manual_CTRLalone_nonimpute.pdf")

View(data)
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
data$normalised_intensity <- 2^(data$normalised_intensity_log2)

# only keep the first pg_protein_groups
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, "^;", "")
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, ";.*", "")

# remove "Frac" from fraction
data$fraction <- str_replace_all(data$fraction, "Frac", "")

# filter for CTRL only
data <- data %>% filter(treatment == "CTRL")
# Rename columns if necessary
#colnames(selected_markers)[colnames(selected_markers) == "uniprot_id"] <- "pg_protein_groups"
colnames(selected_markers)[colnames(selected_markers) == "gene_name"] <- "pg_genes"
# Add deeploc data
data_select <- data %>% left_join(selected_markers, by = "pg_genes")
View(data_select)

# Filter for each deeploc localisation the top100
data_select <- data_select %>% filter(!is.na(hsap_location_marker))
View(data_select)


# 1. Calculate the summed compartment intensities and proportions per fraction
compartment_proportions <- data_select %>%
  group_by(fraction, hsap_location_marker) %>%
  # Sum the linear-scale intensities for each compartment in each fraction
  summarize(compartment_sum = sum(normalised_intensity, na.rm = TRUE), .groups = "drop_last") %>%
  # Divide by the total marker intensity inside that specific fraction
  mutate(proportion = compartment_sum / sum(compartment_sum)) %>%
  ungroup() %>%
  # Ensure fraction is treated as a clean factor for ordering on the X-axis
  mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4", "5", "6")))

# 2. Plot as a Stacked Bar Plot
stacked_proportion_plot <- ggplot(compartment_proportions, 
                                  aes(x = fraction, y = proportion, fill = hsap_location_marker)) +
  # geom_col automatically stacks when a fill aesthetic is provided
  geom_col(width = 0.7, color = "black", linewidth = 0.3) + 
  
  # Styling and axis limits
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  #cale_fill_manual(values = custom_colors) +
  theme_prism(base_size = 11) + # Provides clean, publication-ready borders and text
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Subcellular Compartment Distribution per Fraction",
    x = "Fraction",
    y = "Relative Intensity Proportion (%)",
    fill = "Compartment"
  )

# Display the plot
print(stacked_proportion_plot)


# 1. Calculate the distribution of each compartment across the 6 fractions
compartment_across_fractions <- data_select %>%
  group_by(fraction, hsap_location_marker) %>%
  # Sum raw intensities per compartment per fraction
  summarize(compartment_sum = sum(normalised_intensity, na.rm = TRUE), .groups = "drop") %>%
  
  # Group by the compartment marker to normalize across the fractions
  group_by(hsap_location_marker) %>%
  # Divide the fraction's intensity by the sum of all fractions for that compartment
  mutate(fraction_distribution = compartment_sum / sum(compartment_sum)) %>%
  ungroup() %>%
  
  # Ensure fraction is ordered correctly
  mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4", "5", "6")))

# 2. Visualize as a Faceted Bar Plot to see the profile of each compartment
p2 <- ggplot(compartment_across_fractions, aes(x = fraction, y = fraction_distribution, fill = hsap_location_marker)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3, show.legend = FALSE) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  #scale_fill_manual(values = custom_colors) +
  facet_wrap(~hsap_location_marker, scales = "free_y") +  # Creates a separate profile plot for each organelle
  theme_prism(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Relative Elution Profile of Subcellular Compartments Across Fractions",
    x = "Fraction",
    y = "Proportion of Total Compartment Intensity (%)"
  )

p2

# 1. Calculate the distribution of each compartment across the 6 fractions
compartment_across_fractions <- data_select %>%
  group_by(fraction, hsap_location_marker) %>%
  # Sum raw intensities per compartment per fraction
  summarize(compartment_sum = sum(normalised_intensity, na.rm = TRUE), .groups = "drop") %>%
  
  # Group by the compartment marker to normalize across the fractions
  group_by(hsap_location_marker) %>%
  # Divide the fraction's intensity by the sum of all fractions for that compartment
  mutate(fraction_distribution = compartment_sum / mean(compartment_sum)) %>%
  ungroup() %>%
  
  # Ensure fraction is ordered correctly
  mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4", "5", "6")))

# 2. Visualize as a Faceted Bar Plot to see the profile of each compartment
p3 <- ggplot(compartment_across_fractions, aes(x = fraction, y = fraction_distribution, fill = hsap_location_marker)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3, show.legend = FALSE) +
  scale_y_continuous( expand = c(0, 0)) +
  #scale_fill_manual(values = custom_colors) +
  facet_wrap(~hsap_location_marker, scales = "free_y") +  # Creates a separate profile plot for each organelle
  theme_prism(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Enrichment of Subcellular Compartments Across Fractions",
    x = "Fraction",
    y = "Enrichent factor"
  )
p3

# 1. Calculate both metrics in one data pipeline
compartment_profiles <- data_select %>%
  group_by(fraction, hsap_location_marker) %>%
  summarize(compartment_sum = sum(normalised_intensity, na.rm = TRUE), .groups = "drop") %>%
  
  group_by(hsap_location_marker) %>%
  mutate(
    # Left Axis Metric: Proportion (0 to 1)
    fraction_distribution = compartment_sum / sum(compartment_sum),
    # Right Axis Metric: Enrichment factor
    enrichment_factor = compartment_sum / mean(compartment_sum)
  ) %>%
  ungroup() %>%
  mutate(fraction = factor(fraction, levels = c("1", "2", "3", "4", "5", "6")))

# 2. Define the transformation relationship
# Since we have 6 fractions: enrichment_factor = proportion * 6
scaling_factor <- 6

# 3. Create the Dual Y-Axis Plot with Fixed Grid Dimensions
p_combined <- ggplot(compartment_profiles, aes(x = fraction, y = fraction_distribution, fill = hsap_location_marker)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3, show.legend = FALSE) +
  
  # Dual Y-Axis
  scale_y_continuous(
    name = "Proportion of Total Compartment Intensity (%)",
    labels = scales::percent_format(),
    expand = c(0, 0),
    sec.axis = sec_axis(~ . * scaling_factor, name = "Enrichment Factor")
  ) +
  
  # FIXED DIMENSIONS: Added nrow = 3 and ncol = 5
  facet_wrap(~hsap_location_marker, scales = "free_y", nrow = 3, ncol = 5) +  
  
  theme_prism(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y.right = element_text(vjust = 1.5, size = 9)
  ) +
  labs(
    title = "Subcellular Compartment Profiles & Enrichment Across Fractions",
    x = "Fraction"
  )

p_combined

# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 12, height = 6)
stacked_proportion_plot 
p2
p3
p_combined
dev.off()

