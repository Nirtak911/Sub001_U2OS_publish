# Load libraries
library(tidyverse)
library(protti)
library(pheatmap)
library(metap)
library(EnhancedVolcano)
library(patchwork)
library(gridExtra)
library(grid)
library(purrr)
library(stringr)

# Load data
input <- file.path(data_in, "02_loess_data_log2.csv")
data <- read_csv(input)
input2 <- file.path(data_in, "01_genes.csv")
genes <- read_csv(input2)

# Output data
plot_output <- file.path(data_in, "Plots/02_2_mean_intensity_bargraph.pdf")
plot_output_combined <- file.path(data_in, "Plots/02_2_mean_intensity_bargraph_selected_combined.pdf")

# Define custom color palette
custom_colors <- c("Frac1" = "#E7872B", "Frac2" = "#F2B341", 
                   "Frac3" = "#F0E442", "Frac4" = "#009E73", 
                   "Frac5" = "#56B4E9", "Frac6" = "#0072B2",
                   "CTRL" = "#abdbe3", "Anis" = "#063970", 
                   "p38" = "#f46d43", "Zaki" = "#d73027",
                   "NA" = "lightgrey",
                   "Cytoplasm" = "#33A02C",
                   "Cell membrane" = "#FB9A99",
                   "Mitochondrion" = "#E31A1C",
                   "Endoplasmic reticulum" = "#FF7F00",
                   "Golgi apparatus" = "#FDBF6F",
                   "Peroxisome" = "#CAB2D6",
                   "Lysosome/Vacuole" = "#6A3D9A", 
                   "Nucleus" = "#1F78B4", 
                   "Extracellular" = "#A6CEE3")

# Define a custom theme with the desired font size
custom_theme <- theme(
  text = element_text(size = 14), # Adjust the size as needed
  plot.title = element_text(size = 16, hjust = 0.5),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14)
)

View(data)

data$normalised_intensity <- 2^data$normalised_intensity_log2

#add time column 
data$gene <- map_chr(str_split(data$gene_PTM_key, "_"), 1)

# sort the treatments for the plot in custom order
data <- data %>%
  mutate(treatment = factor(treatment, levels = c("CTRL", "Anis", "p38", "Zaki")))

# generell abundance change
selected_genes_ab <- data %>%
  filter(str_detect(gene_PTM_key, "JUND_S100|MAPK14_T180|LARP1_S766|TRIM28_T541|TRIM28_S473|STAT1_S727|RB1_T356|RB1_S249")) %>%
  pull(gene_PTM_key) %>%
  unique()
print(selected_genes_ab)

# distribution change due to abuncdance change
selected_genes_dis <- data %>%
  filter(str_detect(gene_PTM_key, "LARP1_T526|LMNB1_S23_M1|CENPF_S3023_M1|NONO_T450|IGHMBP2_S716_M1|RPS9_S146_M1|SUN1_S48_M1|ITPRIPL2_S78_M1")) %>%
  pull(gene_PTM_key) %>%
  unique()
print(selected_genes_dis)

# real localisation change
selected_genes_loc <- data %>%
  filter(str_detect(gene_PTM_key, "NUMP_214_S678|LMNA_S390_M2|LMNA_S75_M1|NUMBL_S305|GGT7_S83|EFNB3_S268_M1|NUMP98_S1043")) %>%
  pull(gene_PTM_key) %>%
  unique()
print(selected_genes_loc)

# combine all selections
selected_genes <- c(selected_genes_ab, selected_genes_dis, selected_genes_loc)
selected_genes

#Filter data by treatment
data_filtered <- data
# filter(treatment %in% c("IFNy", "CTRL"))

data_filtered$phosphosite <- map_chr(str_split(data_filtered$gene_PTM_key, "_"), 2)
data_filtered$multi <- map_chr(str_split(data_filtered$gene_PTM_key, "_"), 3)
View(data_filtered)


# --- Pre-filter Data and Identify Keys for Plotting ---
# 1. Filter the main data frame to keep only the rows matching your selected_genes list.
df_pre_filtered <- data_filtered %>%
  filter(gene_PTM_key %in% selected_genes)

# 2. Get the list of unique gene_PTM_keys from this newly filtered data frame.
#    The loop will iterate through these keys.
keys_to_plot <- unique(df_pre_filtered$gene_PTM_key)

# Initialize an empty list to store the generated plots
intensity_plots <- list()

# --- Loop through Genes and Create Plots ---
# This loop now iterates over each unique key found in your pre-filtered data.
for (current_key in keys_to_plot) {
  
  # Filter the pre-filtered data frame to get the data for the current key only.
  # This is a very efficient way to select the data for each plot.
  df_selected <- df_pre_filtered %>%
    filter(gene_PTM_key == current_key)
  
  # A check to ensure data exists (though it should, given the logic).
  if (nrow(df_selected) == 0) {
    print(paste("No data found for key:", current_key)) # This should not happen with the new logic
    next
  }
  
  # Calculate median and standard deviation for the current key
  df_mean_summary <- df_selected %>%
    group_by(r_condition, fraction, treatment) %>%
    summarise(
      mean_intensity = mean(normalised_intensity, na.rm = TRUE),
      sd_intensity = sd(normalised_intensity, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate the maximum intensity for the current key to set the y-axis limit
  max_normalised_intensity <- max(df_selected$normalised_intensity, na.rm = TRUE)
  
  
  # Create the plot for the current key
  p2_mean_intensity <- df_mean_summary %>%
    ggplot(aes(x = r_condition, y = mean_intensity, fill = fraction)) +
    
    # Bar plot representing the median intensity
    geom_col(colour = "black", linewidth = 0.75, alpha = 0.6, width = 0.75) +
    
    # Error bars representing the standard deviation
    geom_errorbar(
      aes(ymin = pmax(0, mean_intensity - sd_intensity), ymax = mean_intensity + sd_intensity),
      width = 0.2, 
      colour = "black", 
      linewidth = 0.75
    ) +
    
    # Jittered points showing all individual data points for the current key
    geom_jitter(
      data = df_selected, 
      aes(x = r_condition, y = normalised_intensity, colour = fraction),
      width = 0.2, 
      size = 2, 
      alpha = 0.6
    ) +
    
    # Facet by treatment to separate the plots
    facet_wrap(~treatment, scales = "free_x", nrow = 1, switch = "x") +
    
    # Apply themes and styling
    theme_bw() +
    custom_theme +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    
    # Apply custom colors
    scale_fill_manual(values = custom_colors) +
    scale_colour_manual(values = custom_colors) +
    
    # # Set y-axis limits and labels
    scale_y_continuous(limits = c(0, max_normalised_intensity * 1.05), expand = c(0, 0)) +
    ylab("Intensity") +
    
    # Add the key name as the title
    ggtitle(paste(current_key))
  
  # plot log2
  # Calculate median and standard deviation for the current key
  df_mean_summary_log2 <- df_selected %>%
    group_by(r_condition, fraction, treatment) %>%
    summarise(
      mean_intensity_log2 = mean(normalised_intensity_log2, na.rm = TRUE),
      sd_intensity_log2 = sd(normalised_intensity_log2, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate the maximum intensity for the current key to set the y-axis limit
  max_normalised_intensity_log2 <- max(df_selected$normalised_intensity_log2, na.rm = TRUE)
  
  # Create the plot for the current key
  p2_mean_intensity_log2 <- df_mean_summary_log2 %>%
    ggplot(aes(x = r_condition, y = mean_intensity_log2, fill = fraction)) +
    
    # Bar plot representing the median intensity
    geom_col(colour = "black", linewidth = 0.75, alpha = 0.6, width = 0.75) +
    
    # Error bars representing the standard deviation
    geom_errorbar(
      aes(ymin = pmax(0, mean_intensity_log2 - sd_intensity_log2), ymax = mean_intensity_log2 + sd_intensity_log2),
      width = 0.2, 
      colour = "black", 
      linewidth = 0.75
    ) +
    
    # Jittered points showing all individual data points for the current key
    geom_jitter(
      data = df_selected, 
      aes(x = r_condition, y = normalised_intensity_log2, colour = fraction),
      width = 0.2, 
      size = 2, 
      alpha = 0.6
    ) +
    
    # Facet by treatment to separate the plots
    facet_wrap(~treatment, scales = "free_x", nrow = 1, switch = "x") +
    
    # Apply themes and styling
    theme_bw() +
    custom_theme +
    theme(
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    
    # Apply custom colors
    scale_fill_manual(values = custom_colors) +
    scale_colour_manual(values = custom_colors) +
    
    # Set y-axis limits and labels
    scale_y_continuous(limits = c(0, max_normalised_intensity_log2 * 1.05), expand = c(0, 0)) +
    ylab("log2(Intensity)") +
    
    # Add the key name as the title
    ggtitle(paste(current_key))
  # Add the newly created plot to our list, using the key as the name
  intensity_plots[[paste0(current_key, "_raw")]] <- p2_mean_intensity
  intensity_plots[[paste0(current_key, "_log2")]] <- p2_mean_intensity_log2
}
print(intensity_plots)

# Combine specific plots into one grid with 2 columns
combined_intensity_plot1 <- grid.arrange(grobs = intensity_plots[c(1,2,3,4)], ncol = 2)
combined_intensity_plot2 <- grid.arrange(grobs = intensity_plots[c(5,6,7,8)], ncol = 2)
combined_intensity_plot3 <- grid.arrange(grobs = intensity_plots[c(9,10,11,12)], ncol = 2)
combined_intensity_plot4 <- grid.arrange(grobs = intensity_plots[c(13,14,15,16)], ncol = 2)
combined_intensity_plot5 <- grid.arrange(grobs = intensity_plots[c(17,18,19,20)], ncol = 2)
combined_intensity_plot6 <- grid.arrange(grobs = intensity_plots[c(21,22,23,24)], ncol = 2)
combined_intensity_plot7 <- grid.arrange(grobs = intensity_plots[c(25,26,27,28)], ncol = 2)
combined_intensity_plot8 <- grid.arrange(grobs = intensity_plots[c(29,30,31,32)], ncol = 2)
combined_intensity_plot9 <- grid.arrange(grobs = intensity_plots[c(33,34,35,36)], ncol = 2)
combined_intensity_plot10 <- grid.arrange(grobs = intensity_plots[c(37,38,39,40)], ncol = 2)
combined_intensity_plot11 <- grid.arrange(grobs = intensity_plots[c(41,42,43,44)], ncol = 2)
combined_intensity_plot12 <- grid.arrange(grobs = intensity_plots[c(45,46,47,48)], ncol = 2)
combined_intensity_plot13 <- grid.arrange(grobs = intensity_plots[c(49,50,51,52)], ncol = 2)
combined_intensity_plot14 <- grid.arrange(grobs = intensity_plots[c(53,1,2,3)], ncol = 2)


# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 20, height = 10)
for (plot in intensity_plots) {
  print(plot)
}
dev.off()

pdf(file = plot_output_combined, width = 10, height = 5)
grid.draw(combined_intensity_plot1)
grid.newpage()
grid.draw(combined_intensity_plot2)
grid.newpage()
grid.draw(combined_intensity_plot3)
grid.newpage()
grid.draw(combined_intensity_plot4)
grid.newpage()
grid.draw(combined_intensity_plot5)
grid.newpage()
grid.draw(combined_intensity_plot6)
grid.newpage()
grid.draw(combined_intensity_plot7)
grid.newpage()
grid.draw(combined_intensity_plot8)
grid.newpage()
grid.draw(combined_intensity_plot9)
grid.newpage()
grid.draw(combined_intensity_plot10)
grid.newpage()
grid.draw(combined_intensity_plot11)
grid.newpage()
grid.draw(combined_intensity_plot12)
grid.newpage()
grid.draw(combined_intensity_plot13)
grid.newpage()
grid.draw(combined_intensity_plot14)
dev.off()
# 
