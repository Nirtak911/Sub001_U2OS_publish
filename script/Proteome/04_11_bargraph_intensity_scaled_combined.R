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


## Load data ---------------------------
input<- file.path(data_in,"/04_scale_not_abs_loess_data_log2.csv")
data<- read_csv(input)
input2<- file.path(data_in,"01_genes.csv")
genes<- read_csv(input2)

## output data------------------
plot_output <- file.path(data_in,"Plots/04_11_barplot_stacked.pdf")
plot_output_combined <- file.path(data_in,"Plots/04_11_barplot_stacked_combined.pdf")

# Define custom color palette
custom_colors <- c("Frac1" = "#E7872B", "Frac2" = "#F2B341", 
                   "Frac3" = "#F0E442", "Frac4" = "#009E73", 
                   "Frac5" = "#56B4E9", "Frac6" = "#0072B2",
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

#left join genes
data <- data%>%
  left_join(genes, by = "pg_protein_groups")
View(data)

# sort the treatments for the plot in custom order
data <- data %>%
  mutate(treatment = factor(treatment, levels = c("CTRL", "Anis", "p38", "Zaki")))

# generell abundance change
selected_genes_ab <- data %>%
  filter(str_detect(pg_genes, "JUND|MAPK14|LARP1|TRIM28|TRIM28|STAT1|RB1")) %>%
  pull(pg_genes) %>%
  unique()
print(selected_genes_ab)

# distribution change due to abuncdance change
selected_genes_dis <- data %>%
  filter(str_detect(pg_genes, "LARP1|LMNB1|CENPF|NONO|IGHMBP2|RPS9|SUN1|ITPRIPL2")) %>%
  pull(pg_genes) %>%
  unique()
print(selected_genes_dis)

# real localisation change
selected_genes_loc <- data %>%
  filter(str_detect(pg_genes, "NUMP|LMNA|NUMBL|GGT7|EFNB3|NUMP98")) %>%
  pull(pg_genes) %>%
  unique()
print(selected_genes_loc)

# combine all selections
selected_genes <- c(selected_genes_ab, selected_genes_dis, selected_genes_loc)
selected_genes


# Initialize an empty list to store the generated plots
stacked_intensity_plots <- list()

  # Loop through the list of selected genes and create plots for each
  for (gene in selected_genes ) {
    # Filter the data for the selected gene
    df_selected <- data %>% filter(pg_genes == gene)
    
    # Check if df_selected is empty
    if (nrow(df_selected) == 0) {
      print(paste("No data found for gene:", gene))
      next
    }
    
    # Use the scaled intensity to make a bargraph where all Fractions are stacked in one bar
    # Calculate the mean and standard deviation for each fraction
    df_scaled <- df_selected %>%
      group_by(r_condition, treatment, fraction) %>%
      summarise(
        mean_intensity = mean(scaled_intensity, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(treatment) %>%
      mutate(
        # Rescale the sum of means for each treatment to equal 100%
        pct_protein = (mean_intensity / sum(mean_intensity, na.rm = TRUE)) * 100
      ) %>%
      ungroup()
    
    # Calculate the maximum scaled_intensity
    max_scaled_intensity <- max(df_selected$scaled_intensity, na.rm = TRUE)
    
    # Create the plot
    p2_scaled <- df_scaled %>%
      ggplot(aes(x = treatment, y = pct_protein, fill = fraction)) +
      geom_col(colour = "black", linewidth = 0.5, alpha = 0.8, width = 1, position = "stack") +
      facet_wrap(~treatment, scales = "free_x", nrow = 1, switch = "x") +
      theme_bw() +
      custom_theme +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(axis.title.x = element_blank()) +
      scale_fill_manual(values = custom_colors) +
      scale_colour_manual(values = custom_colors) +
      scale_y_continuous(limits = c(0, 100.001), expand = c(0, 0)) +  # Set y-axis limit to 100%
      ylab("% Protein") +
      ggtitle(paste(gene)) +
      theme(legend.position = "none") # Omit the legend
    p2_scaled
    
    # Separate plot lists
    stacked_intensity_plots[[gene]] <- p2_scaled
  }
  stacked_intensity_plots
  stacked_intensity_plots[1]

# --- Display Plots ---
# You can now display the plots by printing the list.
print(stacked_intensity_plots)


# Combine specific plots into one grid with 2 columns
combined_intensity_plot1 <- grid.arrange(grobs = stacked_intensity_plots[c(1,2,3,4,5,6,7,8,9,10,11,12)], ncol = 4)
combined_intensity_plot2 <- grid.arrange(grobs = stacked_intensity_plots[c(13,14,15,16,17,18,19,20,21,22,23,24)], ncol = 4)
combined_intensity_plot3 <- grid.arrange(grobs = stacked_intensity_plots[c(25,26,27,28,29,1,2,3,4,5,6,7)], ncol = 4)


# Save the combined plots as a PDF file
## Save plots --------------------------
pdf(file = plot_output, width = 20, height = 10)
for (plot in stacked_intensity_plots) {
  print(plot)
}
dev.off()

pdf(file = plot_output_combined, width = 11, height = 7)
grid.draw(combined_intensity_plot1)
grid.newpage()
grid.draw(combined_intensity_plot2)
grid.newpage()
grid.draw(combined_intensity_plot3)
dev.off()

