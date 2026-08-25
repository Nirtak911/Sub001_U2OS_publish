#load libraries
library(tidyverse)
library(protti)
library(pheatmap)
library(metap)
library(EnhancedVolcano)
library(patchwork)
library(gridExtra)
library(grid)


## Load data ---------------------------
input<- file.path(data_in,"/04_0_scale_not_abs_data_log2.csv")
data<- read_csv(input)
input2<- file.path(data_in,"01_genes.csv")
genes<- read_csv(input2)

## output data------------------
plot_output <- file.path(data_in,"Plots/04_3_barplot_stacked_scaled_combined.pdf")


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



# Define a custom theme with the desired font size
custom_theme <- theme(
  text = element_text(size = 14),  # Adjust the size as needed
  plot.title = element_text(size = 16, hjust = 0.5),
  axis.text.x = element_text(size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14)
)
View(data)

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

# sort the treatments for the plot in custom order
data <- data %>%
  mutate(treatment = factor(treatment, levels = c("CTRL", "Anis", "p38", "Zaki")))

# Initialize an empty list to store the plots
stacked_intensity_plots <- list()

# Loop through the list of selected genes and create plots for each
for (gene in selected_genes) {
  # Filter the data for the selected gene
  df_selected <- data %>% filter(gene_PTM_key == gene)

  # Check if df_selected is empty
  if (nrow(df_selected) == 0) {
    print(paste("No data found for gene:", gene))
    next
  }
  
  # Calculate mean intensity per fraction and convert to 100% relative proportion
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
    ylab("% Phosphosite") +
    ggtitle(paste(gene)) +
    theme(legend.position = "none") # Omit the legend
  p2_scaled
 
   # Separate plot lists
  stacked_intensity_plots[[length(stacked_intensity_plots) + 1]] <- p2_scaled
}
stacked_intensity_plots

# combined_stacked_plot
combined_stacked_intensity_plot1 <- grid.arrange(grobs = stacked_intensity_plots[c(1,2,3,4,5,6,7,8,9, 10, 11, 12)], ncol = 4)
combined_stacked_intensity_plot2 <- grid.arrange(grobs = stacked_intensity_plots[c(13,14,15,16,17,18,19,20,21,22,23,24)], ncol = 4)
combined_stacked_intensity_plot3 <- grid.arrange(grobs = stacked_intensity_plots[c(25,26,27,1,2,3,4,5,6,7,8,9)], ncol = 4)

# Save plots --------------------------
pdf(file = plot_output, width = 11, height = 7)
grid.draw(combined_stacked_intensity_plot1)
grid.newpage()
grid.draw(combined_stacked_intensity_plot2)
grid.newpage()
grid.draw(combined_stacked_intensity_plot3)
dev.off()

