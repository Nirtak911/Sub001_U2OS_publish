#load libraries
library(tidyverse)
library(protti)
library(data.table)
library(readr)
library(EnhancedVolcano)
library(gridExtra)
library(grid)


## Load data ---------------------------
input <- file.path(data_in,"05_0_diff_data_limma_dapar_impute_complete_group_omnibus.csv")
data <- read_csv(input)
input2 <- file.path(data_in,"01_genes.csv")
genes <- read_csv(input2)
input3 <- "raw/Subcell_markers_selected.csv"
selected_markers <- read_protti(input3)
View(data)

## output data  ---------------------------
plot_output <- file.path(data_in,"/Plots/05_1_diff_expression_limma_final.pdf")
diff_matrix <- file.path(data_in,"05_1_diff_expression_data_wide_limma_final.csv")
sig_matrix <- file.path(data_in,"/05_1_diff_expression_data_wide_sig_limma_final.csv")
sig_count_output <- file.path(data_in,"/05_1_diff_expression_sig_count_limma_final.csv")


## Split data frame by Fractions and cell line------------------------------
split_list <- data %>% group_split(comparison)
print(split_list)


# Create an empty list to store plots
split_list_plots <- list()

for (df in split_list) {
  
  # Function to calculate the p-value cutoff line
  calculate_cutoff <- function(df, adj_pval_cutoff = 0.05, min_for_mean = 5) {
    
    # --- Step 1: Count significant proteins ---
    significant_proteins <- df %>%
      dplyr::filter(adj_pval < adj_pval_cutoff)
    
    num_significant <- nrow(significant_proteins)
    
    
    # --- Step 2: Conditional Cutoff Calculation ---
    if (num_significant >= min_for_mean) {
      
      # TECHNIQUE A: Interpolation (Mean) - Used when there are enough significant proteins
      
      # Filter for the points used for interpolation
      below_cutoff <- df %>% filter(adj_pval <= adj_pval_cutoff)
      above_cutoff <- df %>% filter(adj_pval >= adj_pval_cutoff)
      
      # Find the closest adjusted p-values above and below the cutoff
      closest_below <- below_cutoff %>% filter(adj_pval == max(adj_pval))
      closest_above <- above_cutoff %>% filter(adj_pval == min(adj_pval))
      
      # Calculate the mean only if we have a point closest below AND a point closest above
      if (nrow(closest_below) > 0 & nrow(closest_above) > 0) {
        cutoff_pval <- mean(c(closest_below$pval, closest_above$pval))
      } else {
        # Fallback to max pval if interpolation points aren't available (e.g., all points are significant)
        cutoff_pval <- max(significant_proteins$pval, na.rm = TRUE)
      }
      
    } else if (num_significant > 0) {
      
      # TECHNIQUE B: Max P-value - Used when there is 1 significant protein
      
      cutoff_pval <- max(significant_proteins$pval, na.rm = TRUE)
      
    } else {
      
      # No significant proteins (0 significant proteins)
      cutoff_pval <- NA
    }
    
    return(cutoff_pval)
  }
  
  
  # Calculate the cutoff p-value
  cutoff_pval <- calculate_cutoff(df)
  
  cutoff_log2FC <- 0.585
  
  #volcano plot
  keyvals.colour <- ifelse(
    df$diff < -cutoff_log2FC & df$pval <= cutoff_pval, 'royalblue',
    ifelse(df$diff > cutoff_log2FC & df$pval <= cutoff_pval, 'red',
           'darkgrey'))
  keyvals.colour[is.na(keyvals.colour)] <- 'grey'
  names(keyvals.colour)[keyvals.colour == "red"] <- 'high'
  names(keyvals.colour)[keyvals.colour == 'darkgrey'] <- 'mid'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'low'
  
  volcano_plots <- map(unique(df$comparison), function(cont_idx){
    
    res <- df %>% 
      filter(comparison == cont_idx)
    EnhancedVolcano(res,
                    lab = res$pg_genes,
                    x = 'diff',
                    y = 'pval',
                    title = "",
                    subtitle = paste0(cont_idx),
                    # caption = bquote(~Log[2]~ "fold change = 1; adj_pvalue = 0.05"),
                    legendPosition = "none",
                    labSize = 4,
                    pCutoff =  cutoff_pval,
                    FCcutoff = cutoff_log2FC,
                    colCustom = keyvals.colour,
                    colAlpha = 3/5,
                    # boxedLabels = TRUE,
                    # drawConnectors = TRUE,
                    # max.overlaps = 50
    )
    
  })
  
  # volcano_plots_select <- map(unique(df$comparison), function(cont_idx){
  # 
  #   protein_of_interest <- c("OAT",
  #                            "FYCO1",
  #                            "NEMF",
  #                            "CDK16")
  # 
  #   res <- df %>%
  #     filter(comparison == cont_idx)
  #   EnhancedVolcano(res,
  #                   lab = res$pg_genes,
  #                   x = 'diff',
  #                   y = 'pval',
  #                   title = "",
  #                   subtitle = paste0(cont_idx),
  #                   # caption = bquote(~Log[2]~ "fold change = 1; adj_pvalue = 0.05"),
  #                   legendPosition = "none",
  #                   labSize = 4,
  #                   pCutoff =  cutoff_pval,
  #                   FCcutoff = cutoff_log2FC,
  #                   colCustom = keyvals.colour,
  #                   colAlpha = 3/5,
  #                   selectLab = protein_of_interest,
  #                   boxedLabels = TRUE,
  #                   drawConnectors = TRUE,
  #                   max.overlaps = 50
  #   )
  # 
  # })

  # Append the plots to the list
  split_list_plots[[length(split_list_plots) + 1]] <- volcano_plots
  #split_list_plots[[length(split_list_plots) + 1]] <- volcano_plots_select
}
split_list_plots

split_list_sig <- list()

for (df in split_list) {
  # Filter for significant proteins
  df_sig <- df %>% filter(pval <= cutoff_pval & abs(diff) > cutoff_log2FC)
  
  df_wide_sig <- df_sig %>% 
    pivot_wider(names_from = comparison, 
                values_from = c(diff, pval, adj_pval))
  
  # Append the significant proteins to the list
  split_list_sig[[length(split_list_sig) + 1]] <- df_wide_sig
}
View(split_list_sig[[1]])

# filter for significant proteins
data_sig <- data %>% filter(adj_pval <= 0.05 & abs(diff) > cutoff_log2FC)
View(data_sig)

# wide format of the diff data
df_wide <- data %>% 
  pivot_wider(names_from = comparison, 
              values_from = c(diff, pval, adj_pval))
View(df_wide)

# count significant up and down reguluated proteins per comparison
sig_count <- data_sig %>% 
  group_by(comparison) %>% 
  summarise(sig_up = sum(diff > cutoff_log2FC),
            sig_down = sum(diff < -cutoff_log2FC))
View(sig_count)

# Save csv files --------------------------
write_csv(df_wide, diff_matrix)
write_csv(sig_count, sig_count_output)

# Loop through the list and save each data frame with the 4th column name
for (i in seq_along(split_list_sig)) {
  df_wide_sig <- split_list_sig[[i]]
  col_name_sig <- colnames(df_wide_sig)[11]
  write_csv(df_wide_sig, file.path(data_in,"/diff_expression/", paste("sig_matrix_limma_", col_name_sig, ".csv")))
}

## Save plots --------------------------
pdf(file = plot_output, width = 10, height = 10)
# Loop through each plot in the list
for(i in seq_along(split_list_plots)) {
  # Print the plot to the PDF
  print(split_list_plots[[i]])
}
dev.off()
