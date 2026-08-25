#load libraries
library(tidyverse)
library(protti)
library(ggplot2)
library(clusterProfiler)
library(forcats)
library(DOSE)
library(gridExtra)
library(patchwork)

## Load data ---------------------------
input<- file.path(data_in,"/4_1_gsea_cluster_pheatmap_scale_data_log2.csv")
data<- read_protti(input)
View(data)

# output data  ---------------------------
plot_output <- file.path(data_in,"/Plots/04_2_scale_heatmap_GO_terms_bargraph_nonimpute.pdf")

# create a -log10(p_value) column
data$`-log10(p_value)` <- -log10(data$pvalue)
data$`-log10(p_adjust)` <- -log10(data$p_adjust)

# convert gene_ratio to a decimal number
data$gene_ratio <- sapply(data$gene_ratio, function(x) eval(parse(text = x)))

# convert bg_ratio to a decimal number
data$bg_ratio <- sapply(data$bg_ratio, function(x) eval(parse(text = x)))

# Create enrichment factor
data$enrichment_factor <- data$gene_ratio / data$bg_ratio
# Rich factor
data$rich_factor <- data$count / data$bg_ratio

# Filter out GO terms with less than 50 entries
data <- data %>% filter(count >= 20)

# Spit data by cluster
split_list <- data%>% group_split(cluster_no)

# Sort each dataframe by enrichment_factor in descending order
split_list_sort <- lapply(split_list, function(x) {
  x <- arrange(x, p_adjust) %>% #filter for the 10 lowest p_values
  filter(row_number() <= 5) 
  return(x)
})
View(split_list_sort[[2]])

# arrange the category value in the desired order for the plot by sorting the category value according to the enrichment factor in a descending order
split_list_sort <- lapply(split_list_sort, function(x) {
  x$description <- factor(x$description, levels = unique(x$description[order(x$gene_ratio)]), ordered = TRUE)
  return(x)
})

# ceate an emty list to store plots
split_list_plots <- list()

#make a dot plot for each Cluster with enrichment factor as x-axis and the ordered category_values as y-axis
for (df in split_list_sort) {
  
  # define legend
  legend_colors <- c("#b3eebe", "#46bac2", "#371ea3")
  legend_limits <- c(0, 250)
  
  # Wrap long labels
  df$description <- str_wrap(df$description, width = 30)
  
  # Create a bar plot
  bar_plot <- ggplot(df, showCategory=5, aes(gene_ratio, fct_reorder(description, gene_ratio), fill=`-log10(p_adjust)`)) + 
    geom_col() + scale_fill_gradientn(colours= legend_colors, limits = legend_limits, 
                                      guide=guide_colorbar(reverse=TRUE)) + 
    theme_dose(12)+ 
    labs(title = paste("GO CC", df$cluster_no[1]),
         x = "Gene Ratio",  # Add your custom x-axis label here
         fill = "-log10(p_adjust)") + # Custom legend text for color
    ylab(NULL) +
    scale_y_discrete(position="right") 

    
    # Save the plot in list
  split_list_plots[[df$cluster_no[1]]] <- bar_plot
  
  
}
split_list_plots[[1]]

# Combine all plots into one grid with a shared legend
combined_plot <- wrap_plots(split_list_plots) + plot_layout(guides = "collect")
combined_plot


## Save plots --------------------------
pdf(file = plot_output, width = 15, height = 4.5)
combined_plot
dev.off()


