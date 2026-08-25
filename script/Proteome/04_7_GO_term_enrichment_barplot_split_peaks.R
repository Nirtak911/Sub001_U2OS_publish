#load libraries
library(tidyverse)
library(protti)
library(ggplot2)
library(forcats)
library(DOSE)
library(gridExtra)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(purrr)

## Load data ---------------------------
input<- file.path(data_in,"/04_6_Profile_plots_selected_markers_CTRLalone_nonimpute_peaks_split.csv")
data<- read_protti(input)
View(data)

# output data  ---------------------------
plot_output <- file.path(data_in,"/Plots/04_7_scale_split_peaks_GO_terms_bargraph_nonimpute.pdf")

# Create a named list of gene vectors for each custom marker group
gene_lists_by_marker <- data %>%
  filter(!is.na(pg_genes)) %>%
  distinct(hsap_location_marker, pg_genes) %>%
  group_split(hsap_location_marker) %>%
  map(~ .x$pg_genes)

# Name the list elements so clusterProfiler knows which is which
names(gene_lists_by_marker) <- unique(data_prepared$hsap_location_marker)

# Run the comparative cellular component enrichment
go_enrichment_results <- compareCluster(
  geneClusters = gene_lists_by_marker, 
  fun          = "enrichGO", 
  OrgDb        = org.Hs.eg.db, 
  keyType      = "SYMBOL",     # Matches your 'pg_genes' format
  ont          = "CC",         # Cellular Component
  pAdjustMethod = "BH",        # Benjamini-Hochberg FDR control
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Convert results to a standard data frame for saving
go_df <- as.data.frame(go_enrichment_results)

# Generate a comparative dot plot
p_go <- dotplot(
  go_enrichment_results, 
  showCategory = 3,          # Shows top 3 enriched GO terms per marker group
  title        = "GO Cellular Component Enrichment of Subcellular Subpopulations"
) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p_go

# 1. Gather all unique markers from your prepared dataset
unique_markers <- unique(data_prepared$hsap_location_marker)

# Create an empty list to store the final ggplot objects
split_list_plots <- list()

# 2. Run enrichment and plotting one-by-one
for (marker in unique_markers) {
  
  # Extract the genes for the current marker
  marker_genes <- data_prepared %>%
    filter(hsap_location_marker == marker) %>%
    filter(!is.na(pg_genes)) %>%
    distinct(pg_genes) %>%
    pull(pg_genes)
  
  # Skip if the group contains no genes
  if (length(marker_genes) == 0) next
  
  message("Processing: ", marker)
  
  # Run the individual enrichment
  ego <- enrichGO(
    gene          = marker_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "CC", 
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # cutoff = 0.7 is standard. Lower numbers (e.g., 0.5) collapse more terms.
  ego <- clusterProfiler::simplify(ego, 
                                   cutoff = 0.9, 
                                   by = "p.adjust", 
                                   select_fun = min)
  
  # Convert the enrichResult object into a standard data frame
  ego_df <- as.data.frame(ego)
  
  # Check if any significant terms were found; skip if empty
  if (is.null(ego_df) || nrow(ego_df) == 0) {
    message(" -> No significant terms found for: ", marker)
    next
  }
  
  # 3. BASE R COLUMN MAPPING (Immune to package masking & case issues)
  # Find columns matching names regardless of uppercase/lowercase layout
  col_lower <- tolower(colnames(ego_df))
  
  # Directly assign names to columns using their positions
  colnames(ego_df)[col_lower == "p.adjust"]    <- "p_adjust"
  colnames(ego_df)[col_lower == "generatio"]   <- "gene_ratio"
  colnames(ego_df)[col_lower == "bgratio"]     <- "bg_ratio"
  colnames(ego_df)[col_lower == "count"]       <- "count"
  colnames(ego_df)[col_lower == "description"] <- "description"
  
  # Add metadata and transformations
  ego_df <- ego_df %>%
    mutate(
      cluster_no = marker, 
      `-log10(p_adjust)` = -log10(p_adjust)
    )
  
  # Safely parse fraction strings (e.g., "5/100") into decimal numbers
  ego_df$gene_ratio <- sapply(ego_df$gene_ratio, function(x) eval(parse(text = x)))
  ego_df$bg_ratio   <- sapply(ego_df$bg_ratio,   function(x) eval(parse(text = x)))
  
  # Calculate enrichment factor metrics
  ego_df$enrichment_factor <- ego_df$gene_ratio / ego_df$bg_ratio
  ego_df$rich_factor       <- ego_df$count / ego_df$bg_ratio
  
  # 4. Filter and Sort according to your specific display rules
  ego_df_filtered <- ego_df %>% 
    filter(count >= 5) 
  
  if (nrow(ego_df_filtered) == 0) next
  
  # Keep only the top 5 lowest p_adjust values
  ego_df_sorted <- ego_df_filtered %>%
    arrange(p_adjust) %>%
    filter(row_number() <= 5)
  
  # Sort and set factor levels based on gene_ratio magnitude
  ego_df_sorted$description <- factor(
    ego_df_sorted$description, 
    levels = unique(ego_df_sorted$description[order(ego_df_sorted$gene_ratio)]), 
    ordered = TRUE
  )
  
  # 5. Build your exact ggplot layout
  legend_colors <- c("#b3eebe", "#46bac2", "#371ea3")
  max_p_val <- max(ego_df_sorted$`-log10(p_adjust)`, na.rm = TRUE)
  legend_limits <- c(0, max(20, ceiling(max_p_val))) 
  
  # Wrap long labels so they fit nicely on the facet edge
  ego_df_sorted$description <- str_wrap(ego_df_sorted$description, width = 30)
  
  bar_plot <- ggplot(ego_df_sorted, aes(gene_ratio, fct_reorder(description, gene_ratio), fill = `-log10(p_adjust)`)) + 
    geom_col() + 
    scale_fill_gradientn(
      colours = legend_colors, 
      limits = legend_limits, 
      guide = guide_colorbar(reverse = TRUE)
    ) + 
    theme_dose(12) + 
    labs(
      title = paste("GO CC:", marker),
      x = "Gene Ratio",  
      fill = "-log10(p_adjust)"
    ) + 
    ylab(NULL) +
    scale_y_discrete(position = "right") 
  
  # Store the compiled plot object inside our list named by marker
  split_list_plots[[marker]] <- bar_plot
}

# 6. Grid stitching using patchwork
if (length(split_list_plots) > 0) {
  combined_plot <- wrap_plots(split_list_plots) + plot_layout(guides = "collect")
  combined_plot
} else {
  message("No plots generated. Ensure your criteria match the size of your gene pools.")
}

# 6. Define the markers grouped by compartment with increasing Fraction Peaks
target_markers <- c(
  
  # 2. Endoplasmic Reticulum group (increasing fractions)
  "Endoplasmic Reticulum (Frac 3 Peak)",
  "Endoplasmic Reticulum (Frac 4 Peak)",
  
  # 4. Peroxisome group (increasing fractions)
  "Peroxisome (Frac 3 Peak)",
  "Peroxisome (Frac 4 Peak)",
   
  # 3. Lysosome group (increasing fractions)
  "Lysosome (Frac 3 Peak)",
  "Lysosome (Frac 4 Peak)",
  
  # 1. Nucleus group (increasing fractions)
  "Nucleus (Frac 2 Peak)",
  "Nucleus (Frac 5-6 Peak)"
)

# Intersect to drop any that might be missing/empty from the analysis
available_targets <- intersect(target_markers, names(split_list_plots))

if (length(available_targets) > 0) {
  
  # This extracts and forces the exact sequential order defined above
  selected_plots <- split_list_plots[available_targets]
  
  # Combine them. You can optionally set ncol = 2 to make neat side-by-side pairs!
  combined_plot <- wrap_plots(selected_plots, ncol = 2) + 
    plot_layout(guides = "collect")
  
  combined_plot
  
} else {
  message("None of the specified target plots were found in split_list_plots.")
}

## Save plots --------------------------
pdf(file = plot_output, width = 10, height = 10)
combined_plot
dev.off()


