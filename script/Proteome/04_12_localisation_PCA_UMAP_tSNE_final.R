#load libraries
library(tidyverse)
library(protti)
library(umap)
library(Rtsne)
library(plotly)
library(viridis)
library(ggprism)
library(plotly)

## Load data---------------------
input <- file.path(data_in,"04_scale_not_abs_imputeFRsep_data_log2.csv")
data <- read_csv(input)
input2 <- "raw/annotations_hpa_2024-03-05.csv"
deeploc <- read_csv(input2)


## output data------------------
plot_output <- file.path(data_in,"/Plots/04_12_PCA_UMAP_tSNE_localisation.pdf")


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


## UMAPS with max scaled intensity Fraction approach---------------------
# only keep the first pg_protein_groups
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, "^;", "")
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, ";.*", "")

# Rename columns if necessary
colnames(deeploc)[colnames(deeploc) == "PG.ProteinGroups"] <- "pg_protein_groups"
# Add deeploc data
data_loc <- data %>% left_join(deeploc, by = "pg_protein_groups") %>%
  filter(treatment == "CTRL")
View(data_loc)

# Find the Fraction with highest scaled_intensity for each pg_protein_groups
max_fraction <- data_loc %>%
  group_by(pg_protein_groups) %>%
  filter(scaled_intensity == max(scaled_intensity)) %>%
  select(pg_protein_groups, fraction, deeploc_single)
View(max_fraction)


# Convert data to wide format
data_wide <- data_loc %>% 
  pivot_wider(id_cols=c(pg_protein_groups),
              names_from=sample,
              values_from=scaled_intensity) %>%
  drop_na() %>% 
  column_to_rownames("pg_protein_groups")
View(data_wide)

# --- 1. RUN PCA (2D) ---
# scale. = TRUE is highly recommended if data isn't already z-scored/scaled
pca_res <- prcomp(data_wide, center = TRUE, scale. = TRUE)

# Calculate variance explained for axis labels
var_explained <- (pca_res$sdev)^2 / sum((pca_res$sdev)^2) * 100
pc1_label <- paste0("PC1 (", round(var_explained[1], 1), "%)")
pc2_label <- paste0("PC2 (", round(var_explained[2], 1), "%)")

# Convert PCA scores to a tibble and pull in protein IDs from rownames
pca_2D_layout <- as_tibble(pca_res$x) %>%
  mutate(pg_protein_groups = rownames(data_wide))


# --- 2. JOIN METADATA & BUILD INDIVIDUAL PLOTS ---
t_pca <- pca_2D_layout %>% 
  left_join(max_fraction, by = "pg_protein_groups")
View(t_pca)

# Plot 1: PCA Fraction
p_pca_1 <- t_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = fraction), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular PCA Map",
       subtitle = "Fraction Assignment",
       x = pc1_label,
       y = pc2_label,
       color = "Fraction") +
  theme_prism()

# Plot 2: PCA Localisation
p_pca_2 <- t_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = deeploc_single), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular PCA Map",
       subtitle = "DeepLoc Localisation",
       x = pc1_label,
       y = pc2_label,
       color = "DeepLoc") +
  theme_prism()


# --- 3. COMBINE PLOTS SIDE-BY-SIDE ---
combined_pca_plot <- p_pca_1 + p_pca_2
combined_pca_plot <- combined_pca_plot + plot_layout(guides = 'collect') 


# Display the combined layout
print(combined_pca_plot)

# --- 1. RUN UMAP (2D) WITH SPECIFIED SETTINGS ---
set.seed(42) # For reproducibility
umap_2D <- umap(data_wide, n_neighbors = 20, min_dist = 0.1)

# Convert UMAP results to a tibble and rename coordinates for clarity
umap_2D_layout <- as_tibble(umap_2D$layout, rownames = "pg_protein_groups") %>% 
  rename(UMAP1 = V1, UMAP2 = V2)


# --- 2. JOIN METADATA & BUILD INDIVIDUAL PLOTS ---
t1 <- umap_2D_layout %>% 
  left_join(max_fraction, by = "pg_protein_groups") 
View(t1)

# Plot 1: UMAP Fraction
p1 <- t1 %>% 
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = fraction), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular UMAP Map",
       subtitle = "Fraction Assignment",
       color = "Fraction") +
  theme_prism()

# Plot 2: UMAP Localisation
p2 <- t1 %>% 
  ggplot(aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = deeploc_single), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular UMAP Map",
       subtitle = "DeepLoc Localisation",
       color = "DeepLoc") +
  theme_prism()


# --- 3. COMBINE PLOTS SIDE-BY-SIDE ---
combined_umap_plot <- p1 + p2
combined_umap_plot <- combined_umap_plot + plot_layout(guides = 'collect') 

# Display the combined layout
print(combined_umap_plot)

# --- PRE-FLIGHT CHECK: REMOVE DUPLICATES FOR t-SNE ---
# Convert to matrix for duplicate checking
matrix_wide <- as.matrix(data_wide)
duplicate_mask <- duplicated(matrix_wide)

# Optional: print a message to console to know how many were dropped
message(paste("Removing", sum(duplicate_mask), "duplicate protein profiles for t-SNE calculation."))

# Subset the matrix and original rownames to keep them synchronized
matrix_unique <- matrix_wide[!duplicate_mask, ]
rownames_unique <- rownames(data_wide)[!duplicate_mask]


# --- 1. RUN t-SNE (2D) ---
set.seed(42) # For reproducibility
# Pass the unique matrix instead of the full data_wide
tsne_2D <- Rtsne(
  matrix_unique, 
  dims = 2, 
  perplexity = 80,   # Increased from 30 to fix stringiness
  max_iter = 1500,   # Increased to allow better convergence
  check_duplicates = FALSE
)

# Convert to tibble using the unique rownames vector
tsne_2D_layout <- as_tibble(tsne_2D$Y) %>% 
  mutate(pg_protein_groups = rownames_unique) %>% 
  rename(tSNE1 = V1, tSNE2 = V2)


# --- 2. 2D PLOTS & TABLES ---
t_tsne <- tsne_2D_layout %>% 
  left_join(max_fraction, by = "pg_protein_groups") 
View(t_tsne)

# Plot 1: t-SNE Fraction
p_tsne_1 <- t_tsne %>% 
  ggplot(aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = fraction), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular t-SNE Map",
       subtitle = "Fraction Assignment",
       color = "Fraction") +
  theme_prism()
p_tsne_1

# Plot 2: t-SNE Localisation
p_tsne_2 <- t_tsne %>% 
  ggplot(aes(x = tSNE1, y = tSNE2)) +
  geom_point(aes(color = deeploc_single), alpha = 0.6, size = 0.8) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Subcellular t-SNE Map",
       subtitle = "DeepLoc Localisation",
       color = "DeepLoc") +
  theme_prism()
p_tsne_2

# 2. Combine them side-by-side using the '+' operator from patchwork
combined_tsne_plot <- p_tsne_1 + p_tsne_2

# 3. Add a super-title and shared legend (optional, but clean)
combined_tsne_plot <- combined_tsne_plot + 
  # plot_annotation(
  #   title = 'Subcellular Proteomics Baseline t-SNE Comparison',
  #   subtitle = 'Comparing experimental fractionation and deep learning predictions',
  #   tag_levels = 'A' # Adds tags (A, B) to the panels
  # ) +
  plot_layout(guides = 'collect') # Optional: collects legends into one (requires legend to match)

# 4. Display the plot
print(combined_tsne_plot)

## Save plots --------------------------
pdf(file = plot_output, width = 10, height = 4.5)
combined_pca_plot
combined_umap_plot
combined_tsne_plot
dev.off()
