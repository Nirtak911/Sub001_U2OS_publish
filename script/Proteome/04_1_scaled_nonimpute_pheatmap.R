#load libraries
library(tidyverse)
library(protti)
library(pheatmap)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(AnnotationDbi)
library(scales)

## Load data ---------------------------
input<- file.path(data_in,"04_scale_not_abs_loess_data_log2.csv")
data<- read_csv(input)
genes <- file.path(data_in, "01_genes.csv")
genes <- read_csv(genes)
input2 <- "raw/annotations_hpa_2024-03-05.csv"
deeploc <- read_csv(input2)

## output data  ---------------------------
plot_output <- file.path(data_in,"/Plots/04_1_pheatmap_scale_data_log2.pdf")
plot_output_png <- file.path(data_in,"/Plots/04_1_pheatmap_scale_data_log2.png")
gsea_output <- file.path(data_in,"/4_1_gsea_cluster_pheatmap_scale_data_log2.csv")
cluster_output <- file.path(data_in,"/4_1_Cluster_pheatmap_scale_data_log2.csv.csv")

# Filter input data
data <- data %>% 
  filter(cell_line=="U2OS")

# Heatmap -----------------------------------------------------------------
pheatmap_input <- data %>%
  #filter(PG.ProteinGroups %in% tox_significant_all$PG.ProteinGroups) %>% 
  pivot_wider(
    id_cols=c(pg_protein_groups),
    names_from = sample, 
    values_from = scaled_intensity) %>% 
  column_to_rownames("pg_protein_groups")

#replace NA with 0
pheatmap_input <- replace(pheatmap_input, is.na(pheatmap_input), 0)

num_rows <- nrow(pheatmap_input)
num_rows

# Add annotation to the heatmap
# row annotation
pheatmap_row<- data %>%
  distinct(pg_protein_groups,.keep_all = T) %>%
  filter(pg_protein_groups %in% rownames(pheatmap_input)) %>% 
  left_join(deeploc,by=c("pg_protein_groups"="PG.ProteinGroups"))%>%
  mutate(Subcell_location=deeploc_single) %>% 
  column_to_rownames("pg_protein_groups")

# col annotation
pheatmap_col <- data %>%
  dplyr::select(sample,fraction,r_condition,r_file_name, cell_line, treatment, r_replicate) %>%
  distinct(sample,.keep_all = T) %>%
  column_to_rownames("sample")


# define colors for the annotations
annotation_colors=list( treatment = c("CTRL" = "#abdbe3", "Anis" = "#063970", "Zaki" = "#f46d43", "p38" = "#d73027"),
                        fraction = c("Frac1" = "#E7872B", "Frac2" = "#F2B341", 
                                     "Frac3" = "#F0E442", "Frac4" = "#009E73", 
                                     "Frac5" = "#56B4E9", "Frac6" = "#0072B2"),
                        #"total" = "#999999"),
                        Subcell_location = c("NA"  = "lightgrey",                    
                          "Cytoplasm" = "#33A02C",            
                          "Cell membrane"  = "#FB9A99",       
                          "Mitochondrion" = "#E31A1C",         
                          "Endoplasmic reticulum" = "#FF7F00",
                          "Golgi apparatus" = "#FDBF6F",
                          "Peroxisome"   = "#CAB2D6",
                          "Lysosome/Vacuole"  = "#6A3D9A", 
                          "Nucleus" = "#1F78B4",   
                          "Extracellular" = "#A6CEE3"),
                        row_cluster=setNames(hcl.colors(7, palette = "viridis"),c(1:7)))

# define number of clusters
number_cluster_columns<- NA
number_cluster_proteins<-6

p<-pheatmap(pheatmap_input,         
            annotation_row = pheatmap_row %>% dplyr::select(Subcell_location),
            annotation_col = pheatmap_col %>% dplyr::select(fraction, treatment),
            show_colnames = F,
            show_rownames = F,
            #scale = "row",
            clustering_distance_cols = "euclidean",
            clustering_distance_rows = "euclidean",
            cutree_rows = number_cluster_proteins,
            cutree_cols = number_cluster_columns,
            main= paste(pheatmap_col$cell_line[1], " n = ", num_rows),
            ledgend = T,
            fontsize_row = 6,
            border_color = "darkgrey",
            color= RColorBrewer::brewer.pal("YlOrBr", n = 9),
            annotation_colors = annotation_colors,
            filename= NA,
            width = 10,
            height = 8)

p

# # Extract row cluster information for heatmap
# clusters <- cutree(p$tree_row, k=number_cluster_proteins)[p$tree_row[["order"]]]
# 
# # Combine cluster information with pheatmap_col
# annot_row <- data.frame(pg_protein_groups = names(clusters),
#                         row_cluster = as.factor(clusters))

# 1. Get the clusters (this returns a named vector)
cluster_vector <- cutree(p$tree_row, k = number_cluster_proteins)

# 2. Convert to a clear dataframe with ID and Cluster
annot_row <- data.frame(
  pg_protein_groups = names(cluster_vector),
  row_cluster = as.factor(cluster_vector)
)

# 3. Join this to your data (the join will match by ID, preventing shifts)
res <- data %>% 
  left_join(annot_row, by = "pg_protein_groups") %>%
  left_join(genes,by="pg_protein_groups") 
View(res)

pheatmap_row <- pheatmap_row %>% rownames_to_column(var = "pg_protein_groups") %>%
  left_join(annot_row, by = "pg_protein_groups") %>%
  column_to_rownames("pg_protein_groups")


# Heatmap with annotated clusters
p2<-pheatmap(pheatmap_input,
            annotation_row = pheatmap_row %>% dplyr::select(Subcell_location, row_cluster),
            annotation_col = pheatmap_col %>% dplyr::select(fraction, treatment),
            show_colnames = F,
            show_rownames = F,
            #scale = "row",
            clustering_distance_cols = "euclidean",
            clustering_distance_rows = "euclidean",
            cutree_rows = number_cluster_proteins,
            cutree_cols = number_cluster_columns,
            main= paste(pheatmap_col$cell_line[1], " n = ", num_rows),
            ledgend = T,
            fontsize_row = 6,
            border_color = "darkgrey",
            color= RColorBrewer::brewer.pal("YlOrBr", n = 9),
            annotation_colors = annotation_colors,
            filename= NA)
            #width = 10,
            #height = 8)

p2

pheatmap_rows<-tibble(pg_protein_groups=p$tree_row$labels,
row_cluster=cutree(p$tree_row,k=number_cluster_proteins))


# define universe for GO
universe <- res$pg_genes

# define gsea_test function
gse_test <-function(gene_names)   {
  enrichGO(  gene ={{gene_names}},
             ont ="CC",
             keyType = "SYMBOL",
             minGSSize = 10,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2,
             OrgDb = "org.Hs.eg.db",
             pAdjustMethod = "BH",
             universe = universe)}


# save all clusters and their corresponding unique pg_genes into a list
clusters<-list(
  cluster1_genes= res %>% filter(row_cluster==1) %>% pull(pg_genes),
  cluster2_genes= res %>% filter(row_cluster==2) %>% pull(pg_genes),
  cluster3_genes= res %>% filter(row_cluster==3) %>% pull(pg_genes),
  cluster4_genes= res %>% filter(row_cluster==4) %>% pull(pg_genes),
  cluster5_genes= res %>% filter(row_cluster==5) %>% pull(pg_genes),
  cluster6_genes= res %>% filter(row_cluster==6) %>% pull(pg_genes))
  #cluster7_genes= res %>% filter(row_cluster==7) %>% pull(pg_genes))

# run gsea_test on all clusters
clusters_gsea <- lapply(clusters,gse_test)
clusters_gsea[[1]]

# # define simplify function
# simplify_function <- function (data) {
#   simplify(data, cutoff=0.7, by="p.adjust", select_fun=min)}
# 
# 
# ## use simplify to remove redundant terms 
# clusters_gsea <- lapply(clusters_gsea, simplify_function)



# Loop through clusters_gsea and extract the top 5 GO terms for each cluster
gse_comb<-bind_rows(list("Cluster 1"= clusters_gsea$cluster1_genes@result,
                         "Cluster 2"= clusters_gsea$cluster2_genes@result,
                         "Cluster 3"= clusters_gsea$cluster3_genes@result,
                         "Cluster 4"= clusters_gsea$cluster4_genes@result,
                         "Cluster 5"= clusters_gsea$cluster5_genes@result,
                         "Cluster 6"= clusters_gsea$cluster6_genes@result),
                    # "Cluster 7"= clusters_gsea$cluster7_genes@result),
                    .id="cluster_no" )
  
View(gse_comb)

# add the information for which cell_line this is the gsea analysis
gse_comb$cell_line<-pheatmap_col$cell_line[1]


# Save gsea as csv
write.csv(res, cluster_output)
write.csv(gse_comb, gsea_output)

## Save plots --------------------------
pdf(file = plot_output, width = 6, height = 8)
p2
dev.off()

png(file = plot_output_png, width = 6, height = 8, units = "in", res = 600)
p2
dev.off()
