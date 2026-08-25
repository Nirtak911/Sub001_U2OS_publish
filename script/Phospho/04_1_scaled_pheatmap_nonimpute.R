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
input<- file.path(data_in,"/04_0_scale_not_abs_data_log2.csv")
data<- read_csv(input)
genes <- file.path(data_in, "01_genes.csv")
genes <- read_csv(genes)
input2 <- "raw/annotations_hpa_2024-03-05.csv"
deeploc <- read_csv(input2)
View(data)

## output data  ---------------------------
plot_output <- file.path(data_in,"/Plots/04_1_pheatmap_scale_data_log2.pdf")
plot_output_png <- file.path(data_in,"/Plots/04_1_pheatmap_scale_data_log2.png")
gsea_output <- file.path(data_in,"/4_1_gsea_cluster_pheatmap_scale_data_log2.csv")
cluster_output <- file.path(data_in,"/4_1_Cluster_pheatmap_scale_data_log2.csv")

# add ptm_protein_id
data$ptm_protein_id <- map_chr(str_split(data$ptm_collapse_key, "_"), 1)

# Filter input data
data <- data %>% 
  filter(cell_line=="U2OS")

# Heatmap -----------------------------------------------------------------
pheatmap_input <- data %>%
  #filter(PG.ProteinGroups %in% tox_significant_all$PG.ProteinGroups) %>% 
  pivot_wider(
    id_cols=c(ptm_collapse_key),
    names_from = sample, 
    values_from = scaled_intensity) %>% 
  column_to_rownames("ptm_collapse_key")

#replace NA with 0
pheatmap_input <- replace(pheatmap_input, is.na(pheatmap_input), 0)

num_rows <- nrow(pheatmap_input)
num_rows

# Add annotation to the heatmap
# row annotation
pheatmap_row<- data %>%
  distinct(ptm_collapse_key,.keep_all = T) %>%
  filter(ptm_collapse_key %in% rownames(pheatmap_input)) %>% 
  left_join(deeploc,by=c("ptm_protein_id"="PG.ProteinGroups"))%>%
  mutate(Subcell_location=deeploc_single) %>% 
  column_to_rownames("ptm_collapse_key")

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
                        row_cluster=setNames(hcl.colors(6, palette = "viridis"),c(1:6)))

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

# Extract row cluster information for heatmap
clusters <- cutree(p$tree_row, k=number_cluster_proteins)[p$tree_row[["order"]]]

# Combine cluster information with pheatmap_col
annot_row <- data.frame(ptm_collapse_key = names(clusters),
                        row_cluster = as.factor(clusters))

pheatmap_row <- pheatmap_row %>% rownames_to_column(var = "ptm_collapse_key") %>%
  left_join(annot_row, by = "ptm_collapse_key") %>%
  column_to_rownames("ptm_collapse_key")


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

pheatmap_rows<-tibble(ptm_collapse_key=p$tree_row$labels,
row_cluster=cutree(p$tree_row,k=number_cluster_proteins))

# add cluster information to the input data
res <-data %>% left_join(annot_row,by="ptm_collapse_key") %>%
  left_join(genes,by="ptm_collapse_key") 

res$pg_genes <- map_chr(str_split(res$gene_PTM_key, "_"), 1)
View(res)

# define gsea_test function
gse_test <-function(gene_names)   {
  enrichGO(  gene ={{gene_names}},
             ont ="CC",
             keyType = "SYMBOL",
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             OrgDb = "org.Hs.eg.db",
             pAdjustMethod = "BH")}


# save all clusters and their corresponding unique pg_genes into a list
clusters<-list(
  cluster1_genes= res %>% filter(row_cluster==1) %>% pull(pg_genes),
  cluster2_genes= res %>% filter(row_cluster==2) %>% pull(pg_genes),
  cluster3_genes= res %>% filter(row_cluster==3) %>% pull(pg_genes),
  cluster4_genes= res %>% filter(row_cluster==4) %>% pull(pg_genes),
  cluster5_genes= res %>% filter(row_cluster==5) %>% pull(pg_genes),
  cluster6_genes= res %>% filter(row_cluster==6) %>% pull(pg_genes))

# run gsea_test on all clusters
clusters_gsea <- lapply(clusters,gse_test)

# Loop through clusters_gsea and extract the top 5 GO terms for each cluster
gse_comb<-bind_rows(list("Cluster 1"= clusters_gsea$cluster1_genes@result,
                         "Cluster 2"= clusters_gsea$cluster2_genes@result,
                         "Cluster 3"= clusters_gsea$cluster3_genes@result,
                         "Cluster 4"= clusters_gsea$cluster4_genes@result,
                         "Cluster 5"= clusters_gsea$cluster5_genes@result,
                         "Cluster 6"= clusters_gsea$cluster6_genes@result),
                    .id="cluster_no" ) %>%
  group_by(cluster_no) %>%
  slice_min(pvalue,n=10) %>%
  unite(col="GO",ID,Description, sep=" ",remove=F)
View(gse_comb)

# add the information for which cell_line this is the gsea analysis
gse_comb$cell_line<-pheatmap_col$cell_line[1]


# Save gsea as csv
write.csv(res, cluster_output)
write.csv(gse_comb, gsea_output)

## Save plots --------------------------
pdf(file = plot_output, width = 5, height = 8)
p
dev.off()

png(file = plot_output_png, width = 5, height = 8, units = "in", res = 600)
p
dev.off()
