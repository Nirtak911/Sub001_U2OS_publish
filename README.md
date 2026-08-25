# Analysis Pipeline for "Sequential cellular fractionation strategy for sensitive profiling of subcellular proteome and phosphoproteome dynamics"

This repository contains the data analysis workflow and R scripts used for processing, quality control, normalization, imputation, differential expression, profile clustering, differential localization, and figure generation for paired subcellular **Proteomics** and **Phosphoproteomics** experiments.

---

## 📁 Repository Structure

```text
Sub001_U2OS_publish/
├── script/
│   ├── Phospho/                                                                 # Phosphoproteomics Processing Pipeline
│   │   ├── 00_path_definitions.R                                                # Project directory and path definitions
│   │   ├── 01_data_import_con_filtering_qc.R                                    # Import, contaminant filtering, and initial QC
│   │   ├── 02_0_norm_loess_FRsep.R                                              # Cyclic LOESS normalization (fraction-separated)
│   │   ├── 02_1_data_qc_after_normalisation.R                                   # Post-normalization QC diagnostics
│   │   ├── 02_2_bargraph_intensity_Phospho_selection.R                           # Intensity distribution and selection bar graphs
│   │   ├── 03_0_impute_DAPAR_FRsep.R                                            # Missing value imputation using DAPAR
│   │   ├── 03_0_missingness_protti.R                                            # Evaluation of missingness patterns via protti
│   │   ├── 03_1_data_qc_after_imputation.R                                      # Post-imputation QC diagnostics
│   │   ├── 03_1_select_dapar_impute_for_limma.R                                 # Expression matrix compilation for statistical testing
│   │   ├── 04_0_intensity_scaling_not_abs_nonimpute.R                          # Non-imputed relative intensity scaling
│   │   ├── 04_1_scaled_pheatmap_nonimpute.R                                     # Hierarchical clustering heatmaps of scaled profiles
│   │   ├── 04_3_bargraph_intensity_scaled_combined.R                           # Combined scaled intensity bar graphs
│   │   ├── 05_0_limma_differential_expression_dapar_impute_complete_group_omnibus.R # Differential expression via limma
│   │   ├── 05_1_differential_expression_plots_limma_final.R                    # Volcano & MA plots for statistical outputs
│   │   ├── 05_2_diff_mobility_limmainteractions.R                               # Interaction testing for differential translocation
│   │   └── 06_0_final_movement_table_per_rep.R                                 # Replicate-level spatial movement matrix compilation
│   │
│   └── Proteome/                                                                # Total Proteomics Processing Pipeline
│       ├── 00_path_definitions.R                                                # Project directory and path definitions
│       ├── 01_data_import_con_filtering_qc.R                                    # Import, contaminant filtering, and initial QC
│       ├── 02_0_norm_loess_FRsep.R                                              # Cyclic LOESS normalization (fraction-separated)
│       ├── 02_1_data_qc_after_normalisation.R                                   # Post-normalization QC diagnostics
│       ├── 03_0_impute_DAPAR_FRsep.R                                            # Missing value imputation using DAPAR
│       ├── 03_0_missingness_protti.R                                            # Evaluation of missingness patterns via protti
│       ├── 03_1_data_qc_after_imputation.R                                      # Post-imputation QC diagnostics
│       ├── 03_1_select_dapar_impute_for_limma.R                                 # Expression matrix compilation for statistical testing
│       ├── 03_2_bargraph_intensity_combined_selected.R                          # Selected marker intensity distribution plots
│       ├── 04_0_intensity_scaling_not_abs_nonimpute.R                          # Non-imputed relative intensity scaling
│       ├── 04_1_scaled_nonimpute_pheatmap.R                                     # Hierarchical clustering heatmaps
│       ├── 04_2_GO_term_enrichment_barplot_nonimpute.R                          # Gene Ontology (GO) term enrichment analysis
│       ├── 04_3_scaled_intensity_profile_plots_clusters_nonimpute.R             # Subcellular profile cluster trajectory plots
│       ├── 04_4_scaled_intensity_profile_inter_vs_intra_corr.R                 # Inter- vs. intra-treatment correlation analysis
│       ├── 04_5_scaled_intensity_profile_plots_selected_markers_nonimputed.R    # Organelle marker profile plots
│       ├── 04_6_scaled_intensity_profile_plots_selected_markers_nonimputed_split_peaks.R # Marker profiles across bimodal fraction peaks
│       ├── 04_7_GO_term_enrichment_barplot_split_peaks.R                        # GO term enrichment on split peak organelle clusters
│       ├── 04_8_scaled_intensity_profile_plots_selected_markers_nonimputed_RepSep_centroids.R # Centroid trajectories across replicates
│       ├── 04_9_scaled_intensity_profile_plots_selected_markers_nonimputed_Treatments_centroids.R # Centroid trajectories across conditions
│       ├── 04_10_quant_distribution_selected_markers_nonimputed.R              # Quantitative distribution across marker sets
│       ├── 04_11_bargraph_intensity_scaled_combined.R                           # Combined scaled intensity bar graphs
│       ├── 04_12_localisation_PCA_UMAP_tSNE_final.R                            # Dimensionality reduction (PCA, UMAP, t-SNE) plots
│       ├── 05_0_limma_differential_expression_dapar_impute_complete_group_omnibus.R # Differential expression via limma
│       ├── 05_1_differential_expression_plots_limma_final.R                    # Volcano & MA plots for statistical outputs
│       ├── 05_2_diff_mobility_limmainteractions.R                               # Differential translocation interaction testing
│       ├── 06_0_final_movement_table_per_rep.R                                 # Replicate-level spatial movement matrix compilation
│       └── 07_0_scaled_permutation_fast_annotated.R                            # Permutation testing for spatial shifts
│
├── .gitignore
├── CITATION.cff
├── LICENSE
├── README.md
└── Sub001_U2OS_publish.Rproj
