# Run the complete analysis pipeline without RStudio.
#
# Usage (from the project root):
#   Rscript run_pipeline.R                 # Proteome and Phospho, all scripts
#   Rscript run_pipeline.R proteome        # only the Proteome pipeline
#   Rscript run_pipeline.R phospho core    # only the scripts that produce the
#                                          # result tables (no QC/figure scripts)
#   Rscript run_pipeline.R proteome from=03_0  # resume, starting at the first
#                                          # script whose name starts with 03_0
#
# Each script is run in its own R session (Rscript), exactly as if it were run
# on its own, so no script depends on objects left over from another one.
# Output of each script is written to logs/<pipeline>/<script>.log once the
# script has finished (while it runs, the log is kept in the local temp folder).
# The pipeline stops at the first script that fails.
#
# Input and output folders can be changed with the environment variables
# SUB001_DATA_DIR and SUB001_RAW_DIR (see script/*/00_path_definitions.R).

args <- tolower(commandArgs(trailingOnly = TRUE))
which_pipeline <- intersect(args, c("proteome", "phospho"))
if (length(which_pipeline) == 0) which_pipeline <- c("proteome", "phospho")
core_only <- "core" %in% args
from_arg <- sub("^from=", "", grep("^from=", args, value = TRUE))

# Scripts in the order they have to be run. core = produces result tables
# used downstream; the others are QC and figure scripts.
pipelines <- list(
  proteome = list(
    dir = "script/Proteome",
    scripts = c(
      "01_data_import_con_filtering_qc.R"                                        = TRUE,
      "02_0_norm_loess_FRsep.R"                                                  = TRUE,
      "02_1_data_qc_after_normalisation.R"                                       = FALSE,
      "03_0_impute_DAPAR_FRsep.R"                                                = TRUE,
      "03_0_missingness_protti.R"                                                = TRUE,
      "03_1_data_qc_after_imputation.R"                                          = FALSE,
      "03_1_select_dapar_impute_for_limma.R"                                     = TRUE,
      "03_2_bargraph_intensity_combined_selected.R"                              = FALSE,
      "04_0_intensity_scaling_not_abs_nonimpute.R"                               = TRUE,
      "04_1_scaled_nonimpute_pheatmap.R"                                         = FALSE,
      "04_2_GO_term_enrichment_barplot_nonimpute.R"                              = FALSE,
      "04_3_scaled_intensity_profile_plots_clusters_nonimpute.R"                 = FALSE,
      "04_4_scaled_intensity_profile_inter_vs_intra_corr.R"                      = FALSE,
      "04_5_scaled_intensity_profile_plots_selected_markers_nonimputed.R"        = FALSE,
      "04_6_scaled_intensity_profile_plots_selected_markers_nonimputed_split_peaks.R" = FALSE,
      "04_7_GO_term_enrichment_barplot_split_peaks.R"                            = FALSE,
      "04_8_scaled_intensity_profile_plots_selected_markers_nonimputed_RepSep_centroids.R" = FALSE,
      "04_9_scaled_intensity_profile_plots_selected_markers_nonimputed_Treatments_centroids.R" = FALSE,
      "04_10_quant_distribution_selected_markers_nonimputed.R"                   = FALSE,
      "04_11_bargraph_intensity_scaled_combined.R"                               = FALSE,
      "04_12_localisation_PCA_UMAP_tSNE_final.R"                                 = FALSE,
      "05_0_limma_differential_expression_dapar_impute_complete_group_omnibus.R" = TRUE,
      "05_1_differential_expression_plots_limma_final.R"                         = FALSE,
      "05_2_diff_mobility_limmainteractions.R"                                   = TRUE,
      "06_0_final_movement_table_per_rep.R"                                      = TRUE
    )
  ),
  phospho = list(
    dir = "script/Phospho",
    scripts = c(
      "01_data_import_con_filtering_qc.R"                                        = TRUE,
      "02_0_norm_loess_FRsep.R"                                                  = TRUE,
      "02_1_data_qc_after_normalisation.R"                                       = FALSE,
      "02_2_bargraph_intensity_Phospho_selection.R"                              = FALSE,
      "03_0_impute_DAPAR_FRsep.R"                                                = TRUE,
      "03_0_missingness_protti.R"                                                = TRUE,
      "03_1_data_qc_after_imputation.R"                                          = FALSE,
      "03_1_select_dapar_impute_for_limma.R"                                     = TRUE,
      "04_0_intensity_scaling_not_abs_nonimpute.R"                               = TRUE,
      "04_1_scaled_pheatmap_nonimpute.R"                                         = FALSE,
      "04_3_bargraph_intensity_scaled_combined.R"                                = FALSE,
      "05_0_limma_differential_expression_dapar_impute_complete_group_omnibus.R" = TRUE,
      "05_1_differential_expression_plots_limma_final.R"                         = FALSE,
      "05_2_diff_mobility_limmainteractions.R"                                   = TRUE,
      "06_0_final_movement_table_per_rep.R"                                      = TRUE
    )
  )
)

rscript <- file.path(R.home("bin"), "Rscript")

for (p in which_pipeline) {
  pl <- pipelines[[p]]
  scripts <- names(pl$scripts)
  if (core_only) scripts <- scripts[pl$scripts]
  if (length(from_arg) == 1) {
    start <- which(startsWith(tolower(scripts), from_arg))[1]
    if (is.na(start)) stop("No script in ", p, " starts with '", from_arg, "'", call. = FALSE)
    scripts <- scripts[start:length(scripts)]
  }

  log_dir <- file.path("logs", p)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  for (s in scripts) {
    path <- file.path(pl$dir, s)
    log_file <- file.path(log_dir, sub("\\.R$", ".log", s))
    message(sprintf("[%s] %s ...", format(Sys.time(), "%H:%M:%S"), path))
    t0 <- Sys.time()
    # write the log to the local temp folder first and copy it afterwards:
    # the scripts print large tables, and writing them line by line to a
    # network drive can slow a script down by hours
    tmp_log <- tempfile(fileext = ".log")
    status <- system2(rscript, c("--no-save", shQuote(path)),
                      stdout = tmp_log, stderr = tmp_log)
    file.copy(tmp_log, log_file, overwrite = TRUE)
    unlink(tmp_log)
    if (!identical(status, 0L)) {
      stop(sprintf("%s failed (exit status %s). See %s", path, status, log_file),
           call. = FALSE)
    }
    message(sprintf("    done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

message("Pipeline finished.")
