# Define common paths for all scripts of the Phospho pipeline.
# This file is sourced at the top of every script, so each script can be run
# on its own (RStudio or Rscript). Paths are relative to the project root
# (open Sub001_U2OS_publish.Rproj, or run Rscript from the project root).
#
# The base folders can be changed without editing this file by setting
# environment variables before running the scripts, e.g. to run the pipeline
# on the example dataset or to write results to another folder:
#   SUB001_DATA_DIR : folder for intermediate and result files (default "data")
#   SUB001_RAW_DIR  : folder with input reports and annotations (default "raw")

data_root <- Sys.getenv("SUB001_DATA_DIR", unset = "data")
raw_root  <- Sys.getenv("SUB001_RAW_DIR",  unset = "raw")

# Data in and output (intermediate and result tables, plots in data_in/Plots)
data_in <- file.path(data_root, "Phospho")

# raw data (Spectronaut PTM report)
raw_data <- file.path(raw_root, "Phospho")

# annotation files shared by both pipelines (marker sets, HPA annotation)
annot_dir <- raw_root

# script path
script <- "script/Phospho"

# create output folders if they do not exist yet
dir.create(file.path(data_in, "Plots"), recursive = TRUE, showWarnings = FALSE)

