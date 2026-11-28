# Analysis pipeline for "Sequential cellular fractionation strategy for sensitive profiling of subcellular proteome and phosphoproteome dynamics"

This repository contains the R scripts used to analyse the subcellular proteome and phosphoproteome data of the protocol (U2OS cells, anisomycin with or without p38 or ZAK inhibition, six fractions, four biological replicates). The scripts follow the data analysis steps of the Procedure (Step 10) and generate the result tables and the figures of the manuscript.

The pipeline is fully unsupervised: it quantifies changes in the relative distribution of proteins and phosphosites across fractions between conditions, but does not assign proteins to compartments or predict their localization.

**Data availability**
- Input data and result tables: Zenodo, DOI [10.5281/zenodo.24631325](https://doi.org/10.5281/zenodo.24631325)
- Raw MS data and Spectronaut search results: ProteomeXchange / PRIDE, [PXD067608](https://www.ebi.ac.uk/pride/archive/projects/PXD067608)
- Example dataset and vignette (Proteome): `example_data/` and `vignettes/example_workflow.html`, see [Example dataset and vignette](#example-dataset-and-vignette)

---

## Contents

1. [Requirements](#requirements)
2. [Input data](#input-data)
3. [Running the pipeline](#running-the-pipeline)
4. [Example dataset and vignette](#example-dataset-and-vignette)
5. [Pipeline overview](#pipeline-overview)
6. [Statistical approach](#statistical-approach)
7. [Main result tables](#main-result-tables)
8. [Scripts and figures](#scripts-and-figures)
9. [Repository structure](#repository-structure)
10. [Citation and license](#citation-and-license)

---

## Requirements

- **R 4.5.x** (the analysis was run with R 4.5.0). With other R versions (e.g. R 4.6), `renv::restore()` cannot use the recorded package versions; install R 4.5 alongside and select it (RStudio: *Tools → Global Options → General → R version*).
- **Build tools**, because `renv::restore()` compiles packages for which no pre-built binary of the recorded version is available any more:
  - Windows: [Rtools 4.5](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) (check in R with `Sys.which("make")`, which should return a path in the Rtools folder)
  - macOS: Xcode command line tools (`xcode-select --install`)
  - Linux: a C/C++/Fortran compiler toolchain (e.g. `r-base-dev` on Debian/Ubuntu)
- All package versions are recorded in `renv.lock`. To install them into a project library, open the project (`Sub001_U2OS_publish.Rproj`) or start R in the project root and run:

```r
install.packages("renv")   # if renv is not installed yet
renv::restore()
```

Answer `y` to install the packages. If renv offers to install the latest available versions of packages that failed, answer `n` and check the build tools instead: newer versions would no longer match the environment of the analysis. Installing all packages takes about 30–60 minutes. Run `renv::restore()` in a fresh R session (*Session → Restart R*), because packages that are already loaded cannot be replaced on Windows.

Main packages: tidyverse, protti, limma, DAPAR / imp4p, qvalue, pheatmap, clusterProfiler, org.Hs.eg.db, Rtsne, umap, patchwork.

---

## Input data

The scripts expect the following files (not included in the repository, see `.gitignore`):

| File | Description | Zenodo file |
|---|---|---|
| `raw/Proteome/20250221_110126_ProteinGroups_KS_normal (Normal).tsv` | Spectronaut protein group report (directDIA+, fractions searched separately and combined) | `input_Proteome.zip` |
| `raw/Phospho/20250220_171101_PTM Report_normal_KS (Normal).tsv` | Spectronaut PTM site report | `input_Phospho.zip` |
| `raw/Subcell_markers_selected.csv` | Curated subcellular marker proteins (Supplementary Table S1) | `input_annotations.zip` |
| `raw/annotations_hpa_2024-03-05.csv` | DeepLoc 2.0 subcellular localization prediction per protein group (columns `PG.ProteinGroups`, `PG.Genes`, `deeploc_single`), used to annotate heatmaps and dimensionality reduction plots | `input_annotations.zip` |

These files and the main result tables of both pipelines are archived on Zenodo (DOI [10.5281/zenodo.24631325](https://doi.org/10.5281/zenodo.24631325)). Each ZIP file unpacks into the folder structure above; unpack them in the project root. The `results_*.zip` files contain the result tables reported in the manuscript (in `data/Proteome/` and `data/Phospho/`).

The raw MS data and the Spectronaut search results are available from ProteomeXchange via PRIDE (PXD067608).

The folders can be changed without editing the scripts by setting environment variables (see `script/*/00_path_definitions.R`):

| Variable | Default | Content |
|---|---|---|
| `SUB001_RAW_DIR` | `raw` | input reports (`<dir>/Proteome`, `<dir>/Phospho`) and annotation files |
| `SUB001_DATA_DIR` | `data` | intermediate and result tables (`<dir>/Proteome`, `<dir>/Phospho`), plots in `<dir>/<pipeline>/Plots` |

---

## Running the pipeline

All paths are relative to the project root. Every script loads its paths from `00_path_definitions.R` and reads its inputs from files written by earlier scripts, so each script can be run on its own, as long as the scripts before it have been run once.

**Without RStudio (complete pipeline):**

```sh
Rscript run_pipeline.R                 # Proteome and Phospho, all scripts
Rscript run_pipeline.R proteome        # only the Proteome pipeline
Rscript run_pipeline.R phospho core    # only the scripts that produce the result tables
Rscript run_pipeline.R proteome from=03_0   # resume at the first script starting with 03_0
```

`run_pipeline.R` runs each script in its own R session in the order given below, writes the console output of each script to `logs/<pipeline>/<script>.log`, and stops at the first error. After fixing the error, the run can be resumed with `from=<script prefix>`.

**In RStudio:** open `Sub001_U2OS_publish.Rproj` and run the scripts of `script/Proteome` and `script/Phospho` in the order of their numbers. Data tables are shown with `View()` only in interactive sessions.

Runtime: the complete pipeline takes several hours on a standard desktop computer, mostly for imputation (DAPAR) and the limma models. Keep the project on a local disk: on a network drive, reading the reports and the package library is much slower.

---

## Example dataset and vignette

`example_data/input/` contains a small Proteome dataset to try the pipeline without downloading the full data: a subset of the Spectronaut protein group report with all identified marker proteins of `Subcell_markers_selected.csv` (639) plus 500 randomly selected proteins, for all 96 runs (6 fractions × 4 conditions × 4 replicates). The report is stored gzipped (4 MB); `00_path_definitions.R` uses the `.tsv.gz` file when no uncompressed report is present.

The vignette `vignettes/example_workflow.Rmd` runs the Proteome pipeline (`01`–`06`) on this subset, reproduces the marker profile plot of Figure 3c and builds a ranked table of candidate redistributions. Render it from the project root (about 15–30 min):

```r
rmarkdown::render("vignettes/example_workflow.Rmd")
```

The rendered result is included as `vignettes/example_workflow.html`. Results of the example run are written to `example_data/results/` (not tracked). Adjusted p-values and the set of significant redistributions differ from the full dataset, because the Benjamini–Hochberg correction covers far fewer tests and normalisation and imputation are estimated from fewer proteins.

To run the scripts on the example data without the vignette, set the input and output folders first, e.g. in PowerShell:

```powershell
$env:SUB001_RAW_DIR = "example_data/input"; $env:SUB001_DATA_DIR = "example_data/results"
Rscript run_pipeline.R proteome core
```

---

## Pipeline overview

Scripts marked **core** produce the result tables; the other scripts generate quality control plots and figures. The same steps are used for the Proteome (`pg_` columns, protein groups) and the Phospho pipeline (`ptm_` columns, phosphosites).

| Script | Procedure step | Input | Output | Scale |
|---|---|---|---|---|
| `01_data_import_con_filtering_qc.R` (core) | 10.1 | Spectronaut report | `01_data_log2.csv`, `01_genes.csv`, `01_meta_data.csv` | log2 of the Spectronaut quantity; contaminants and entries without gene names removed |
| `02_0_norm_loess_FRsep.R` (core) | 10.2, 10.3 | `01_data_log2.csv` | `02_loess_data_log2.csv` | log2; valid-value filter (≥ 3 of 4 replicates in at least one condition per fraction), cyclic LOESS normalization separately per fraction |
| `02_1_data_qc_after_normalisation.R` | QC | `02_loess_data_log2.csv` | plots, `02_counts.tsv` | log2 |
| `03_0_missingness_protti.R` (core) | 10.4 | `02_loess_data_log2.csv` | `03_0_missingness_protti_data_log2.csv` | log2; classification of missing values per fraction and comparison (complete / MAR / MNAR) |
| `03_0_impute_DAPAR_FRsep.R` (core) | 10.4 | `02_loess_data_log2.csv` | `03_imputeFRsep_data_log2.csv` | log2; imputed candidate values per fraction (see [Imputation](#missing-values-and-imputation)) |
| `03_1_select_dapar_impute_for_limma.R` (core) | 10.4.1 | outputs of the two `03_0` scripts | `03_1_dapar_limma_input_data_log2.csv` | log2; imputed values are used only for MNAR, all other values are the normalized values |
| `03_1_data_qc_after_imputation.R` | QC | `03_imputeFRsep_data_log2.csv` | plots | log2 |
| `04_0_intensity_scaling_not_abs_nonimpute.R` (core) | 10.5 | `02_loess_data_log2.csv` (not imputed) | `04_scale_not_abs_loess_data_log2.csv` (Phospho: `04_0_scale_not_abs_data_log2.csv`) | normalized log2 values back-transformed to linear, missing values set to 0, scaled intensity = proportion of the summed signal across the six fractions (0–1); delta scaled intensity (treatment − control) |
| `04_*` (other) | 10.5, QC | scaled intensities | plots and tables | proportions (0–1) |
| `05_0_limma_differential_expression_dapar_impute_complete_group_omnibus.R` (core) | 10.6, 10.7.3, 10.7.4 | `03_1_dapar_limma_input_data_log2.csv` | `05_0_diff_data_limma_dapar_impute_complete_group_omnibus.csv` | log2; limma contrasts (treatment effect per fraction and fraction-pair interaction terms), BH correction per contrast |
| `05_1_differential_expression_plots_limma_final.R` | 10.6 | `05_0_…csv` | volcano plots, `05_1_diff_expression_sig_count_limma_final.csv` | log2 fold changes |
| `05_2_diff_mobility_limmainteractions.R` (core) | 10.7.3 | `05_0_…csv` | `05_2_diff_mobility_limma_interaction.csv` | interaction results labelled with their loss-to-gain direction (`Frac X->Frac Y`) |
| `06_0_final_movement_table_per_rep.R` (core) | 10.7.1, 10.7.2, 10.7.5 | `05_2_…csv`, scaled intensities from `04_0` | `06_0_diff_and_scaled_mobility_per_rep_final.csv`, `06_0_diff_and_scaled_mobility_top_movements_median.csv` | proportions (0–1) for the scores, log2 for the interaction effect (`diff`) |

---

## Statistical approach

### Missing values and imputation

- Imputation is applied separately within each fraction and is used **only for statistical testing** (Steps 10.6, 10.7). The relative distribution profiles (scaled intensities, Step 10.5) are calculated from the normalized, non-imputed values, with undetected fractions set to zero.
- **MNAR** (`03_0_missingness_protti.R`, `protti::assign_missingness`, `completeness_MNAR = 0.2`, `completeness_MAR = 0.8`): for a comparison of two conditions within a fraction, values are classified as MNAR when a protein or phosphosite is not detected in any replicate of one condition while it is detected in all four replicates of the other condition.
- MNAR values are replaced with a constant low value per sample, the 2.5 % quantile of that sample's intensity distribution (`DAPAR::wrapper.impute.detQuant`, `qval = 0.025`, `factor = 1`).
- The classification is made per comparison and then combined per sample (`03_1_select_dapar_impute_for_limma.R`): if a sample is classified as MNAR in at least one comparison, its imputed values are used in all contrasts in which this sample is included.
- Partially observed conditions are **not imputed** and enter limma with their observed values only. `03_0_impute_DAPAR_FRsep.R` also computes KNN (Proteome) or SLSA (Phospho) values for partially observed conditions, as needed for older-generation instruments (Step 10.4.2), but these values are not used for the Orbitrap Astral data analysed here.
- **Run-to-run variation:** `DAPAR::wrapper.impute.KNN` draws a new random seed for `impute::impute.knn` on every run, and `wrapper.impute.detQuant` computes its quantile after the KNN step. The KNN values themselves are not used, but the MNAR values shift slightly between runs, and with them the limma results. Repeated runs can therefore differ in individual borderline cases around `adj_pval = 0.05`.

### Differential abundance and differential distribution

- Differential abundance within a fraction (Step 10.6) and differential distribution between fractions (Step 10.7) are tested with limma moderated t-tests on the normalized and imputed log2 intensities.
- Differential distribution is tested with interaction terms for all predefined fraction pairs: the interaction tests whether the treatment effect differs between Fraction X and Fraction Y.
- **Multiple-testing correction:** Benjamini–Hochberg, applied **separately within each contrast** (`add_grouped_qvalue_and_bh()` in `05_0`, grouped by `comparison`). The family of hypotheses is all proteins (or phosphosites) tested in a given contrast, so the FDR is controlled per contrast. The list of proteins significant in at least one fraction-pair contrast is therefore **not FDR-controlled as a whole** and should be regarded as a descriptive, prioritized list of candidate redistribution events. For FDR control of the final list, a single BH correction across all protein × fraction-pair tests per treatment can be applied instead (remove the grouping by `comparison`), at the cost of sensitivity.

### Scores (`06_0_final_movement_table_per_rep.R`)

All scores are calculated per protein (or phosphosite), treatment and replicate from the scaled intensities, and reported as the median across replicates.

- **Delta scaled intensity:** ΔS<sub>i</sub> = S<sub>Treatment,i</sub> − S<sub>CTRL,i</sub> for each fraction i. Negative values are losses, positive values gains; all loss-to-gain pairs define the possible trajectories (`Frac X->Frac Y`).
- **Pair mobility score** (`pair_mobility_score`): (ΔS<sub>gain,Y</sub> − ΔS<sub>loss,X</sub>) / 2, the average proportion of signal shifted from Fraction X to Fraction Y (0–1). It is specific to one trajectory and therefore corresponds to the fraction-pair contrast of the interaction test.
- **Mobility score** (`mobility_score`): Σ<sub>i</sub> |ΔS<sub>i</sub>| / n<sub>detected</sub>, where n<sub>detected</sub> is the number of fractions in which the protein is detected under treatment (`count_frac_changing`). It summarizes the overall redistribution of a protein and is shared by all of its trajectories.
- Only trajectories observed in at least three replicates are kept.

**Prioritization (Step 10.7.5):** keep trajectories with a significant interaction (`adj_pval < 0.05`) and sort them by `median_pair_mobility_score`.

---

## Main result tables

Written to `data/Proteome` and `data/Phospho`.

| File | Content |
|---|---|
| `05_0_diff_data_limma_dapar_impute_complete_group_omnibus.csv` | All limma results: `comparison`, `Treatment`, `Contrast_Type` (treatment effect within a fraction or `Differential_Movement` between two fractions), `Fraction_A`, `Fraction_B`, `diff` (log2), `t_statistic`, `pval`, `q_value`, `adj_pval` (BH per contrast) |
| `05_2_diff_mobility_limma_interaction.csv` | Interaction results with direction (`movement`) |
| `06_0_diff_and_scaled_mobility_per_rep_final.csv` | Per replicate and trajectory: loss and gain fraction, delta scaled intensities, pair mobility score, mobility score, interaction statistics |
| `06_0_diff_and_scaled_mobility_top_movements_median.csv` | Per protein, treatment and trajectory: medians across replicates (`median_pair_mobility_score`, `median_mobility_score`) and interaction statistics (`diff`, `pval`, `adj_pval`) |

---

## Scripts and figures

| Figure | Script(s) |
|---|---|
| Fig. 3a, b (identifications per fraction) | `Proteome/02_1`, `Phospho/02_1` |
| Fig. 3c (marker profiles per replicate with centroids) | `Proteome/04_8` |
| Supplementary Fig. S1a (marker profiles) | `Proteome/04_5` |
| Fig. 3d–f (heatmap, cluster profiles, GO CC enrichment) | `Proteome/04_1`, `Proteome/04_3`, `Proteome/04_2` |
| Fig. 4 (t-SNE, PCA) | `Proteome/04_12` |
| Fig. 5a (phosphosite heatmap) | `Phospho/04_1` |
| Fig. 5b–d, Supplementary Fig. S5 (protein and phosphosite intensity and scaled intensity bar plots) | `Proteome/03_2`, `Proteome/04_11`, `Phospho/02_2`, `Phospho/04_3` |
| Supplementary Fig. S1b (split peaks and GO CC) | `Proteome/04_6`, `Proteome/04_7` |
| Supplementary Fig. S1c (marker distribution and enrichment factors) | `Proteome/04_10` |
| Supplementary Fig. S1d (table: intra- vs inter-fraction correlation) | `Proteome/04_4` |
| Supplementary Fig. S3 (marker profiles per treatment) | `Proteome/04_9` |

---

## Repository structure

```text
Sub001_U2OS_publish/
├── run_pipeline.R            # runs all scripts in order without RStudio
├── example_data/input/       # example Proteome subset (report, marker set)
├── vignettes/                # example workflow (Rmd and rendered HTML)
├── script/
│   ├── Proteome/             # protein-level pipeline (scripts 00–06)
│   └── Phospho/              # phosphosite-level pipeline (scripts 00–06)
├── raw/                      # input reports and annotation files (not tracked)
├── data/                     # intermediate and result tables, plots (not tracked)
├── logs/                     # logs written by run_pipeline.R (not tracked)
├── renv.lock                 # package versions
├── renv/                     # renv infrastructure
├── CITATION.cff
├── LICENSE.md
└── Sub001_U2OS_publish.Rproj
```

---

## Citation and license

If you use this workflow, please cite the protocol (Stuber et al., Nature Protocols, in revision), the archived code release on Zenodo (linked under "Related works" of the data record) and, if you use the data, the data record (Zenodo DOI [10.5281/zenodo.24631325](https://doi.org/10.5281/zenodo.24631325)); see `CITATION.cff`. The code is released under the license in `LICENSE.md`.
