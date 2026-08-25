#load libraries
library(tidyverse)
library(protti)
library(limma)
library(qvalue) # Need to load qvalue explicitly

## Load data ---------------------------
input <- file.path(data_in,"/03_1_dapar_limma_input_data_log2.csv")
meta <- file.path(data_in,"01_meta_data.csv")
genes <- file.path(data_in,"01_genes.csv")
data <- read_protti(input)
meta <- read_csv(meta)
genes <- read_csv(genes)
View(data) # Commenting out for execution flow

## output data  ---------------------------
diff_output <- file.path(data_in,"/05_0_diff_data_limma_dapar_impute_complete_group_omnibus.csv")

# cleanup r_condition
data <- data %>% dplyr::mutate(
  r_condition = stringr::str_remove(r_condition, "^U2OS_")) %>%
  dplyr::mutate(
    parts = stringr::str_split(r_condition, pattern = "_", n = 2, simplify = TRUE),
    r_condition = paste(parts[, 2], parts[, 1], sep = "_") # Swap to Frac_Treat
  ) %>%
  dplyr::select(-parts)

# 1. Create Input Data Matrix (Intensity)
moderated_t_test_input <- data %>%
  dplyr::distinct(ptm_collapse_key, sample, final_intensity_log2) %>%
  tidyr::drop_na(final_intensity_log2) %>%
  dplyr::arrange(sample) %>%
  tidyr::pivot_wider(names_from = sample, values_from = final_intensity_log2)

moderated_t_test_input <- column_to_rownames(moderated_t_test_input, var = "ptm_collapse_key") %>%
  as.matrix()

# 2. Create Map and Block Factor (Ensuring alignment)
moderated_t_test_map <- data %>%
  dplyr::distinct(sample, r_condition, r_replicate) %>%
  dplyr::arrange(sample)

# 3. Align the intensity matrix columns (in case samples were dropped by drop_na)
sample_order <- moderated_t_test_map$sample
moderated_t_test_input <- moderated_t_test_input[, sample_order]

# Create Design Matrix (Cell Means Model)
moderated_t_test_design <- stats::model.matrix(~0 + factor(moderated_t_test_map$r_condition))

# Rename columns to remove the prefix
colnames(moderated_t_test_design) <- str_replace_all(colnames(moderated_t_test_design), "factor\\(moderated_t_test_map\\$r_condition\\)", "")

# Create Replicate Block Factor (from the aligned map)
rep <- factor(moderated_t_test_map$r_replicate)

# estimate intra-replicate correlation (assumes fractions from same replicate are correlated)
corfit <- duplicateCorrelation(moderated_t_test_input, moderated_t_test_design, block=rep)
# corfit$consensus.correlation  # useful to inspect

# Fit linear model
limma_fit <- limma::lmFit(moderated_t_test_input, moderated_t_test_design,  block=rep, correlation=corfit$consensus.correlation)

# --- Contrast Generation ---

# 1. Get Fraction and Treatment names from the design matrix columns
all_conditions <- colnames(moderated_t_test_design)
parsed_conditions <- str_split(all_conditions, "_", simplify = TRUE)

fractions <- unique(parsed_conditions[, 1]) # e.g. Frac1, Frac2, ...
treatments <- unique(parsed_conditions[, 2]) # e.g. Anis, CTRL, p38, Zaki
treatments_non_ctrl <- treatments[treatments != "CTRL"]

all_contrasts <- list()

# --- 1. Treatment Effect within Each Fraction (Direct Effect) ---
for (frac in fractions) {
  for (trt in treatments_non_ctrl) {
    contrast_name <- paste0(trt, "_in_", frac)
    term_treat <- paste0(frac, "_", trt)
    term_ctrl <- paste0(frac, "_CTRL")
    all_contrasts[[contrast_name]] <- paste(term_treat, "-", term_ctrl)
  }
}

# --- 2. Differential Treatment Effect / Movement (Interaction Contrast) ---
frac_pairs <- combn(fractions, 2, simplify = FALSE)

for (trt in treatments_non_ctrl) {
  for (pair in frac_pairs) {
    frac_A <- pair[1]
    frac_B <- pair[2]
    
    # Contrast: (FracA_Treat - FracA_CTRL) - (FracB_Treat - FracB_CTRL)
    contrast_name <- paste0(trt, "_diff_", frac_A, "vs", frac_B)
    
    # The algebraic expression for differential movement:
    contrast_formula <- paste0("(", frac_A, "_", trt, " - ", frac_A, "_CTRL) - (",
                               frac_B, "_", trt, " - ", frac_B, "_CTRL)")
    
    all_contrasts[[contrast_name]] <- contrast_formula
  }
}

# 3. Omnibus Treatment Effect (Independent of Fraction)
# Tests the average effect of treatment (Trt - CTRL) across all fractions.
for (trt in treatments_non_ctrl) {
  contrast_formula <- ""
  
  for (frac in fractions) {
    term_treat <- paste0(frac, "_", trt)
    term_ctrl <- paste0(frac, "_CTRL")
    
    # Add positive treatment term and negative control term
    contrast_formula <- paste0(contrast_formula, " + ", term_treat)
    contrast_formula <- paste0(contrast_formula, " - ", term_ctrl)
  }
  
  # Remove the leading ' + ' from the final string
  contrast_formula <- sub("^\\s\\+\\s", "", contrast_formula)
  
  contrast_name <- paste0(trt, "_Omnibus")
  
  # Add to the list
  all_contrasts[[contrast_name]] <- contrast_formula
}


# --- Apply Contrasts ---
contrast_matrix <- do.call(makeContrasts,
                           c(all_contrasts,
                             list(levels = moderated_t_test_design)))

# Fit contrasts
limma_fit2 <- contrasts.fit(limma_fit, contrast_matrix)
limma_fit3 <- eBayes(limma_fit2)

# --- Result Extraction ---

# Get the names of the contrasts you want to process
contrast_names <- colnames(contrast_matrix)

# Extract results for each comparison and combine into one dataframe
limma_result <- purrr::map_dfr(contrast_names, function(cont_name) {
  
  # Use topTable to extract results for the current contrast
  topTable(limma_fit3,
           coef = cont_name,
           number = Inf,
           adjust.method = "BH") %>%
    
    # Add the name of the comparison as a new column
    tibble::add_column(comparison = cont_name) %>%
    
    # Move the row names (protein groups) to a column
    tibble::rownames_to_column("ptm_collapse_key")
})

# Rename and select final columns
limma_result <- limma_result %>%
  dplyr::rename(
    diff = "logFC",
    t_statistic = "t",
    avg_abundance = "AveExpr",
    pval = "P.Value",
    adj_pval_limma = "adj.P.Val"
  ) %>%
  # Separate the components of the 'comparison' column for easier filtering/analysis
  tidyr::separate(
    col = comparison,
    into = c("Treatment_Type", "Fraction_Comparison"),
    sep = "_diff_|_in_|_Omnibus", # Split by all three types
    remove = FALSE,
    extra = "merge"
  )

# Separate the Fraction_Comparison further (for 'diff' results only)
limma_result <- limma_result %>%
  dplyr::mutate(
    # Identify contrast type based on suffix
    Contrast_Type = case_when(
      grepl("_diff_", comparison) ~ "Differential_Movement",
      grepl("_Omnibus", comparison) ~ "Omnibus_Treatment_Effect",
      TRUE ~ "Direct_Effect"
    ),
    
    # For differential movement, separate the two fractions
    Fraction_A = ifelse(Contrast_Type == "Differential_Movement",
                        stringr::str_extract(Fraction_Comparison, "^Frac\\d+"), NA),
    Fraction_B = ifelse(Contrast_Type == "Differential_Movement",
                        stringr::str_extract(Fraction_Comparison, "Frac\\d+$"), NA),
    
    # Final Treatment name
    Treatment = stringr::str_extract(comparison, "^(Anis|p38|Zaki)")
  ) %>%
  dplyr::select(
    ptm_collapse_key, comparison, Treatment, Contrast_Type,
    Fraction_A, Fraction_B, diff, avg_abundance, t_statistic,
    pval, adj_pval_limma, B
  )


# Limma ----------------------------------------------------------------------
# Function to add grouped q-value and BH-adjusted p-value
add_grouped_qvalue_and_bh <- function(df, pval_col = "pval", comparison_col = "comparison") {
  
  pval_sym <- sym(pval_col)
  comparison_sym <- sym(comparison_col)
  
  df_corrected <- df %>%
    group_by(!!comparison_sym) %>%
    
    # Calculate q_value
    mutate(
      q_value = {
        p_vec <- !!pval_sym
        valid_indices <- !is.na(p_vec) & p_vec > 0
        p_filtered <- p_vec[valid_indices]
        q_result <- rep(NA_real_, length(p_vec))
        
        if (length(p_filtered) >= 2) {
          q_values_obj <- qvalue(p_filtered)
          q_result[valid_indices] <- q_values_obj$qvalues
        }
        q_result
      },
      
      # Calculate BH-adjusted p-value (recalculating the FDR for consistency/verification)
      adj_pval = {
        p_vec <- !!pval_sym
        valid_indices <- !is.na(p_vec) & p_vec > 0
        p_filtered <- p_vec[valid_indices]
        bh_result <- rep(NA_real_, length(p_vec))
        
        if (length(p_filtered) >= 1) {
          bh_calculated <- p.adjust(p_filtered, method = "BH")
          bh_result[valid_indices] <- bh_calculated
        }
        bh_result
      }
    ) %>%
    ungroup()
  
  return(df_corrected)
}

limma_result <- add_grouped_qvalue_and_bh(limma_result, pval_col = "pval", comparison_col = "comparison")

limma_result <- limma_result %>%
  filter(!is.na(diff), !is.na(adj_pval_limma)) %>%
  left_join(genes, by = c("ptm_collapse_key" = "ptm_collapse_key"))
limma_result %>% count(Contrast_Type)
View(limma_result)

## Save data ---------------------------
write_csv(limma_result, diff_output)