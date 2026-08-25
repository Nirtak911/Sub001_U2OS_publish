# 1. Load libraries ---------------------------
library(tidyverse)
library(protti)
library(reshape2)
library(stringr)

# Ensure data_in is defined
# data_in <- "path/to/data"

## Load data ---------------------------
input <- file.path(data_in, "04_scale_not_abs_loess_data_log2.csv")
meta  <- file.path(data_in, "01_meta_data.csv")
genes <- file.path(data_in, "01_genes.csv")

data  <- read_protti(input)
meta  <- read_csv(meta)
genes <- read_csv(genes)

perm_output <- file.path(data_in, "07_0_scaled_permutation_fast_annotated.csv")

# Clean up r_condition
data <- data %>% 
  dplyr::mutate(r_condition = stringr::str_remove(r_condition, "^U2OS_")) %>%
  dplyr::mutate(
    parts = stringr::str_split(r_condition, pattern = "_", n = 2, simplify = TRUE),
    r_condition = paste(parts[, 2], parts[, 1], sep = "_") # Swap to Frac_Treat
  ) %>%
  dplyr::select(-parts)

perm_data <- data %>%
  dplyr::select(
    Protein.Group = pg_protein_groups,
    value = scaled_intensity,
    replicate = r_replicate,
    condition = r_condition
  ) %>%
  tidyr::separate(
    condition,
    into = c("fraction", "treatment"),
    sep = "_"
  )

# Assign unique integer replicate IDs per condition
perm_data <- perm_data %>%
  mutate(
    replicate_perm = case_when(
      treatment == "CTRL" ~ as.integer(replicate),
      TRUE ~ as.integer(replicate) + max(as.integer(replicate))
    )
  )

perm_selected <- perm_data %>%
  filter(treatment %in% c("CTRL", "Anis"))

# Get observed replicate sets
control_obs <- perm_selected %>%
  filter(treatment == "CTRL") %>%
  pull(replicate_perm) %>%
  unique()

treat_obs <- perm_selected %>%
  filter(treatment == "Anis") %>%
  pull(replicate_perm) %>%
  unique()


# 2. Define Helper Functions ---------------------------

compute_stat <- function(df, treat_reps, control_reps) {
  df$value <- 2^(df$value)
  df$value[is.na(df$value)] <- 0.01
  
  # Use replicate_perm to match treat_reps / control_reps
  value_mat <- xtabs(value ~ fraction + replicate_perm, data = df) 
  
  # Cast to character to avoid matrix index positional errors
  treat_cols   <- as.character(treat_reps)
  control_cols <- as.character(control_reps)
  
  mean_treat <- rowMeans(value_mat[, treat_cols, drop = FALSE])
  mean_ctrl  <- rowMeans(value_mat[, control_cols, drop = FALSE])
  
  scaled_treat   <- mean_treat / sum(mean_treat)
  scaled_control <- mean_ctrl / sum(mean_ctrl)
  
  diff <- scaled_treat - scaled_control
  fractions <- rownames(value_mat)
  
  up_idx   <- which.max(diff)
  down_idx <- which.min(diff)
  
  T_stat <- diff[up_idx] - diff[down_idx]
  return(list(T_stat, fractions[up_idx], fractions[down_idx]))
}

compute_stat_fc <- function(df, up, down, treat_reps) {
  default_t <- 0
  df <- df %>%
    mutate(condition = ifelse(replicate_perm %in% treat_reps, "treatment", "control"))
  
  temp_table_up <- df %>% filter(fraction == up)
  if (any(tapply(!is.na(temp_table_up$value), temp_table_up$condition, sum) <= 1)) {
    t_res_result_up <- default_t
  } else {
    t_res_up <- lm(value ~ condition, data = temp_table_up) 
    t_res_result_up <- summary(t_res_up)$coefficients["conditiontreatment", "t value"]
  }
  
  temp_table_down <- df %>% filter(fraction == down)
  if (any(tapply(!is.na(temp_table_down$value), temp_table_down$condition, sum) <= 1)) {
    t_res_result_down <- default_t
  } else {
    t_res_down <- lm(value ~ condition, data = temp_table_down) 
    t_res_result_down <- summary(t_res_down)$coefficients["conditiontreatment", "t value"]
  }
  
  T_stat <- sqrt(t_res_result_up^2 + t_res_result_down^2)
  return(list(T_stat, temp_table_up, temp_table_down))
}

prot_sig <- function(df_protein, prot) {
  all_reps <- unique(df_protein$replicate_perm)
  n_treat <- length(treat_obs)
  treat_sets <- combn(all_reps, n_treat, simplify = FALSE)
  
  T_obs    <- compute_stat(df_protein, treat_obs, control_obs)[[1]]
  Up_obs   <- compute_stat(df_protein, treat_obs, control_obs)[[2]]
  Down_obs <- compute_stat(df_protein, treat_obs, control_obs)[[3]]
  
  T_obs_fc       <- compute_stat_fc(df_protein, Up_obs, Down_obs, treat_obs)[[1]]
  fc_matrix_up   <- compute_stat_fc(df_protein, Up_obs, Down_obs, treat_obs)[[2]]
  fc_matrix_down <- compute_stat_fc(df_protein, Up_obs, Down_obs, treat_obs)[[3]]
  
  T_perm    <- numeric(length(treat_sets))
  Up_perm   <- numeric(length(treat_sets))
  Down_perm <- numeric(length(treat_sets))
  T_perm_fc <- numeric(length(treat_sets))
  
  for (j in seq_along(treat_sets)) {
    treat_i   <- treat_sets[[j]]
    control_i <- setdiff(all_reps, treat_i)
    
    T_perm[j]    <- compute_stat(df_protein, treat_i, control_i)[[1]]
    Up_perm[j]   <- compute_stat(df_protein, treat_i, control_i)[[2]]
    Down_perm[j] <- compute_stat(df_protein, treat_i, control_i)[[3]]
    
    T_perm_fc[j] <- compute_stat_fc(df_protein, Up_perm[j], Down_perm[j], treat_i)[[1]]
  }
  
  p_value    <- mean(T_perm >= T_obs) 
  p_value_fc <- mean(T_perm_fc >= T_obs_fc)
  
  return(tibble(
    pvalue    = p_value,
    pvalue_fc = p_value_fc,
    Up_obs    = Up_obs,
    Down_obs  = Down_obs
  ))
}

# 3. Execution Block ---------------------------

results <- perm_selected %>%
  group_by(Protein.Group) %>%
  group_modify(~ prot_sig(.x, .y$Protein.Group))

# Annotate with Gene Names
results <- results %>%
  left_join(genes, by = c("Protein.Group" = "pg_protein_groups"))
View(results)

# Export results
write_csv(results, perm_output)