#load libraries
library(tidyverse)
library(protti)
library(plotly)
library(viridis)
library(ggprism)
library(gridExtra)
library(patchwork)

## Load data---------------------
input <- file.path(data_in,"/04_scale_not_abs_loess_data_log2.csv")
data <- read_csv(input)
input2 <- "raw/Subcell_markers_selected.csv"
selected_markers <- read_protti(input2)

## output data------------------
corr_output <- file.path(data_in,"/04_4_scaled_intensity_profile_inter_vs_intra_corr.csv")
corr_output_combined <- file.path(data_in,"/04_4_scaled_intensity_profile_inter_vs_intra_corr_combined.csv")

View(data)
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

## Add localisation data---------------------
# only keep the first pg_protein_groups
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, "^;", "")
data$pg_protein_groups <- str_replace_all(data$pg_protein_groups, ";.*", "")

# 1. Reshape data incorporating Treatment, Fraction, and r_replicate
data_wide <- data %>%
  # Creates an ID like "CTRL_Frac1_Rep1" or "TREATED_Frac1_Rep1"
  mutate(sample = paste0(treatment, "_", fraction, "_", r_replicate)) %>%
  select(pg_protein_groups, sample, scaled_intensity) %>%
  pivot_wider(names_from = sample, values_from = scaled_intensity) %>%
  drop_na() 

# 2. Compute the global correlation matrix
cor_matrix <- data_wide %>%
  select(-pg_protein_groups) %>%
  cor(method = "pearson")

# 3. Convert matrix to a long data frame and parse the IDs
cor_df <- as.data.frame(cor_matrix) %>%
  rownames_to_column(var = "Sample_A") %>%
  pivot_longer(cols = -Sample_A, names_to = "Sample_B", values_to = "Correlation") %>%
  # Split the column names back into distinct metadata variables
  separate(Sample_A, into = c("Treatment_A", "Frac_A", "Rep_A"), sep = "_") %>%
  separate(Sample_B, into = c("Treatment_B", "Frac_B", "Rep_B"), sep = "_") %>%
  # CRUCIAL: Only compare samples within the same treatment group
  filter(Treatment_A == Treatment_B) %>%
  rename(Treatment = Treatment_A) %>%
  # Remove self-correlations
  filter(!(Frac_A == Frac_B & Rep_A == Rep_B))

# 4. Extract r_replicate Correlations (Same Treatment, Same Fraction, Diff r_replicate)
corr_r_replicates <- cor_df %>%
  filter(Frac_A == Frac_B & Rep_A != Rep_B) %>%
  rename(Corr_Rep = Correlation) %>%
  select(Treatment, Fraction = Frac_A, Corr_Rep) %>%
  distinct() 

# 5. Extract Fraction Correlations (Same Treatment, Diff Fraction, Same r_replicate)
corr_fractions <- cor_df %>%
  filter(Frac_A != Frac_B & Rep_A == Rep_B) %>%
  group_by(Treatment, Frac_A, Frac_B) %>%
  summarise(Corr_Frac = mean(Correlation), .groups = 'drop') 

# 6. Combine and calculate the final ratios per Treatment
final_ratios <- corr_fractions %>%
  left_join(corr_r_replicates, by = c("Treatment", "Frac_A" = "Fraction")) %>%
  mutate(Ratio = Corr_Rep / abs(Corr_Frac))

View(final_ratios)

cor_df_combined <- as.data.frame(cor_matrix) %>%
  rownames_to_column(var = "Sample_A") %>%
  pivot_longer(cols = -Sample_A, names_to = "Sample_B", values_to = "Correlation") %>%
  separate(Sample_A, into = c("Treatment_A", "Frac_A", "Rep_A"), sep = "_") %>%
  separate(Sample_B, into = c("Treatment_B", "Frac_B", "Rep_B"), sep = "_") %>%
  # Remove identical sample self-correlations
  filter(!(Frac_A == Frac_B & Rep_A == Rep_B & Treatment_A == Treatment_B))

corr_r_replicates_combined <- cor_df_combined %>%
  # Same fraction, but eliminate the exact same sample identity
  filter(Frac_A == Frac_B) %>% 
  group_by(Frac_A) %>% 
  summarise(
    Corr_Rep = mean(Correlation), 
    .groups = 'drop'
  ) %>%
  rename(Fraction = Frac_A)

corr_fractions_combined <- cor_df_combined %>%
  filter(Frac_A != Frac_B) %>%  
  group_by(Frac_A, Frac_B) %>%  
  summarise(
    Corr_Frac = mean(Correlation), 
    .groups = 'drop'
  )

final_ratios_combined <- corr_fractions_combined %>%
  left_join(corr_r_replicates_combined, by = c("Frac_A" = "Fraction")) %>%
  mutate(Ratio = Corr_Rep / abs(Corr_Frac))

View(final_ratios_combined)

# save files
write_csv(final_ratios, corr_output)
write_csv(final_ratios_combined, corr_output_combined)