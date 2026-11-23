
#load the libraries and etc 
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')
library(tidyverse)

#graphing for function f, naive vs RF ---- 
folder1 <- 'outputs/statistical'
combined_data <- list.files(
  path = folder1,
  pattern = "\\.csv$",
  full.names = TRUE
) %>%
  set_names(nm = basename(.)) %>%
  map_dfr(
    ~ read_csv(.x, show_col_types = FALSE),
    .id = "dataset"
  )
#graphing for desirability ---- 
#load all the datasets and compile into one dataset 
# Set folder containing CSV files
folder2 <- "outputs_desirability"

# Read all CSV files and combine them
combined_data <- list.files(
  path = folder2,
  pattern = "\\.csv$",
  full.names = TRUE
) %>%
  set_names(nm = basename(.)) %>%
  map_dfr(
    ~ read_csv(.x, show_col_types = FALSE),
    .id = "dataset"
  )

cleaned_names <- combined_data %>%
  mutate(dataset = str_remove(dataset, "results_")) %>%
  mutate(dataset = str_remove(dataset, "\\.csv"))

#okay let's keep it for training size 10 
df_20 <- cleaned_names %>% 
  filter(training_size == 10)
