#Ria Raut, Jan 20th, 2026 
#plotting script, remade 

#housekeeping ---- 
library(tidyverse) 
library(ggplot2)

#Data cleaning and processing - mostly stolen from Michael  ----
#don't question it when you're adding new data, just run the code and call it a day 

#First let's read in all the results files
path <- "C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW/outputs/statistical"
files <- list.files(
  path = path,
  pattern = "^results.*\\.csv$",
  full.names = TRUE)
results <- lapply(files, read_csv)
#some of them have num_losses_mean -- let's just delete that column from everyone 
results <- lapply(
  results,
  \(df) dplyr::select(df, -dplyr::any_of("num_losses_mean"))
)
#let's combine all into one big dataset, with a 'name' column showing the dataset 
combined_results <- bind_rows(results, .id = "name")

#Some stolen Michael coding trying to organize the names all nice 
fn_outputs = dir('outputs/statistical',pattern="result.*\\.csv",full.names = TRUE)

truncate_name <- function(fn) {gsub('\\.csv','',gsub('outputs/statistical/results_','',fn))}
names_nice = truncate_name(fn_outputs)

#data cleaning up -- organie the names and treatments, add some additional information, do the 99.5 quantile trip 
names_nice_orig = names_nice
names_nice = str_to_sentence(gsub("_"," ",names_nice))
names_nice = gsub("Grassland ","Grassland\n",names_nice)
names(names_nice) = names_nice_orig

varnames_nice = c(experimental_design='Experimental design',
                  n_sp='Dataset, total # species',
                  type='Dataset, type',
                  #num_losses_mean='Outcome, mean # species lost',
                  n_env_levels='Dataset, # of environments',
                  deterministic='Dataset, deterministic dynamics',
                  nice_name='Dataset')
experimental_design_nice = c(mixed='Mixed') 
methods_nice = c(
  rf='Random forest',
  naive=' Naïve (mean abundance)')

# add NA removal counts
source('utils/quantile_trim.R')
source('src/dataset_stats.R')
fns = sprintf('data/%s/data_%s.csv',names(names_nice),names(names_nice))
fns = gsub('data/human_gut','data/human_and_mouse_gut',fns) #YOU NEED THESE TWO LINES BECAUSE MOUSE AND HUMAN DATA ARE IN THE SAME FILE 
fns = gsub('data/mouse_gut','data/human_and_mouse_gut',fns)


dataset_stats_all = rbindlist(lapply(fns, get_dataset_stats))

# add some additional info
df_all = combined_results %>% 
  mutate(name_nice = names_nice[name]) %>%
  mutate(experimental_design_nice = factor(experimental_design_nice[experimental_design],levels=experimental_design_nice,ordered=TRUE)) %>%
  mutate(method_nice = factor(methods_nice[method],levels=methods_nice,ordered=TRUE)) %>%
  left_join(dataset_stats_all,by='name') %>%
  # mutate(deterministic = name %in% c('grassland_annual_plants','grassland_annual_plants_drought','human_gut','mouse_gut')) %>%
  # mutate(empirical = name %in% c("ciliates","fly_gut","fruit_flies","prairie_plants","soil_bacteria")) %>% # if data was empirical or modelled 
  mutate(abundance_mae_mean_test_scaled = 
           abundance_mae_mean_test / q_995) %>%
  mutate(abundance_mae_mean_test_scaled_clipped = 
           ifelse(abundance_mae_mean_test_scaled > 10, NA, abundance_mae_mean_test_scaled)) %>%
  mutate(name_nice = factor(name_nice)) %>%
  mutate(name = factor(name))

#now prune down to only the vital information required for 
df_cleaned <- df_all %>% 
  select(name_nice, method, experimental_design_nice, replicate_index, num_train, num_test, 
         abundance_mae_mean_test, n_cases, n_na, n_sp, n_env, n_env_levels, n_combos, 
         q_995, abundance_mae_mean_test_scaled, abundance_mae_mean_test_scaled_clipped)


#now filter by the training data size -- can change it up 
df_21_only <- df_cleaned %>%
  filter(num_train == 21)

write.csv(df_21_only, 'outputs/statistical/21_training_final_dataset.csv')

#Plotting graphs, unstolen ---- 

df_21 <- read.csv('outputs/statistical/21_training_final_dataset.csv')

graph_21 <- ggplot(df_21, aes(x =name_nice , y = abundance_mae_mean_test_scaled, color = method)) + geom_point() + geom_jitter(width = 0.05)
graph_21

#when you get better graph ideas you can export them in the output/figures folder 