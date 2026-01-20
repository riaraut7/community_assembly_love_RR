#Ria Raut, Jan 20th, 2026 
#plotting script, remade 

#housekeeping ---- 
library(tidyverse) 
library(ggplot2)

#Data cleaning and processing - mostly stolen from Michael  ----

fn_outputs = dir('outputs/statistical',pattern="result.*\\.csv",full.names = TRUE)

df_all_raw = rbindlist(lapply(1:length(fn_outputs), function(i) {
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  return(df_this)
}))

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
                  num_losses_mean='Outcome, mean # species lost',
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
df_all = df_all_raw %>% 
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

df_43_only <- df_cleaned %>%
  filter(num_train == 43)