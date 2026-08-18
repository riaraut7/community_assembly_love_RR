#Ria Raut, Jan 20th, 2026 
#plotting script, remade 

#housekeeping ---- 
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(wesanderson)
library(tibble)
library(ggpubr)
library(MuMIn)
library(caret)
library(visreg)
library(lme4)
library(vegan)
library(RColorBrewer)
library(MuMIn)
library(conflicted)
library(ggbiplot)
library(ggrepel)
library(car)
library(sjPlot)
library(ggheatmap)
library(DHARMa)
library(e1071)
library(stringr)
library(ggsci)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "plyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarize", "dplyr")
conflicts_prefer(reshape::melt)


palette_11 = pal_d3(palette="category20")(11)

getwd()

#Michael's nicer script ---- 

###
if (!file.exists('outputs/figures')){
  dir.create('outputs/figures')  
}

###

fn_outputs = dir('outputs/statistical',pattern="result.*\\.csv",full.names = TRUE)

df_all_raw = rbindlist(lapply(1:length(fn_outputs), function(i){
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  return(df_this)
}), 
fill = TRUE

)

truncate_name <- function(fn) {gsub('\\.csv','',gsub('outputs/statistical/results_','',fn))}

names_nice = truncate_name(fn_outputs)
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
                  nice_name='Dataset'
)

experimental_design_nice = c(mixed='Mixed',
                             `low-2`='Doublets',
                             `low-3`='Doublets, triplets',
                             `high-1`='1-dropouts',
                             `high-2`='1-dropouts, 2-dropouts',
                             prior='Doublets, 1-dropouts; then mixed')

methods_nice = c(
  rf='Random forest',
  # glv='GLV predictions',
  # glv_rf='Random forest on GLV residuals',
  # glv_rf_full='Random forest + GLV predictions',
  # sequential_rf='Random forest, sequential',
  naive=' Naïve (mean abundance)')

# add NA removal counts
source('utils/quantile_trim.R')
fns = sprintf('data/%s/data_%s.csv',names(names_nice),names(names_nice))
fns = gsub('data/human_gut','data/human_and_mouse_gut',fns)
fns = gsub('data/mouse_gut','data/human_and_mouse_gut',fns)

source('src/dataset_stats.R')

dataset_stats_all = rbindlist(lapply(fns, get_dataset_stats))

# add some additional info
df_all = df_all_raw %>% 
  mutate(name_nice = names_nice[name]) %>%
  mutate(experimental_design_nice = factor(experimental_design_nice[experimental_design],levels=experimental_design_nice,ordered=TRUE)) %>%
  mutate(method_nice = factor(methods_nice[method],levels=methods_nice,ordered=TRUE)) %>%
  left_join(dataset_stats_all,by='name') %>%
  mutate(deterministic = name %in% c('grassland_annual_plants','grassland_annual_plants_drought','human_gut','mouse_gut')) %>%
  mutate(empirical = name %in% c("ciliates","fly_gut","fruit_flies","prairie_plants","soil_bacteria")) %>%
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
  dplyr::filter(num_train == 21)
write.csv(df_21_only, 'outputs/statistical/21_training_final_dataset.csv')

  #also going to try for ten number of training states (?)) 
df_10_only <- df_cleaned %>%
  dplyr::filter(num_train == 10)
write.csv(df_10_only, 'outputs/statistical/10_training_final_dataset.csv')

#Plotting graphs, unstolen ---- 

df_21 <- read.csv('outputs/statistical/21_training_final_dataset.csv')
graph_21 <- ggplot(df_21, aes(x =name_nice , y = abundance_mae_mean_test_scaled, color = method)) + geom_point() + geom_jitter(width = 0.05)
graph_21

df_10 <- read.csv('outputs/statistical/10_training_final_dataset.csv')
graph_10 <- ggplot(df_10, aes(x =name_nice , y = abundance_mae_mean_test_scaled, color = method)) + geom_point() + geom_jitter(width = 0.05)
graph_10

colnames(df_10)
data_and_env_number <- df_10 %>% 
  dplyr::select(c(name_nice, n_env_levels, n_env))
