
#load the libraries and etc 
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')
library(tidyverse)

#Function f: Data loading and cleaning ---- 
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

cleaned_names <- combined_data %>% 
  mutate(dataset = str_remove(dataset, "results_")) %>%
  mutate(dataset = str_remove(dataset, "\\.csv")) %>% 
  select(c(dataset, method, num_train, abundance_mae_mean_test)) %>% 
  rename(mae_test = abundance_mae_mean_test)

#Let's write this into a csv so it's easier to use for next time 
unique(cleaned_names$dataset)  
write.csv(cleaned_names, 'combined_function_f_results.csv', row.names = F)

#let's get only those w training set of 10 
#need to run mouse_gut, human_gut, and forest_trees again 

#Function f: graph 1: MAE for naive vs RF ---- 
full_dataset <- read.csv('combined_function_f_results.csv', stringsAsFactors = T)
df_train_21 <- full_dataset %>% 
  filter(num_train == 21) %>% 
  filter(dataset != 'fly_gut') %>% 
  filter(dataset != 'grassland_annual_plants_drought')
unique(df_train_21$dataset) #12 datasets present -- will add at least 3 more 

rf_vs_naive <- ggplot(df_train_21, aes(x = dataset, y = mae_test, color = method)) + 
  geom_jitter(
    size = 4.5,     
    alpha = 0.5,   
    width = 0.1,    
    height = 0) + 
  scale_color_manual (values = c( 'naive' = 'burlywood', 'rf' = 'darkolivegreen')) + 
  labs(x = 'Dataset names', 
       y = 'Standardized mean absolute error') + 
  theme_classic() + 
  theme(
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
        ) #+ 
#  coord_cartesian(ylim = c(0, 20))
rf_vs_naive

  #statistical test ---- 
#for each dataset, pair the naive w rf 

stat_test <- df_train_21 %>%
  group_by(dataset, num_train, method) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = method,
    values_from = mae_test
  ) %>%
  select(-replicate)


results <- stat_test %>%
  group_by(dataset) %>%
  summarise(
    t_test = list(t.test(rf, naive, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, \(x) x$p.value),
    t_statistic = sapply(t_test, \(x) x$statistic)
  )

total_p <- mean(results$p_value)


#Function f: graph 2 MAE by training size ---- 
only_rf <- full_dataset %>% 
  filter(method == 'rf') %>% 
  filter(dataset != 'fly_gut')

function_f_mae <- ggplot(only_rf, aes(x = num_train, y = mae_test, color = dataset)) + 
  geom_jitter(
    size = 4,     
    alpha = 0.6,   
    width = 0.1,    
    height = 0)   + 
  coord_cartesian(xlim = c(0, 150), 
                  ylim = c(0, 30)) + 
  labs(x = 'Training size', 
       y = 'Standardized mean absolute error predicting community outcomes') + 
  scale_color_viridis_d(option = 'D', 
                        direction = 1) + 
  #  scale_color_binned(palette = 'GnBu') +  
  theme_classic() + 
  theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14)
  ) 
function_f_mae

#Function f: graph 3: env types and MAE on x and y ---- 
#let's first edit the full dataset 
full_dataset <- read.csv('combined_function_f_results.csv', stringsAsFactors = T)
unique(full_dataset$dataset)

edited_database <- full_dataset %>%
  mutate(`num_environments` = case_when(
    dataset == "california_grasses" ~ 3,
    dataset == "ciliates" ~ 6,
    dataset == "fly_gut" ~ 1,
    dataset == "forest_trees" ~ 1,
    dataset == "fruit_flies" ~ 4,
    dataset == "grassland_annual_plants" ~ 1,
    dataset == "grassland_annual_plants_drought" ~ 1,
    dataset == "human_gut" ~ 1,
    dataset == "mouse_gut" ~ 1,
    dataset == "jena_wildflowers" ~ 2,
    dataset == "prairie_plants" ~ 1,
    dataset == "soil_bacteria" ~ 1,
    dataset == "wildflowers" ~ 4,
    dataset == "grassland_diversity" ~ 4,
    dataset == "parasite_host_diversity" ~ 1,
    dataset == "tree_colonization" ~ 4,
    dataset == "netherland_grasses" ~ 2,
    TRUE ~ NA_real_
  )) %>% 
  filter(method != 'naive')

edited_database$num_environments <- as.factor(edited_database$num_environments)

mae_by_env_types <- ggplot(edited_database, aes(x = num_train, y = mae_test, color = num_environments)) + 
  geom_jitter(
    size = 4,     
    alpha = 0.6,   
    width = 0.1,    
    height = 0)   + 
  coord_cartesian(xlim = c(0, 150), 
                  ylim = c(0, 30)) + 
  labs(x = 'Number of experiments (communities) used for training', 
       y = 'Standardized mean absolute error') + 
# scale_color_binned(palette = 'Dark2') +  
  scale_color_manual (values = c( '6' = '#512a18', '4' = '#833717',
                                 # '3' = '#e07e12', 
                                  '2' = '#e07e12', '1' = '#c4ae3d')) + 
  
#  scale_color_viridis_b(option = "D") + 
  theme_classic() + 
  theme(
#    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) 
  
mae_by_env_types  

#combine all datasets, but then have a new column w num environment types. num on x, MAE on y, 
#colour code by dataset 



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
unique(cleaned_names$dataset)
#function g graph 1: graph MAE on y, training size on x, colour coe by dataset ---- 
function_G_mae <- ggplot(cleaned_names, aes(x = training_size, y = diversity.desirability, color = dataset)) + 
  geom_jitter(
    size = 4,     
    alpha = 0.6,   
    width = 0.1,    
    height = 0)   + 
  coord_cartesian(xlim = c(0, 150)) + 
  labs(x = 'Training size', 
       y = 'Standardized mean absolute error predicting diversity') + 
  scale_color_viridis_d(option = 'D', 
                        direction = 1) + 
#  scale_color_binned(palette = 'GnBu') +  
  theme_classic() + 
  theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12)
  ) 
function_G_mae

#function g graph 2: direct vs indirect ----- 
#have the columns: dataset, desir type, MAE, indirect vs direct. All for trainig size 10. Then 
#graph direct vs indirect, comparison of values 
#do a statistical test 
#okay let's keep it for training size 10 
df_20 <- cleaned_names %>% 
  filter(training_size == 20)

df_20_long <- pivot_longer(df_20, 
                           cols = contains('.desirability'), 
                           names_to = 'community_characteristics', 
                           values_to = 'MAE')
cleaned_attributes <- df_20_long %>%
  mutate(community_characteristics = str_remove(community_characteristics, "\\.desirability"))
unique(cleaned_attributes$community_characteristics)


edited_database_characteristics <- cleaned_attributes %>%
  mutate(`indirect_vs_direct` = case_when(
    community_characteristics == "alaria" ~ 'IND',
    community_characteristics == "cephalo" ~ 'IND',
    community_characteristics == "ribei" ~ 'IND',
    community_characteristics == "diversity" ~ 'DIR',
    community_characteristics == "total_seeded" ~ 'IND',
    community_characteristics == "introduced_species" ~ 'IND',
    community_characteristics == "parasite_richness" ~ 'IND',
    community_characteristics == "percent_cover" ~ 'IND',
    community_characteristics == "average_height" ~ 'IND',
    community_characteristics == "earliest_establishment" ~ 'IND',
    community_characteristics == "herb_colonization" ~ 'IND',
    community_characteristics == "herbivory" ~ 'DIR',
    community_characteristics == "infection" ~ 'IND',
    community_characteristics == "direct_sown_sla" ~ 'IND', 
    TRUE ~ NA_character_ 
    )) %>% 
  drop_na() %>% 
  filter(MAE < 10)

colnames(edited_database_characteristics)

characteristics_final <- edited_database_characteristics %>%
  group_by(dataset) %>%
  filter(all(c("IND", "DIR") %in% indirect_vs_direct)) %>%
  ungroup()

character_final_2 <- characteristics_final %>% 
  filter(dataset != 'parasite_host_diversity')


direct_vs_indirect <- ggplot(character_final_2, 
                             aes(x = indirect_vs_direct, y = MAE, fill = dataset, color = dataset)) + 
  geom_boxplot(
    alpha = 0.6,          
    outlier.size = 3,      
    outlier.alpha = 0.8
  ) +  
  labs(x = 'Type of community characteristic (direct vs indirect)', 
       y = 'Mean absolute error', 
       legend = 'Dataset') + 
  scale_color_viridis_d(option = 'D', alpha = 0.9) +
  scale_fill_viridis_d(option = 'D', alpha = 0.6)+ 
  coord_cartesian(ylim = c(0, 4)) + 
  theme_classic() + 
  theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) 
direct_vs_indirect

  #statistical test ----- 
#let's pivot wider first


wider_cc <- character_final_2 %>%
  group_by(dataset) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  pivot_wider(
    names_from = indirect_vs_direct,
    values_from = MAE
  ) %>%
  select(-replicate)

results_g <- wider_cc %>%
  group_by(dataset) %>%
  summarise(
    t_test = list(t.test(IND, DIR, paired = F)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, \(x) x$p.value),
    t_statistic = sapply(t_test, \(x) x$statistic)
  )
g_p_value <- mean(results_g$p_value)
g_p_value

#function g, graph 3: MAE by environment het ----- 

edited_database <- cleaned_names %>%
  mutate(`num_environments` = case_when(
    dataset == "california_grasses" ~ 3,
    dataset == "ciliates" ~ 6,
    dataset == "fly_gut" ~ 1,
    dataset == "forest_trees" ~ 1,
    dataset == "fruit_flies" ~ 4,
    dataset == "grassland_annual_plants" ~ 1,
    dataset == "grassland_annual_plants_drought" ~ 1,
    dataset == "human_gut" ~ 1,
    dataset == "mouse_gut" ~ 1,
    dataset == "jena_wildflowers" ~ 2,
    dataset == "prairie_plants" ~ 1,
    dataset == "soil_bacteria" ~ 1,
    dataset == "wildflowers" ~ 4,
    dataset == "grassland_diversity" ~ 4,
    dataset == "parasite_host_diversity" ~ 1,
    dataset == "tree_colonization" ~ 4,
    dataset == "netherland_grasses" ~ 2,
    TRUE ~ NA_real_
  )) %>% 
  select(c(dataset, training_size, diversity.desirability, num_environments)) %>% 
  drop_na()

edited_database$num_environments <- as.factor(edited_database$num_environments)
colnames(edited_database)
mae_by_env_types_g <- ggplot(edited_database, aes(x = training_size, y = diversity.desirability, color = num_environments)) + 
  geom_jitter(
    size = 4,     
    alpha = 0.6,   
    width = 0.6,    
    height = 0)   + 
  coord_cartesian(xlim = c(0, 150), 
                  ylim = c(0, 5)) + 
  labs(x = 'Number of experiments (communities) used for training', 
       y = 'Standardized mean absolute error for diversity') + 
  scale_color_manual (values = c( '6' = '#512a18', '4' = '#833717',
                                   '3' = '#e07e12', 
                                  '2' = '#c4ae3d', '1' = 'burlywood4')) + 
  
  #  scale_color_viridis_b(option = "D") + 
  theme_classic() + 
  theme(
    #    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) 

mae_by_env_types_g  
