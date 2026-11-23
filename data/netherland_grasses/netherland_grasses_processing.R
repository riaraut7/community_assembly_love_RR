#Processing netherland_grasses data 
#Ria Raut, Nov 20th 

#Basic housekeeping and loading data ---- 
library(tidyverse)
getwd()
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW/data/netherland_grasses')

#this is the .action/starting dataset: 
starting_dataset <- read.csv('De_long_et_al_planting_data.csv', stringsAsFactors = T)
treatment_1 <- read.csv('Heinen_et_al_treatment_1.csv')
desirability <- read.csv('Heinen_et_al_desirability_data.csv')

#some cleaning for each dataset ---- 

#cleaning column names for the starting dataset
starting_dataset <- starting_dataset %>%
  rename(community = From.DeLong.et.al.) %>% 
  select(-count.of.species)

action_values <- starting_dataset %>% 
  rename_with(~paste0(gsub(" ", "_", .x)), 3:26)
  

#cleaning column names for the outcome values 
outcome_values <- treatment_1 %>%
  rename_with(~ gsub(" ", "_", .x), 8:81) %>% 
  select(-c(block, speed, forb.1, grass, total_cover))
  
#In outcomes, seeing how many species are introduced, then summing them all as 'introduced' 
missing_cols <- setdiff(names(outcome_values), names(starting_dataset))
missing_cols #all introduced species -- columns missing from the action -- to be summed 
missing_cols_from_outcome <- setdiff(names(action_values), names(outcome_values))
missing_cols_from_outcome #species that were introduced but aren't found there anymore -- to be added 

cols_to_ignore <- c('plot', 'subplot', 'legacy')
cols_to_sum <- setdiff(missing_cols, cols_to_ignore)

outcome_values1 <- outcome_values %>% 
  mutate(introduced_species.desirability = rowSums(select(., all_of(cols_to_sum)), na.rm = T)) %>% 
  select(-all_of(cols_to_sum)) 

#let's add the species we missed that should have columns 
outcome_values2 <- outcome_values1 %>%
  mutate(!!!setNames(rep(list(0), length(missing_cols_from_outcome)), missing_cols_from_outcome))

#now let's organize it by alphabetical order, for both 
exclude <- c('plot', 'subplot', 'forb', 'community', 'legacy', 'introduced_species.desirability')
outcome_values3 <- outcome_values2 %>%
  select(all_of(exclude), sort(setdiff(names(.), exclude)))
#  select(1:5, sort(names(.)[-(1:5)]))

action_values1 <- action_values %>%
  select(1:2, sort(names(.)[-(1:2)]))

#now lets see what needs to be added or removed 
colnames(action_values1) 
colnames(outcome_values3)

outcome_values3 <- outcome_values3 %>% 
  rename_with(~ paste0(.x, ".outcome"), 7:30) %>% 
  mutate(legacy.initial = legacy) %>% 
  select(-legacy)

#in action column, adding the .action 
action_values2 <- action_values1 %>% 
  rename_with(~ paste0(.x, '.action'), 3:26)

colnames(action_values2) 
colnames(outcome_values3)

#now let's join both based on the values found in the forb and community columns 
combined <- action_values2 %>%
  right_join(outcome_values3, by = c("community", "forb"))
#to see how many unique columns w actions and outcomes -- should be the same! 
names(combined)[str_detect(names(combined), "\\.action")]
names(combined)[str_detect(names(combined), "\\.outcome")]


#Let's add a desirability section 
desirability2 <- desirability %>% 
  mutate(AMF.desirability = AMF, 
         Saprotrophs.desirability = Saprotrophs, 
         total_pathogen.desirability = TotalPathogens) %>% 
  select (c(plot, subplot, community, forb, legacy, 
             AMF.desirability, Saprotrophs.desirability, total_pathogen.desirability))
  
total_combined <- combined %>% 
  right_join(desirability2, by = c('community', 'forb', 'plot', 'subplot'))

final_netherland_grasses <- total_combined %>% 
  select (c (plot, subplot, community, forb, 
             contains('.action'), 
             contains('.initial'), 
             contains('.outcome'), 
             contains('.desirability')
             ))

colnames(final_netherland_grasses)

write.csv(final_netherland_grasses, 'data_netherland_grasses.csv', row.names = F)
