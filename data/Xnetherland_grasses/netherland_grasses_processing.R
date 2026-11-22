#Processing netherland_grasses data 
#Ria Raut, Nov 20th 

#Basic housekeeping and loading data ---- 
library(tidyverse)
getwd()
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW/data/Xnetherland_grasses')

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
cols_to_ignore <- c('plot', 'subplot', 'legacy')
cols_to_sum <- setdiff(missing_cols, cols_to_ignore)

outcome_values1 <- outcome_values %>% 
  mutate(introduced_species.desirability = rowSums(select(., all_of(cols_to_sum)), na.rm = T)) %>% 
  select(-all_of(cols_to_sum))

outcome_values2 <- outcome_values1 %>% 
  rename_with(~ paste0(.x, ".outcome"), 6:23) %>% 
  mutate(legacy.initial = legacy) %>% 
  select(-legacy)


#in action column, adding the .action 
action_values <- action_values %>% 
  rename_with(~ paste0(.x, '.action'), 3:26)

#now let's join both based on the values found in the forb and community columns 
combined <- action_values %>%
  right_join(outcome_values2, by = c("community", "forb"))

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


write.csv(final_netherland_grasses, 'netherland_grasses_data.csv', row.names = F)
