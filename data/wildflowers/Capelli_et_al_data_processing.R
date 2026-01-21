#Cappelli et al. data processing 
#Ria Raut, August 23rd 

library(tidyverse)
library(janitor)


setwd("C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/LOVE_git_scripts/new_data/Cappelli_et_al_2022_data") 

      

#let's read in all the semicolon-separted value files first ---- 
abundance_semicolon_df <- data.frame(str_split(scan(file = "Abundance.csv", character()), ";", simplify = T))
original_species_abundance <- abundance_semicolon_df %>% 
  janitor::row_to_names(row_number = 1)

biomass_semicolon_df <- data.frame(str_split(scan(file = "Biomass.csv", character()), ";", simplify = T))
species_biomass <- biomass_semicolon_df %>% 
  janitor::row_to_names(row_number = 1)

herbivory_semicolon_df <- data.frame(str_split(scan(file = "Insect_herbivory.csv", character()), ";", simplify = T), stringsAsFactors = F)
species_herbivory <- herbivory_semicolon_df %>% 
  janitor::row_to_names(row_number = 1)

pathogen_semicolon_df <- data.frame(str_split(scan(file = "Pathogen_infection.csv", character()), ";", simplify = T))
species_pathogen <- pathogen_semicolon_df %>% 
  janitor::row_to_names(row_number = 1)

plot_info <- data.frame(str_split(scan(file = 'Plot_information.csv', character()), ';', simplify = T), stringsAsFactors = F)
plot_info <- plot_info %>% 
  janitor::row_to_names(row_number = 1)

#remove some junk 
rm(abundance_semicolon_df, biomass_semicolon_df, herbivory_semicolon_df, pathogen_semicolon_df)

#Original species abundance list cleanup ----- 
wider_abundance <- original_species_abundance %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>%
  rename_with(~ paste0(.x, ".outcome")) %>% 
  rename(Plot = Plot.outcome) %>% 
  rename(Harvest = Harvest.outcome)

#this syntax is suggested by chatGPT, sorry Ben I tried by myself first and couldn't figure it out 
for (col in grep("\\.outcome$", names(wider_abundance), value = TRUE)) { #for each column  name that contains .outcome
  species <- sub("\\.outcome$", "", col) #get the species name, which is antyhign other than the .outcome 
  wider_abundance[[paste0(species, ".action")]] <- ifelse(is.na(wider_abundance[[col]]), 0, 1) #and make another column with species + .action, where the value depends on df's outcome col? 
}

wider_abundance <- wider_abundance %>% 
  dplyr::select(c(Plot, Harvest, 
         contains('action'), 
         contains('outcome')))

#let's just take this for now as abndance to biomass, start to end state dataset
#write.csv(wider_abundance, 'data_function1_cappelli_grasses.csv', row.names = F)


#herbivory dataset cleanup ---- 
#first let's deal with the extra random rows 
clean_species_herbivory <- species_herbivory %>% 
  filter(!is.na(Species)) %>% #this removes all the random columns when notes were read in 
  filter(Species != '') %>% 
  mutate(Herbivory = na_if(Herbivory, 'NA')) %>%
  mutate(Notes  = ifelse(Notes == 'set','set to associated monoculture value', NA))


#calculate herbivory per plot per year 
clean_species_herbivory$Herbivory <- as.numeric(clean_species_herbivory$Herbivory)

estimated_herbivory <- clean_species_herbivory %>% 
  group_by(Plot, Harvest) %>% 
  summarise(Calculated_herbivory = sum(Herbivory, na.rm = T), 
            n_herb_occurences = sum(!is.na(Herbivory)), 
            .groups = 'drop')  


#pathogen dataset cleanup ----- 
clean_species_pathogen <- species_pathogen %>% 
  filter(!is.na(Species)) %>% 
  filter(Species != '') %>% 
  mutate(Infection = na_if(Infection, 'NA')) %>%
  mutate(Notes  = ifelse(Notes == 'set','set to associated monoculture value', NA))

clean_species_pathogen$Infection <- as.numeric(clean_species_pathogen$Infection)

estimated_infection <- clean_species_pathogen %>% 
  group_by(Plot, Harvest) %>% 
  summarise(Measured_infection = sum(Infection, na.rm = T), 
            n_infect_occurences = sum(!is.na(Infection)), 
            .groups = 'drop')

#plot_info cleanup and everything joinup  ---- 
clean_plot_info <- plot_info %>% 
  filter(!is.na(Nitrogen)) %>% 
  filter(Nitrogen != '')

total_cappelli_processed <- left_join(wider_abundance, estimated_herbivory, by = c('Plot', 'Harvest')) %>% 
  left_join(., estimated_infection, by = c('Plot', 'Harvest')) %>% 
  left_join(., clean_plot_info, by = 'Plot')

replaced_NAs <- total_cappelli_processed %>%
  mutate(across(
    ends_with(".outcome"),
    ~ replace_na(as.numeric(.x), 0)
  ))

write.csv(replaced_NAs, 'LOVE_updated_data_cappelli_grasses.csv', row.names = F)

#just remember that when you're using this data to train LOVE, you gotta split it up into 4 groups (by enviromental condition)
#OR include Nitrogen and fungal treatments as training variables! Each group is !84 plots iirc 