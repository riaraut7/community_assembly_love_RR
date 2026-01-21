#Stubble et al. data processing 
#Ria Raut, August 23rd 

library(tidyverse)
library(janitor)


setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/LOVE_git_scripts/new_data/Stubble_et_al_2017_PRYER_info/actual_data')

Pryer_data <- read.csv('PRYER_STG.csv', stringsAsFactors = F)

#Cleaning and columns prep ----
#First take out useless columns that you won't need, ex. exp, ref, data.collection.year, growing.seasons, rep.5,
#rep.10, n.s, bare, unseeded.forbs, unseeded.grasses, Total..data.check
Pryer_cleaned <- Pryer_data %>% select(-c("exp", "ref",
                                          "data.collection.year", "site", "growing.seasons", "plant.year",    
                                          "rep.5", "rep.10", "n.s", "bare", "unseeded.forbs", "unseeded.grasses", 
                                          "Total..data.check."
                                          ))

#All the data you see is actual the .outcomes results. So let's add that to all the species 
Pryer_cleaned <- Pryer_cleaned %>% 
  rename_with(~ paste0(.x, ".outcome")) %>% 
#fix some new names -- ignore for now, can always fix later 
  rename(water = water.outcome) %>% 
  rename(treatment = treatment.outcome) %>% 
  rename(total.n = total.n.outcome) %>% 
  rename(total.i = total.i.outcome) %>% 
  rename(total.seeded = total.seeded.outcome) %>% 
  rename(problem.row = problem.row.outcome) %>% 
  rename(Comments = Comments.outcome)

species_names <- c("Brca", #Bromus carinatus  
                   "Elgl", #Elymus glaucus 
                   "Hobr", #Hordeum brachyantherum 
                   "Stpu", 
                   "Avsp", #Avena species, including barbata or fatua at Davis vs Hop and Mclaughlin respectively 
                   "Brho", #Bromus hordeaceus 
                   "Homu", #Hordeum murinum 
                   "Vusp")
#define native plants (make a list of column names) + invasives 
native_spp <- c("Brca", #Bromus carinatus  
                "Elgl", #Elymus glaucus 
                "Hobr", #Hordeum brachyantherum 
                "Stpu") #Stipa puchra 

exotic_spp <- c("Avsp", #Avena species, including barbata or fatua at Davis vs Hop and Mclaughlin respectively 
                "Brho", #Bromus hordeaceus 
                "Homu", #Hordeum murinum 
                "Vusp") #Vulpia species, including myuros and bromoides at Dais and Mclaughlin vs Hopland respectively 
#We're going to ingore Napu, which I typed somewhere but now forgot? 


#let's get .action, .mid_state, and .late_state columns 
for(spp in species_names) {
  Pryer_cleaned[[paste0(spp, ".action")]] <- 0 
  Pryer_cleaned[[paste0(spp, ".mid_state")]] <- 0 
  Pryer_cleaned[[paste0(spp, ".late_state")]] <- 0 
}

#Big step! Distribution of 1's based on treatment conditions ---- 
#let's add 1's and 0's conditionally 
#I don't know why it says all of my paramtheses and quotations are not right?? 
Pryer_allinfo <- Pryer_cleaned %>%
  rowwise() %>%
  mutate(
    # case: N
    across(ends_with(".action"),
           ~ if (treatment == "N" && str_remove(cur_column(), "\\..*") %in% native_spp) 1 else .x),
    
    # case: NE
    across(ends_with(".action"),
           ~ if (treatment == "NE" && str_remove(cur_column(), "\\..*") %in% species_names) 1 else .x), #i.e. all species 
    
    #case: NtE
    across(ends_with(".action"),
           ~ if (treatment == "NtE" && str_remove(cur_column(), "\\..*") %in% native_spp) 1 else .x),
    across(ends_with(".mid_state"),
           ~ if (treatment == "NtE" && str_remove(cur_column(), "\\..*") %in% exotic_spp) 1 else .x),
    
    #case: NttE 
    across(ends_with(".action"),
           ~ if (treatment == "NttE" && str_remove(cur_column(), "\\..*") %in% native_spp) 1 else .x),
    across(ends_with(".late_state"),
           ~ if (treatment == "NttE" && str_remove(cur_column(), "\\..*") %in% exotic_spp) 1 else .x),
    
    #case: tE 
    across(ends_with(".mid_state"),
           ~ if (treatment == "tE" && str_remove(cur_column(), "\\..*") %in% exotic_spp) 1 else .x)
    
    #there's also 'tE?' and 'NtE?' but we'll ignore that for now 
    
  ) %>%
  ungroup() 



#let's organize the dataframe columns 
Pryer_final <- Pryer_allinfo %>% 
  dplyr::select(c(water, treatment, 
                  contains('action'),
                  contains('mid_state'), 
                  contains('late_state'),
                  contains('outcome'), 
                  total.n, 
                  total.i, 
                  total.seeded, 
                  problem.row, 
                  Comments))

#Writeup :) ---- 
write.csv(Pryer_final, 'LOVE_updated_data_pryer.csv', row.names = F)
  