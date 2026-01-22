#Ria Raut
#Jan 21st, 2026
#Heckman et al., 2017 og data - let's LOVE-ify this >:) 

setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/LOVE_git_scripts/new_data_initial_processing')
library(tidyverse)


fletcher_data <- read.csv('Heckman_et_al_2017_Ecol_society_america/All_DRE_cover_2015.csv')


#Making starting values for 2011 ---- 
fletcher_2012_only <- fletcher_data %>% 
  filter(year == 2012) %>%
  filter(blk == 1)
unique(fletcher_2012_only$species.id) #oup, def not just 6 
unique(fletcher_2012_only$month) #only measured once 

plot_number_slices <- fletcher_2012_only %>%
  group_by(plot.num) %>%
  slice(1) %>%
  ungroup()

pruned <- plot_number_slices %>% 
  select(plot.num, diversity, spp) %>% 
  mutate(value = 1)

pruned_pivoted <- pruned %>% 
  pivot_wider(names_from = spp, 
              names_glue = "{spp}.action",
              values_from =  value, 
              values_fill = 0)

#flipping values for polycultures 
pruned_flip_polyculture <- pruned_pivoted %>% 
  mutate(
  across(
    -c(plot.num, diversity),
    ~ if_else(diversity %in% c("Monoculture", "K"), .x, 1 - .x)
  )
)
final_2011 <- pruned_flip_polyculture %>% 
  rename(NEWSPS.action = BLAN.action) %>% 
  mutate(NEWSPS.action = 0)
#  select(-BLAN.action) 

#write this up! 
write.csv(final_2011, 'Heckman_et_al_2017_Ecol_society_america/grassland_diversity_intial_state.csv', row.names = F)

#Getting final composition + percov data ---- 
fletcher_2015_only <- fletcher_data %>% 
  filter(year == 2015) %>% 
  select(-c(plot, provenance, func.2, month, year))
unique(fletcher_2015_only$species.id)

og_species_list <- c('Setaria parviflora', 'Scutellaria integrifolia', 'Tridens flavus', 
                     'Andropogon virginicus', 'Packera anonyma','Solidago pinetorum')

non_og_summed <- fletcher_2015_only %>% 
  mutate(
    species_group = if_else(
      species.id %in% og_species_list, 
      species.id, 
      'non_og'
    )
  ) %>% 
  group_by(plotid, species_group, plot.num, blk, enemy, resources) %>% 
  summarise(
    percov = sum(percov, na.rm = TRUE), 
    .groups = 'drop'
)

#renaming things 
species_renamed <- non_og_summed %>%
  mutate(species_group = case_when(
    species_group == "Setaria parviflora" ~ "SEPA.outcome",
    species_group == "Tridens flavus" ~ "TRFL.outcome",
    species_group == "Scutellaria integrifolia" ~ "SCIN.outcome",
    species_group == "Andropogon virginicus" ~ "ANVI.outcome",
    species_group == "Packera anonyma" ~ "PAAN.outcome",
    species_group == "Solidago pinetorum" ~ "SOPI.outcome",
    species_group == "non_og" ~ "NEWSPS.outcome",
    TRUE ~ species_group # This keeps any other values (like "maybe") unchanged
  ))

pivoted_outcomes <- species_renamed%>%
  pivot_wider(names_from =  species_group, values_from = percov)

replaced_NAs <- pivoted_outcomes %>%
  mutate(across(
    ends_with(".outcome"),
    ~ replace_na(as.numeric(.x), 0)
  ))

write.csv(replaced_NAs, 'Heckman_et_al_2017_Ecol_society_america/grassland_diversity_outcomes.csv', row.names = F)

#Combining initial and final states into one dataset ---- 
initial_state <- read.csv('Heckman_et_al_2017_Ecol_society_america/grassland_diversity_intial_state.csv')
final_state <- read.csv('Heckman_et_al_2017_Ecol_society_america/grassland_diversity_outcomes.csv')

#change some column names 
final_state<- final_state %>% 
  rename(weeded = enemy, 
         fertilized = resources) 


final_dataset_almost <- final_state %>% 
  left_join(initial_state, by = 'plot.num')

final_dataset <- final_dataset_almost %>% 
  dplyr::select(c(#plotid, 
                  plot.num, 
                  #blk, 
                  contains('action'),
                  weeded, fertilized,
                  contains('outcome')
                  ))

#write it up!
write.csv(final_dataset, 'Heckman_et_al_2017_Ecol_society_america/data_grassland_diversity.csv', row.names = F)


