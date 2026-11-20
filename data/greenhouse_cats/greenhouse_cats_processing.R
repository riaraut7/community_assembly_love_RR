getwd()
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')

#housekeeping and data loading ---- 
library(tidyverse) 
cats_data_raw <- read.csv('data/Xgreenhouse_cats/greenhouse_cats_raw_data.csv')

new_cats <- cats_data_raw%>% 
  rename(species = Species.replicates.in.limestone.treatments) %>% 
  rename(species2 = Species.replicates.in.basalt.treatments)
colnames(new_cats)

#let's split into two dataset 
limestone_cats <- new_cats %>% 
  select(c(species, shade.limestone, sun.limestone))

basalt_cats <- new_cats %>% 
  select(c(species2, shade.basalt, sun.basalt))

#clean up the limestone data ---- 
limestone_cats <- limestone_cats %>%
  mutate(species = case_when(
    species == "B. gracilis" ~ "BOGR",
    species == "E. elymoides" ~ "ELEL",
    species == "F. arizonica" ~ "FEAR",
    species == "M. montana" ~ "MUMO",
    species == "P. fendleriana" ~ "POFE",
    TRUE ~ species
  ))

limestone_sun <- limestone_cats %>% select(c (species, sun.limestone)) %>% 
  mutate(plotID = 1:85)
class(limestone_sun$sun.limestone)
limestone_sun_wider <- limestone_sun %>% 
  pivot_wider(names_from = species, values_from = sun.limestone)
limestone_sun_wider <- limestone_sun_wider %>% 
  mutate(light.initial = 'high') %>% 
  mutate(soil.initial = 'limestone')


limestone_shade <- limestone_cats %>% select(c(species, shade.limestone)) %>% 
  mutate(plotID = 1:85)
limestone_shade_wider <- limestone_shade %>% 
  pivot_wider(names_from = species, values_from = shade.limestone)
limestone_shade_wider <- limestone_shade_wider %>% 
  mutate(light.initial = 'low') %>% 
  mutate(soil.initial = 'limestone')


#clean up the basalt data ---- 
basalt_cats <- basalt_cats %>%
  mutate(species2 = case_when(
    species2 == "B. gracilis" ~ "BOGR",
    species2 == "E. elymoides" ~ "ELEL",
    species2 == "F. arizonica" ~ "FEAR",
    species2 == "M. montana" ~ "MUMO",
    species2 == "P. fendleriana" ~ "POFE",
    TRUE ~ species2
  ))

#sun 
basalt_sun <- basalt_cats %>% select(c (species2, sun.basalt)) %>% 
  drop_na() %>% 
  mutate(plotID = 1:80)
basalt_sun_wider <- basalt_sun %>% 
  pivot_wider(names_from = species2, values_from = sun.basalt)
basalt_sun_wider <- basalt_sun_wider %>% 
  mutate(light.initial = 'high') %>% 
  mutate(soil.initial = 'basalt')

#shade 
basalt_shade <- basalt_cats %>% select(c (species2, shade.basalt)) %>% 
  drop_na() %>% 
  mutate(plotID = 1:80)
basalt_shade_wider <- basalt_shade %>% 
  pivot_wider(names_from = species2, values_from = shade.basalt)
basalt_shade_wider <- basalt_shade_wider %>% 
  mutate(light.initial = 'low') %>% 
  mutate(soil.initial = 'basalt')

#now combine all 
combined_data <- bind_rows(
  limestone_sun_wider, 
  limestone_shade_wider, 
  basalt_sun_wider, 
  basalt_shade_wider
)

final_combined <- combined_data %>% 
  select(-plotID) %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  filter(!if_all(where(is.numeric), ~ .x == 0))

write.csv(final_combined, 'data/Xgreenhouse_cats/data_greenhouse_cats.csv')

#1, wait....maybe abort? Bc we don't know which of these measured abundances come from the same pot..
#what is this experimental design bro... 


