#Ria Raut 
#January 23rd, 2026 

setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/LOVE_git_scripts/new_data_initial_processing/JenaExperiment_biomass_all')
getwd()

library(readr) 
library(tidyverse)
library(stringr)

#Cleaning 2002 initial action state data ---- 
#load stupid yumcky messy dataset 
jena_dom_2002 <- as.data.frame(
  read_tsv(
  "datasets/JenExp_biomass_2002_trial.txt",
  col_names = TRUE,
  na = c("", "NA")
))

#clean column titles, select out useless columns 
jena_2002_cleaned <- jena_dom_2002 %>% 
  rename_with(
    ~.x %>% 
      str_replace_all(' ', '_') %>% #replace all columns' spaces w underscores
      str_replace_all("_biom_\\[g/m\\*\\*2\\]", "")  %>% 
      str_replace_all("_biom_dw_\\[g/m\\*\\*2\\]", "") %>% 
      str_replace_all('_\\[m\\]', '') %>% 

      str_replace_all('._', '_') 
  )
jena_2002_cleaned <- jena_2002_cleaned[ ,-(5:6)] #taking out stupid date columns 
colnames(jena_2002_cleaned)

#you only want to keep the mean and not the 1st and second replicate, becaues youre converting all of this to presence/abs anyway 
#how many unique plots do you have? 
length(unique(jena_2002_cleaned$Experimenta_plot)) #206 -- that's how many rows you should have 

jena_2002_MEANS_ONLY <- jena_2002_cleaned %>% 
  filter(Repl == 'mean')

jena_2002_pruned <- jena_2002_MEANS_ONLY%>% 
  select(-c(Heigh_aboveg_min, Heigh_aboveg_max, Height, Repl, Sow_plant)) 
colnames(jena_2002_pruned)

#replace all -9999
replace_9999 <- jena_2002_pruned %>% 
  mutate(
    across(
      where(is.numeric), 
      ~if_else(.x == -9999, 0, .x)
    )
  )
summary(replace_9999)

jean_pruned_cleaned <- replace_9999 %>% 
  mutate(NEWSPS = Uni_plan_mat+Weeds) %>% 
  select(-c(Uni_plan_mat, Weeds))

#replace all numbers with 1 (metric of 'present )
presence_only <- jean_pruned_cleaned %>%
  mutate (
    across(
      where(is.numeric), #go through all columns that are numeric 
      ~if_else(is.na(.x), NA_real_, 1) #and if there's an NA, then keep the real NA, else replace with 1 
    )
  )

#replace every NA w 0 
presence_absence <- presence_only %>%
  mutate(across(
    where(is.numeric), #across all numeric columns 
    ~ replace_na(as.numeric(.x), 0) #replace NAs w 0? 
  ))

#And lastly, even though new species were harvested, they weren't initially planted. So let's change all 1's there to 0
presence_absence$NEWSPS <- 0

#NICE! That's the cleaned data, so now let's write it up 
write.csv(presence_absence, 'jena_2002_initial_state.csv')


#Cleaning 2015 final state data ---- 
jena_dom_2015 <- as.data.frame(
  read_tsv(
    "datasets/JenExp_biomass_2015_trial.txt",
    col_names = TRUE,
    na = c("", "NA")
  ))

#clean column titles, select out useless columns 
jena_dom_2015 <- jena_dom_2015 %>% 
  rename_with(
    ~.x %>% 
      str_replace_all(' ', '_') %>% #replace all columns' spaces w underscores
      str_replace_all("_biom_\\[g/m\\*\\*2\\]", "")  %>% 
      str_replace_all("_biom_dw_\\[g/m\\*\\*2\\]", "") %>% 
      str_replace_all('_\\[m\\]', '') %>% 
      
      str_replace_all('._', '_') 
  )

#I'm making the conscious decision here to keep all the dates and harvests in both may and august 
#to have an increased # of replicates. Let's instead keep date as a dif env 
jena_2015_by_months <- jena_dom_2015 %>% 
  mutate(
    month = jena_dom_2015[,5]
  )

jena_2015_by_months_cleaner <- jena_2015_by_months %>%
  mutate(
    month = case_when(
      month == as.Date("2015-05-27") ~ "May",
      month == as.Date("2015-08-25") ~ "August",
      month == as.Date("2015-08-26") ~ "August",
      TRUE ~ NA_character_
    )
  )


jena_2015_pruned <- jena_2015_by_months_cleaner%>% 
  select(-c(Heigh_aboveg_min, Heigh_aboveg_max, Height, Sow_plant, Dea_plan_mat)) 
colnames(jena_2015_pruned)


jena_2015_pruned <- jena_2015_pruned[ ,-(2:3)] #taking out stupid date columns 
colnames(jena_2015_pruned)

#replace all -9999
replace_9999 <- jena_2015_pruned %>% 
  mutate(
    across(
      where(is.numeric), 
      ~if_else(.x == -9999, 0, .x)
    )
  )
summary(replace_9999)

jena_2015_cleaned <- replace_9999 %>% 
  mutate(NEWSPS = Uni_plan_mat+Weeds) %>% 
  select(-c(Uni_plan_mat, Weeds))
colnames(jena_2015_cleaned)

#replace every NA w 0 
jena_2015_final <- jena_2015_cleaned %>%
  mutate(across(
    where(is.numeric), #across all numeric columns 
    ~ replace_na(as.numeric(.x), 0) #replace NAs w 0? 
  ))

final_2015_rightcolumns <- jena_2015_final %>% 
  select(Experimenta_plot, month, T_repens, T_pratense, P_trivialis, P_pratense, 
         G_pratense, D_glomerata, A_elatius, A_sylvestris, A_pratensis, NEWSPS)
colnames(final_2015_rightcolumns)

#let's write this up 
write.csv(final_2015_rightcolumns, 'jena_2015_final_state.csv')

#Let's join intial and final states ---- 

jena_2002 <- read.csv('jena_2002_initial_state.csv')
jena_2015 <- read.csv('jena_2015_final_state.csv')

#Let's change the names to include action and output 
jena_2002_renamed <- jena_2002 %>% 
  rename(
    TREP.action = T_repens,
    TPRA.action = T_pratense,
    PTRI.action = P_trivialis,
    PPRA.action = P_pratense,
    GPRA.action = G_pratense,
    DGLO.action = D_glomerata,
    AELA.action = A_elatius,
    ASYL.action = A_sylvestris,
    APRA.action = A_pratensis,
    NEWSPS.action = NEWSPS
    
  ) %>% 
  select(-X)
colnames(jena_2002_renamed)

jena_2015_renamed <- jena_2015 %>% 
  rename(
    TREP.outcome = T_repens,
    TPRA.outcome = T_pratense,
    PTRI.outcome = P_trivialis,
    PPRA.outcome = P_pratense,
    GPRA.outcome = G_pratense,
    DGLO.outcome = D_glomerata,
    AELA.outcome = A_elatius,
    ASYL.outcome= A_sylvestris,
    APRA.outcome = A_pratensis,
    NEWSPS.outcome = NEWSPS
    
  ) %>% 
  select(-X)

colnames(jena_2015_renamed)

final_dataset <- jena_2002_renamed %>% 
  left_join(jena_2015_renamed, by = 'Experimenta_plot')

#Okay, good lord I think it's ready 
#do I need to take out all rows that have an NA in them? 
final_data_no_na <- na.omit(final_dataset)

write.csv(final_data_no_na, 'data_jena_wildflowers.csv', row.names = F)

