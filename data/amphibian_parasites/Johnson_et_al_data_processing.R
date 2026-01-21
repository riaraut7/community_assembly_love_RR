#Johson et al., 2019 data processing 
#Ria Raut, Sept 11th 

library(tidyverse)

setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/LOVE_git_scripts/new_data/Johnson_et_al_2019_royal_society_data')

infections_data <- read.csv('experiment_infection_totals_UPLOAD.csv', stringsAsFactors = F)

#see if p/a of amphibians (and thus community?) influences p/a of parasite -- not quite the same as persistence 
#of species across time? 
#take out host count 
#Make it so that you're multipling the B, P, etc. host columns by the number of P, B, etc. in the 'treatment' column 
#these counts are your action columns 
#then the outcome cols are just...the Count? Change that title to parasite_count 

cleaned <- infections_data %>% 
  select(-c(HostCount, CommFreq, B, P, R, G, T))


#Make individual counts 
cleaned$B.action <- str_count(cleaned$Treatment, 'B')
cleaned$P.action <- str_count(cleaned$Treatment, 'P')
cleaned$R.action <- str_count(cleaned$Treatment, 'R')
cleaned$G.action <- str_count(cleaned$Treatment, 'G')
cleaned$T.action <- str_count(cleaned$Treatment, 'T')
cleaned$row_ID <- row_number(cleaned)

wider_cleaned <- cleaned %>% 
  pivot_wider(names_from = Parasite, values_from = Count) %>% 
  rename(Carriers_added = TotalAdded)

final <- wider_cleaned %>% 
  select( c (row_ID, 
             Treatment, 
             HostRichness, 
             Carriers_added,
             contains('action'), 
             Alaria, 
             Cephalogonimus, 
             Ribeiroia))

write.csv(final, 'LOVE_updated_Johnson2019_parasites.csv', row.names = F)
