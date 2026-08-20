getwd()

library(tidyverse)

# Remove unused columns
DRE_Cleaned_1 <- DRE_Diversity_DF_25_March_2020 %>% 
  select(-any_of(c(
    "S2", "S3", "S4", "S", "mpd.obs.z",
    "dummy_1", "dummy_2", "dummy_3", "dummy_4", "dummy_5",
    "resources", "diversity", "month"
  )))

# Remove 2013 data
DRE_Filtered <- DRE_Cleaned_1 %>%
  filter(Year != 2013)

# Remove Year, Plot, PlotID, and S5m, then rename uniqueid to Plot
DRE_PreFinal <- DRE_Filtered %>%
  select(-any_of(c("Year", "Plot", "PlotID", "S5m"))) %>%
  rename(Plot = uniqueid)

# Create 0/1 species columns
DRE_PreFinal <- DRE_PreFinal %>%
  mutate(
    ANVI = if_else(Div == 0, as.integer(spp == "ANVI"), as.integer(spp != "ANVI")),
    SEPA = if_else(Div == 0, as.integer(spp == "SEPA"), as.integer(spp != "SEPA")),
    TRFL = if_else(Div == 0, as.integer(spp == "TRFL"), as.integer(spp != "TRFL")),
    PAAN = if_else(Div == 0, as.integer(spp == "PAAN"), as.integer(spp != "PAAN")),
    SCIN = if_else(Div == 0, as.integer(spp == "SCIN"), as.integer(spp != "SCIN")),
    SOPI = if_else(Div == 0, as.integer(spp == "SOPI"), as.integer(spp != "SOPI"))
  )

# Creating Action and Outcome Columns
DRE_PreFinal_2 <- DRE_PreFinal %>%
  rename(
    ANVI.action = ANVI,
    SEPA.action = SEPA,
    TRFL.action = TRFL,
    PAAN.action = PAAN,
    SCIN.action = SCIN,
    SOPI.action = SOPI
  ) %>%
  mutate(
    ANVI.outcome = ANVI.action,
    SEPA.outcome = SEPA.action,
    TRFL.outcome = TRFL.action,
    PAAN.outcome = PAAN.action,
    SCIN.outcome = SCIN.action,
    SOPI.outcome = SOPI.action
  )

# Grouping Desirability Indexes, Then Get Rid of Them For Time 1
desirability <- c(
  'PlantS', 'S5', 'percov', 'Div', 'Res', 'divXres', 'n.plants'
)

DRE_PreFinal_3 <- DRE_PreFinal_2 %>%
  relocate(
    all_of(desirability), 
    .after = SOPI.outcome
  ) %>%
  mutate(
    across(all_of(desirability),
      ~ if_else(time == 1, NA_real_, .) 
    )
  )
# Final draft sent for review
  write_csv(DRE_PreFinal_3, "DRE_Final.csv")
  
# Deleting all 2012 rows to reduce redundancy, then removing the time column
  DRE_PreFinal_4 <- DRE_PreFinal_3 %>%
    filter(time != 1) %>%
    select(-any_of(c('time')))
  
# Getting rid of the 2014 part of plot ID and moving it to the first column spot
  DRE_Final_CSV <- DRE_PreFinal_4 %>%
    relocate(
      Plot,
      .before = Blk
    ) %>%
    mutate(
      Plot = str_remove_all(Plot, " 2014")
    )

# Writing Final CSV
  write_csv(DRE_Final_CSV, 'DRE_Final_CSV.csv')

    

