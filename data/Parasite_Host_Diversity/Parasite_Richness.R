getwd()

library(readr) 
library(tidyverse)
library(stringr)

"remove useless columns"
DRE_Cleaned_1 <- DRE_Diversity_DF_25_March_2020 %>% 
  select(-c("S2", "S3", "S4", "S", "mpd.obs.z", "dummy_1", "dummy_2", "dummy_3", "dummy_4", "dummy_5", "resources", "diversity", "month"))


"remove all 2013 data"
DRE_Filtered <- DRE_Cleaned_1 %>%
    filter(Year != 2013)

"now that we removed 2013 get rid of year all together and another useless column"
DRE_Filtered <- DRE_Filtered %>% select(-c("Year", "S5m"))

"replace the plot and plot id columns w uniqueid"
DRE_PreFinal <- DRE_Filtered %>%
    select(-c("Plot", "PlotID", "S5m")) %>%
    rename(Plot = uniqueid)

"need to create columns for each species with entries of 0/1 and convert composition column data"
DRE_PreFinal <- DRE_PreFinal %>%
    mutate(
      ANVI = if_else(Div == 0, as.integer(spp == "ANVI"), as.integer(spp != "ANVI")),
      SEPA = if_else(Div == 0, as.integer(spp == "SEPA"), as.integer(spp != "SEPA")),
      TRFL = if_else(Div == 0, as.integer(spp == "TRFL"), as.integer(spp != "TRFL")),
      PAAN = if_else(Div == 0, as.integer(spp == "PAAN"), as.integer(spp != "PAAN")),
      SCIN = if_else(Div == 0, as.integer(spp == "SCIN"), as.integer(spp != "SCIN")),
      SOPI = if_else(Div == 0, as.integer(spp == "SOPI"), as.integer(spp != "SOPI"))
    )

