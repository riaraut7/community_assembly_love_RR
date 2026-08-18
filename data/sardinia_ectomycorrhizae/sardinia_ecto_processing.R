#Ria Raut, July 28th, 2026 
#Steps: 
#1.load both datasets 
#change names of the columns in your metadata 
#2.Go through the metadata and add environmental variables to the other dataset by sample 
#3. Go through the raw data and compile species pools 

#Basic housekeeping and loading ---- 
library(tidyverse)
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')

metadata <- read.csv('data/sardinia_ectomycorrhizae/metadata_sardinia.csv')
colnames(metadata)

raw_data <- read.csv('data/sardinia_ectomycorrhizae/sardinia_raw_assembly.csv')
colnames(raw_data)

#check for names in metadata that aren't in raw_data 
setdiff(metadata$camp, raw_data$Genere) #everything in metadata is in raw data 
#vice versa 
setdiff(raw_data$Genere, metadata$camp)
#There's stuff in raw data that isn't in the metadata 

