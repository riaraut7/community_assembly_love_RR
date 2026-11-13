#September 10th, 2026 
#Ria Raut 

#Basic housekeeping + loading your datasets ---- 
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')
library(data.table)
library(tidyverse)
library(stringr)

fn_outputs = dir('outputs_desirability',pattern="results.*\\.csv",full.names = TRUE)
df_all_raw = rbindlist(lapply(1:length(fn_outputs), function(i){
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  return(df_this)
}), 
fill = TRUE

)
all_des <- df_all_raw %>% 
  select(c (name, 
            training_size, 
            contains(".desirability"))
         )

#STEP 1: compile into direct and non-direct outcomes ---- 

#first filter by the training size you want -- maybe let's keep the lowest, 20? 
  #maybe keep 2-3 levels, actually? So filter out anything bigger than a certain amount 
#first pivot longer so that the des traits become row values (des_trait)
#then pivot wider and take values from name to making column 
#then add a column called des_trait_type and sort by if a trait is direct or indirect 

#then graph so your y axis is MAE values, x is direct vs indirect 
