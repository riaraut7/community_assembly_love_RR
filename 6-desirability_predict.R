#Ria Raut 
#August 21st, 2026 

setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')

# Basic housekeeping. Run every time ----- 

try(dir.create(file.path(getwd(), 'outputs_desirability/figures'), recursive = TRUE))
try(dir.create(file.path(getwd(), 'outputs_desirability/statistical'), recursive = TRUE))
directory_string = file.path(getwd(), 'outputs_desirability/statistical')

# Load helpers and settings
DEBUG_MODE = FALSE
source('src/configs.R')
source('src/desirability_love.R')

# if on cluster
# CORES <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))

# Perform analyses
getwd()


# Getting shannon diversity for all datasets ---- 


  #wildflowers ---- 
set.seed(1)
data_wildflowers = read.csv('data/wildflowers/data_wildflowers.csv') %>% 
  mutate(Block = factor(Block)) %>%
  mutate(nitrogen.initial = factor(nitrogen.initial)) %>% 
  mutate(fungicide.initial = factor(fungicide.initial))

data_wildflowers <- data_wildflowers %>% 
  select(-c("Plot", 'Harvest', 'Calculated_herbivory', 'n_herb_occurences', 'Measured_infection', 
            'n_infect_occurences', "composition","Species_diversity", "Functional_composition" , "Sown_sla", "N", 
            "Sown_mpd_sla", "Notes", "Block"))

wildflowers_w_diversity <- add_shannon_diversity(data_wildflowers)
