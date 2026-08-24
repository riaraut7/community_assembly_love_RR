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
