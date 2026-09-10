#Ria Raut 
#August 21st, 2026 

setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')

# Basic housekeeping -- packages, variables, source scripts. Run every time ----- 
try(dir.create(file.path(getwd(), 'outputs_desirability'), recursive = TRUE))
source('src/desirability_love.R')

# Load helpers and settings
DEBUG_MODE = FALSE

#processing datasets ---- 
#forest_trees 
forest_trees <- read.csv('data/forest_trees/data_forest_trees.csv', stringsAsFactors =  T)
dataset_name <- string_this(forest_trees)
function_g(forest_trees, training_sizes_list, 10, dataset_name)

#tree colonization 
tree_colonization <- read.csv('data/tree_colonization/data_tree_colonization.csv', stringsAsFactors =  T)
dataset_name <- string_this(tree_colonization)
function_g(tree_colonization, training_sizes_list, 10, dataset_name)

#wildflowers 
wildflowers <- read.csv('data/wildflowers/data_wildflowers.csv', stringsAsFactors =  T)
dataset_name <- string_this(wildflowers)
function_g(wildflowers, training_sizes_list, 10, dataset_name)

#soil_bacteria
soil_bacteria <- read.csv('data/soil_bacteria/data_soil_bacteria.csv', stringsAsFactors =  T)
dataset_name <- string_this(soil_bacteria)
function_g(soil_bacteria, training_sizes_list, 10, dataset_name)

#prairie_plants 
prairie_plants <- read.csv('data/prairie_plants/data_prairie_plants.csv', stringsAsFactors =  T)
dataset_name <- string_this(prairie_plants)
function_g(prairie_plants, training_sizes_list, 10, dataset_name)

#parasite_host_diversity
parasite_host_diversity <- read.csv('data/parasite_host_diversity/data_parasite_host_diversity.csv', stringsAsFactors =  T)
dataset_name <- string_this(parasite_host_diversity)
function_g(parasite_host_diversity, training_sizes_list, 10, dataset_name)

#jena_wildflowers
jena_wildflowers <- read.csv('data/jena_wildflowers/data_jena_wildflowers.csv', stringsAsFactors =  T)
dataset_name <- string_this(jena_wildflowers)
function_g(jena_wildflowers, training_sizes_list, 10, dataset_name)

#grassland_diversity
grassland_diversity <- read.csv('data/grassland_diversity/data_grassland_diversity.csv', stringsAsFactors =  T)
dataset_name <- string_this(grassland_diversity)
function_g(grassland_diversity, training_sizes_list, 10, dataset_name)

#grassland_annual_plants_drought
grassland_annual_plants_drought <- read.csv('data/grassland_annual_plants_drought/data_grassland_annual_plants_drought.csv', stringsAsFactors =  T)
dataset_name <- string_this(grassland_annual_plants_drought)
function_g(grassland_annual_plants_drought, training_sizes_list, 10, dataset_name)

#fruit_flies 
fruit_flies <- read.csv('data/fruit_flies/data_fruit_flies.csv', stringsAsFactors =  T)
dataset_name <- string_this(fruit_flies)
function_g(fruit_flies, training_sizes_list, 10, dataset_name)

#fly_gut
fly_gut <- read.csv('data/fly_gut/data_fly_gut.csv', stringsAsFactors =  T)
dataset_name <- string_this(fly_gut)
function_g(fly_gut, training_sizes_list, 10, dataset_name)

#ciliates
ciliates <- read.csv('data/ciliates/data_ciliates.csv', stringsAsFactors =  T)
dataset_name <- string_this(ciliates)
function_g(ciliates, training_sizes_list, 10, dataset_name)






