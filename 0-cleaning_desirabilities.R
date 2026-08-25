#Ria Raut 
#August 24th, 2026 
#one-time use 
#Cleaning up all datasets and adding actual desirability columns 

#Basic housekeeping -----
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW/data')

library(tidyverse)

# Function to compute Shannon diversity for a numeric vector of abundances
shannon_index <- function(columns) {
  columns <- columns[!is.na(columns)]        # drop NAs
  columns <- columns[columns > 0]            # zeros contribute nothing, so let's ignore those 
  if (length(columns) == 0) return(NA_real_)
  
  p <- columns / sum(columns)
  -sum(p * log(p))
}
# applying the calculate shannon index function to an actual dataset -- previous function nested in here 
add_shannon_diversity <- function(dataset) {
  dataset %>%
    mutate(
      diversity.desirability = pmap_dbl(
        select(., contains(".outcome")),
        ~ shannon_index(c(...))
      )
    )
}


#amphibian_parasites ---- 
amphibian_parasites <- read.csv('amphibian_parasites/data_amphibian_parasites.csv', stringsAsFactors = T)
colnames(amphibian_parasites)
amphibian_parasites <- amphibian_parasites %>% 
  select(Treatment, B.action, P.action, R.action, G.action, T.action, 
         Alaria, Cephalogonimus, Ribeiroia)

# columns that define a "block" (everything that's identical across
# the rows you want to compact together) - adjust to match your actual names
group_cols <- c("Treatment", 
                "B.action", "P.action", "R.action", "G.action", "T.action")

parasite_cols <- c("Alaria", "Cephalogonimus", "Ribeiroia")
#taken from chat gpt, sorry 
df_compact <- amphibian_parasites %>%
  group_by(across(all_of(group_cols))) %>%
  group_modify(~ {
    vecs <- lapply(parasite_cols, function(col) {
      .x[[col]][!is.na(.x[[col]])]   # pull just the non-NA values, in row order
    })
    names(vecs) <- parasite_cols
    
    max_len <- max(lengths(vecs))
    vecs <- lapply(vecs, function(v) { length(v) <- max_len; v })  # pads end with NA
    
    as_tibble(vecs)
  }) %>%
  ungroup()

#let's also rename some things 
colnames(df_compact)
amp_parasites_new <- df_compact %>% 
  rename(ANAB.outcome = B.action, 
         TAGR.outcome = G.action, 
         PSRE.outcome = P.action,
         RACA.outcome = R.action, 
         TATO.outcome = T.action, 
         alaria.desirability = Alaria, 
         cephalo.desirability = Cephalogonimus, 
         ribei.desirability = Ribeiroia) 
amp_parasites_new <- amp_parasites_new %>% select(-Treatment)
amp_parasites_new <- add_shannon_diversity(amp_parasites_new)
write.csv(amp_parasites_new, 'amphibian_parasites/data_amphibian_parasites.csv', row.names = F)

#california_grasses ---- 
california_grasses <- read.csv('california_grasses/data_california_grasses.csv', stringsAsFactors = T)
colnames(california_grasses)
#change values for the water column, change name to water.initial 
#only keep total_seeded.desirability column 
#add shannon diversity index 

unique(california_grasses$water)
cal_grasses_treatments <- california_grasses %>% 
  mutate(water = recode(water, 
                        'No' = 'No', 
                        'Fall 2011 1.25cm' = 'Yes', 
                        'Winter 2013 1.45cm' = 'Yes', 
                        'Winter 2014 1.27cm' = 'Yes', 
                        'Winter 2014 1.27cm?' = 'Yes'
                        ))
cal_grasses_treatments <- cal_grasses_treatments %>% 
  rename('water.initial' = 'water',
         'total_seeded.desirability' = 'total.seeded')  %>% 
  select(-c(treatment, total.n, total.i, problem.row, Comments))
new_calgrasses_updated <- add_shannon_diversity(cal_grasses_treatments)

write.csv(new_calgrasses_updated, 'california_grasses/data_california_grasses.csv', row.names = F)

#ciliates ---- 
ciliates <- read.csv('ciliates/data_ciliates.csv', stringsAsFactors = T)
new_ciliates <- add_shannon_diversity(ciliates)

write.csv(new_ciliates, 'ciliates/data_ciliates.csv', row.names = F)

#fly_gut ---- 
fly_gut <- read.csv('fly_gut/data_fly_gut.csv', stringsAsFactors = T)
new_fly_gut <- add_shannon_diversity(fly_gut)
write.csv(new_fly_gut, 'fly_gut/data_fly_gut.csv', row.names = F)

#forest_trees ---- 
forest_trees <- read.csv('forest_trees/data_forest_trees.csv', stringsAsFactors = T)
new_forest_trees <- add_shannon_diversity(forest_trees)
write.csv(new_forest_trees, 'forest_trees/data_forest_trees.csv', row.names = F)

#fruit_flies ---- 
fruit_flies <- read.csv('fruit_flies/data_fruit_flies.csv', stringsAsFactors = T)
new_fruit_flies <- add_shannon_diversity(fruit_flies)
write.csv(new_fruit_flies, 'fruit_flies/data_fruit_flies.csv', row.names = F)

#grassland_annual_plants ---- 
grassland_annual_plants <- read.csv('grassland_annual_plants/data_grassland_annual_plants.csv', stringsAsFactors = T)
new_grassland_annual_plants <- add_shannon_diversity(grassland_annual_plants)
write.csv(new_grassland_annual_plants, 'grassland_annual_plants/data_grassland_annual_plants.csv', row.names = F)

#grassland_annual_plants_drought ---- 
grassland_annual_plants_drought <- read.csv('grassland_annual_plants_drought/data_grassland_annual_plants_drought.csv', stringsAsFactors = T)
new_grassland_annual_plants_drought <- add_shannon_diversity(grassland_annual_plants_drought)
write.csv(new_grassland_annual_plants_drought, 'grassland_annual_plants_drought/data_grassland_annual_plants_drought.csv', row.names = F)
#tbh I'm not sure what's happening with this dataset.... 

#grassland_diversity ---- 
grassland_diversity <- read.csv('grassland_diversity/data_grassland_diversity.csv', stringsAsFactors = T)
grassland_diversity <- grassland_diversity %>% 
  select(-plot.num) %>% 
  rename('weeded.initial' = 'weeded',
         'fertilized.initial' = 'fertilized')
new_grassland_diversity <- add_shannon_diversity(grassland_diversity)
write.csv(new_grassland_diversity, 'grassland_diversity/data_grassland_diversity.csv', row.names = F)

#human_gut ---- 
human_gut <- read.csv('human_and_mouse_gut/data_human_gut.csv', stringsAsFactors = T)
new_human_gut <- add_shannon_diversity(human_gut)
write.csv(new_human_gut, 'human_and_mouse_gut/data_human_gut.csv', row.names = F)

#mouse_gut ---- 
mouse_gut <- read.csv('human_and_mouse_gut/data_mouse_gut.csv', stringsAsFactors = T)
new_mouse_gut <- add_shannon_diversity(mouse_gut)
write.csv(new_mouse_gut, 'human_and_mouse_gut/data_mouse_gut.csv', row.names = F)


#jena_flowers ---- 
jena_wildflowers <- read.csv('jena_wildflowers/data_jena_wildflowers.csv', stringsAsFactors = T)
jena_wildflowers <- jena_wildflowers %>% select(-Experimenta_plot)
colnames(jena_wildflowers)
new_jena_wildflowers <- add_shannon_diversity(jena_wildflowers)
write.csv(new_jena_wildflowers, 'jena_wildflowers/data_jena_wildflowers.csv', row.names = F)

#parasite_host_diversity ---- 
#cols to keep: 
#S5: Parasite richness standardized to 5 random plants sampled within the plot (Desirability Variable)
#percov: percent cover within the plot (Desirability Varaibility)
parasite_host_diversity <- read.csv('parasite_host_diversity/data_parasite_host_diversity.csv', stringsAsFactors = T)
colnames(parasite_host_diversity)
new_parhost_diversity <- parasite_host_diversity %>% 
  rename(parasite_richness.desirability = S5, 
         percent_cover.desirability = percov)
new_parhost_diversity <- new_parhost_diversity %>% 
  select(ANVI.action, SEPA.action, TRFL.action, PAAN.action, SCIN.action, SOPI.action,
         ANVI.outcome, SEPA.outcome, TRFL.outcome, PAAN.outcome, SCIN.outcome, SOPI.outcome, 
         parasite_richness.desirability, percent_cover.desirability)
new_parhost_diversity <- add_shannon_diversity(new_parhost_diversity)
colnames(new_parhost_diversity)

write.csv(new_parhost_diversity, 'parasite_host_diversity/data_parasite_host_diversity.csv', row.names = F)

#prairie_plants ---- 
prairie_plants <- read.csv('prairie_plants/data_prairie_plants.csv', stringsAsFactors = T)
colnames(prairie_plants)
new_prairie_plants <- add_shannon_diversity(prairie_plants)
write.csv(new_prairie_plants, 'prairie_plants/data_prairie_plants.csv', row.names = F)

#soil_bacteria ---- 
soil_bacteria <- read.csv('soil_bacteria/data_soil_bacteria.csv', stringsAsFactors = T)
colnames(prairie_plants)
new_soil_bacteria <- add_shannon_diversity(soil_bacteria)
write.csv(new_soil_bacteria, 'soil_bacteria/data_soil_bacteria.csv', row.names = F)

#tree_colonization ----- 
#desirability cols to keep: 
#AvgHeight = Average height of individual trees in each plot (Desirability), 
#Trees 2015 = # of Individual trees/plot in 2015 (Desirability),
#rich_colherbs = Colonizing herbaceous plants in 2015 (Desirability) 

tree_colonization <- read.csv('tree_colonization/data_tree_colonization.csv', stringsAsFactors = T)
colnames(tree_colonization)

new_tree_colonization <- tree_colonization %>% 
  rename(average_height.desirability = AvgHeight, 
         earliest_establishment.desirability = Duration, 
         herb_colonization.desirability = rich_colherbs) 
new_tree_colonization <- new_tree_colonization %>% 
  select(c(
    Pesticide.initial, Nitrogen.initial, ANVI.outcome, PAAN.outcome,
    SCIN.outcome, SEPA.outcome, SOPI.outcome, TRFL.outcome, 
    average_height.desirability, earliest_establishment.desirability, herb_colonization.desirability
  ))

new_tree_colonization <- add_shannon_diversity(new_tree_colonization)
write.csv(new_tree_colonization, 'tree_colonization/data_tree_colonization.csv', row.names = F)

#wildflowers ---- 
wildflowers <- read.csv('wildflowers/data_wildflowers.csv', stringsAsFactors = T)
colnames(wildflowers)
 
#keep calculated_herbivory as herbivory.desirability 
#keep measured_infection as infection.desirability 
#keep sown_sla as direct_sown_sla.desirability
#take out plot and harvest 
#take out species diversity, species diversity, functional composition
new_wildflowers <- wildflowers %>% 
  rename(
    'herbivory.desirability' = 'Calculated_herbivory', 
    'infection.desirability' = 'Measured_infection', 
    'direct_sown_sla.desirability' = 'Sown_sla'
  )

colnames(new_wildflowers)

new_wildflowers <- new_wildflowers %>% 
  select( -c (Plot, Harvest, n_herb_occurences, n_infect_occurences, 
              Block, composition, Species_diversity, N, Functional_composition, Sown_mpd_sla, Notes))
new_wildflowers <- add_shannon_diversity(new_wildflowers)
write.csv(new_wildflowers, 'wildflowers/data_wildflowers.csv', row.names = F)
