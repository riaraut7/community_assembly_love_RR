#Ria Raut 
#August 24th, 2026 
#one-time use 
#Cleaning up all datasets and adding actual desirability columns 

# Basic housekeeping -----
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

write.csv(amp_parasites_new, 'amphibian_parasites/data_amphibian_parasites.csv', row.names = F)

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
