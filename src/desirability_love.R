#Ria Raut 
#August 25th, 2026 
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')
getwd()

#Parameters to define and libraries to load ----- 
directory_string <- 'C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW/outputs_desirability'
training_sizes_list <- c(20, 50, 80, 150, 400, 800, 1200, 1500)
library(randomForestSRC)
library(tidyverse)

#STEP 1a: Spitting the data and fitting an rf model ----- 

#subsect data -- takes in full dataset and size of training data, spits out training and testing data. Can contain replicates 
split_dataset <- function(data, training_size) {
  n <- nrow(data)
  
  if (training_size > n) {
    stop("training_size cannot be larger than the number of rows in the dataset")
  }
  
  train_indices <- sample(seq_len(n), size = training_size, replace = FALSE)
  
  training_data <- data[train_indices, ]
  testing_data <- data[-train_indices, ]
  
  list(training_data = training_data, testing_data = testing_data)
}
#actual function g function 
fit_rf_desirability <- function(training_dataset) { #fughh this doesn't work 
  #get outcomes cols 
  data_outcomes <- training_dataset %>% 
    dplyr::select(contains(".outcome"))
  dout_colnames <- names(data_outcomes)
  
  #get desirability cols 
  data_desirability <- training_dataset %>% 
    dplyr::select(contains(".desirability"))
  ddesir_colnames <- names(data_desirability)
  
  #define your formula 
  formula_rf_model <- formula(sprintf(
    "Multivar(%s) ~ %s", 
    paste(ddesir_colnames, collapse = ", "), 
    paste(dout_colnames, collapse = " + ")    
  ))
  
  #fit the rf model! 
  rf_outcome_desirability <- rfsrc(
    formula = formula_rf_model, 
    data = training_dataset,
    forest = TRUE,
    ntree = 500
  )
  
  return(rf_outcome_desirability)
}

#STEP 1b: testing the 1a functions ---- 
forest_trees <- read.csv('data/forest_trees/data_forest_trees.csv', stringsAsFactors = T)
forest_datasets <- split_dataset(forest_trees, 500)
training_dataset <- forest_datasets$training_data

forest_rf <- fit_rf_desirability(training_dataset)

#STEP 2a: Getting predictions, groundtruth vals, and mae ---- 
#get comparison tables -- both functions output a single table with outcome columns and desirabilities 
get_predictions_table <- function(rf_model, testing_dataset) {
  predictions_raw <- randomForestSRC::predict.rfsrc(rf_model, testing_dataset) 
  predicted_actual <- as.data.frame(predictions_raw$predicted) 
  colnames(predicted_actual) <- rf_model$yvar.names
  return(predicted_actual)
}

#rn it gets all desirability columns -- change it to take desirability, then select 
#for string(desirability) in the col name. Do that once you change the rf to fit per 
#des index, but dwai for now 
get_groundtruth_table <- function (testing_dataset) {
  data_desirability <- testing_dataset %>% 
    dplyr::select(contains(".desirability")) %>% 
    tidyr::drop_na()
  return(data_desirability)
} 

#intakes two datasets (in this case, just the desirability column) -- just use mean_absolute error casewise mean 
mean_absolute_error <- function(pred, obs) {
  pred <- pred[[1]]
  obs  <- obs[[1]]
  if (length(as.numeric(pred)) == length(as.numeric(obs)))
  {
    mae = mean(abs(as.numeric(pred) - as.numeric(obs)), na.rm=T)
  }
  else
  {
    mae = NA
  }
  
  return(mae)
}

#STEP 2b: testing the 2a functions ---- 
testing_dataset <- forest_datasets$testing_data

forest_predictions <- get_predictions_table(forest_rf, testing_dataset)
class(forest_predictions)
forest_ground_truth <- get_groundtruth_table(testing_dataset)

forest_mae <- mean_absolute_error(forest_predictions, forest_ground_truth)

#IT WORKED! 

#STEP 3a: Replicating per desirability trait and also n times ---- 
replicate_per_des_trait <- function(dataset, training_size) {
  
  desirability_cols <- names(dataset)[grepl(".desirability", names(dataset))]
  
  results <- purrr::map_dfr(desirability_cols, function(trait) {
    
    # keep this trait's column, drop the other .desirability columns
    other_traits <- setdiff(desirability_cols, trait)
    trait_dataset <- dataset %>% dplyr::select(-dplyr::all_of(other_traits))
    
    split_data <- split_dataset(trait_dataset, training_size) 
    training_data <- split_data$training_data 
    testing_data <- split_data$testing_data 
    
    specific_rf_model <- fit_rf_desirability(training_data) 
    
    predictions_table <- get_predictions_table(specific_rf_model, testing_data)
    groundtruth_table <- get_groundtruth_table(testing_data) 
    
    specific_mae <- mean_absolute_error(predictions_table, groundtruth_table)
    
    data.frame(trait = trait, mae = specific_mae)
  })
  wide_results <- tidyr::pivot_wider(results, names_from = trait, values_from = mae)
  return(wide_results)
}

replicate_n_times <- function(dataset, training_size, n_reps) {
  
  #getting a results table 
  results <- purrr::map_dfr(1:n_reps, function(i) {
    replicate_per_des_trait(dataset, training_size)
  })
  
  results <- results %>% 
    dplyr::mutate(training_size = training_size)
  
  return(results)
}


#STEP 3b: testing the 3a functions ---- 
forest_trees <- read.csv('data/forest_trees/data_forest_trees.csv', stringsAsFactors = T)
wildflowers <- read.csv('data/wildflowers/data_wildflowers.csv', stringsAsFactors = T)

forest_mae <- replicate_per_des_trait(forest_trees, 500)
wildflowers_mae <- replicate_per_des_trait(wildflowers, 500)
class(forest_mae)

replicated_forest_mae <- replicate_n_times(forest_trees, 500, 15)
replicated_wildflowers_mae <- replicate_n_times(wildflowers, 500, 15)

#STEP 4a: Replicate multiple times and across training sizes and make it into a csv ---- 

#has both functions for 3a nested within in 
replicate_all_across_training_sizes <- function(dataset, training_sizes, n_reps) {
  
  # only keep training sizes smaller than the dataset itself
  valid_sizes <- training_sizes[training_sizes < nrow(dataset)]
  
  results <- purrr::map_dfr(valid_sizes, function(size) {
    replicate_n_times(dataset, size, n_reps)
  })
  
  return(results)
}

#to get the dataset name for when you write up your results as a csv 
string_this <- function(x) {
  deparse(substitute(x))
}

#make your results into a csv! 
csv_this <- function (results_df, dataset_name) {
  if (!is.null(results_df)) {
    write.csv(
      results_df, file=sprintf('%s/results_%s.csv', directory_string, dataset_name), 
      row.names=FALSE)
  }
}  

#STEP 4b: testing 4a functions ----- 
forest_mae_table <- replicate_all_across_training_sizes(forest_trees, training_sizes_list, 10)

dataset_name <- string_this(forest_trees)
csv_this(forest_mae_table, dataset_name)

#STEP 5a: pulling all together into a final function ---- 
function_g <- function (dataset, training_sizes, n_reps, dataset_name) {
  results_df <- replicate_all_across_training_sizes(dataset, training_sizes, n_reps)
  csv_this(results_df, dataset_name)
}

#STEP 5b: testing 5a functions ----- 
tree_colonization <- read.csv('data/tree_colonization/data_tree_colonization.csv', stringsAsFactors =  T)

dataset_name <- string_this(tree_colonization)
function_g(tree_colonization, training_sizes_list, 10, dataset_name)
