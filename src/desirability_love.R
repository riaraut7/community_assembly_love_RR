#Ria Raut 
#August 25th, 2026 
setwd('C:/Users/riara/OneDrive/All Documents/UCBerk Personal research work/LOVE/community_assembly_love_RR-NEW')
getwd()

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
fit_rf_desirabilities <- function(training_dataset) { #fughh this doesn't work 
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

forest_rf <- fit_rf_desirabilities(training_dataset)

#STEP 2a: Getting predictions, groundtruth vals, and mae ---- 
#get comparison tables -- both functions output a single table with outcome columns and desirabilities 
get_predictions_table <- function(rf_model, testing_dataset, des_type) {
  predictions_raw <- randomForestSRC::predict.rfsrc(rf_model, testing_dataset) 
  predicted_actual <- as.data.frame(predictions_raw$predicted) 
  colnames(predicted_actual) <- rf_model$yvar.names
  return(predicted_actual)
}

#rn it gets all desirability columns -- change it to take desirability, then select 
#for string(desirability) in the col name. Do that once you change the rf to fit per 
#des index, but dwai for now 
get_groundtruth_table <- function (testing_dataset, des_type) {
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


#STEP 3a: Replicating multiple times and across training sizes ---- 
replicate_ten_times <- function(data, training_size) {
  
}
  
change_training_sizes <- function () { 
  }
#STEP 4a: Replicating 
  

#Ignore these functions for now and first let's test the top ---------- 
repeat_across_training_datasets <- function (data) { 
    outcome_cols <- grep("\\.outcome$", names(data), value = TRUE)
    combo_key <- do.call(paste, c(data[, outcome_cols, drop = FALSE], sep = "___"))
    n_unique <- length(unique(combo_key))

    n_list <- c(10, 14, 21, 30, 43, 62, 89)

    for (i in n_list[n_list < n_unique]) {
      # fit training RF model, test against testing dataset, get MAE using mean_absolute_error_casewise_mean 
      #add this to a collective dataset 
      #do this ten times 
    }
  
  #for each split, builds rf model and gets MAE_values 
  #replicates it ten times and makes a .matrix of desirability you tested, n used, n cols actually used, MAE for desirability 1 (with names), MAE for desirability 2 
  
  
}
get_final_results_csv <- function (dataset) {
  #get value from repeat_across_training_datasets and put into a csv  
}



