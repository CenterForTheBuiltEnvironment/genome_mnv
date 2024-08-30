library(tidyverse)
library(lubridate)
library(stats)
library(dplyr)
source("../functions/create_temp_matrix.R")

model_pred <- function(dataframe, model){
  
  #' TOWT model prediction
  #' @description This is a function defined to use fitted TOWT model to make predictions on the given dataset. 
    #' Part of the code is extracted from nmecr package: https://github.com/kW-Labs/nmecr
  #' @param dataframe dataframe that stores the prediction info and results
  #' @param model fitted towt model
  #' @usage model_pred(dataframe, model)
  #' @return The function returns a dataframe containing the prediction datetime, towt(predicted energy) and eload (if contained in the input dataset)

  # get all info from input
  amod <- model$model_occupied 
  bmod <- model$model_unoccupied
  time <- dataframe$time
  temp <- dataframe$temp
  temp_knots <- model$model_input_options$calculated_temp_knots
  occ_info <- model$model_input_options$occupancy_info
  interval_minutes <- model$model_input_options$interval_minutes
  
  # prepare interval of week for training data
  minute_of_week <- (lubridate::wday(time) - 1) * 24 * 60 + 
    lubridate::hour(time) * 60 + lubridate::minute(time)
  
  interval_of_week <- 1 + floor(minute_of_week / interval_minutes)
  
  # which time intervals are 'occupied'?
  occ_intervals <- occ_info[occ_info[, 2] == 1, 1] 
  
  # create an occupancy vector for training dataset
  occ_vec <- rep(0, nrow(dataframe))
  for (i in 1 : length(occ_intervals)) {
    occ_vec[interval_of_week == occ_intervals[i]] <- 1
  }
  
  ok_occ <- occ_vec == 1
  ok_occ[is.na(ok_occ)] <- TRUE
  
  ftow <- factor(interval_of_week)
  
  temp_mat <- create_temp_matrix(temp, temp_knots)
  
  dframe <- data.frame(dataframe, ftow, temp_mat)
  
  # make predictions for occupied
  if (sum(ok_occ)){
    amod_towt <- stats::predict(amod, select(dframe[ok_occ, ], -eload))
    
    amod_results <- data.frame(dframe[ok_occ, ]$time, dframe[ok_occ, ]$eload, amod_towt)
    names(amod_results) <- c('time','eload','towt')
  } else {
    amod_results <- data.frame()
  }

  
  # make predictions for unoccupied
  if (sum(!ok_occ)){
    bmod_towt <- stats::predict(bmod, select(dframe[!ok_occ, ], -eload))
    
    bmod_results <- data.frame(dframe[!ok_occ, ]$time, dframe[!ok_occ, ]$eload, bmod_towt)
    names(bmod_results) <- c('time','eload','towt')
  } else {
    bmod_results <- data.frame()
  }
  
  # return results
  results <- rbind(amod_results, bmod_results) %>% 
    arrange(time)
  
  return(results)
}