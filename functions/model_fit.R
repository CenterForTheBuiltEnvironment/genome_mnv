library(tidyverse)
library(lubridate)
library(nmecr)
library(dplyr)
source("./functions/resample.R")
source("./functions/find_occ_unocc.R")

model_fit <- function(train_data){
  
  #' @description This is a function defined to fit TOWT model and predict energy consumption given datetime and outdoor temperature. 
  #' Part of the code is extracted from nmecr package: https://github.com/kW-Labs/nmecr/tree/master
  #' @usage model_fit(dataframe, method)
  #' @param dataframe Training dataset consisting of datetime, measured energy consumption and outdoor temperature
  #' @return The function will return a fitted energy prediction model and a dataframe
    
  model_input_options <- nmecr::assign_model_inputs(regression_type = "TOWT", occupancy_threshold = 0.65)
  
  df_towt <- train_data %>% 
    rename(time = datetime, 
           temp = t_out)
  
  # do baseline model to estimate occupancy
  towt_base <- nmecr::model_with_TOWT(training_data = df_towt,
                                      model_input_options = model_input_options)
  
  # Add occupancy info to model specification
  base_model_input_options <- map(model_input_options, ~.) %>% 
    append(list(find_occ_unocc(training_data = df_towt, model_input_options = towt_base$model_input_options))) %>% 
    set_names(., c(names(model_input_options), "occupancy_info"))
  
  # re-estimate baseline with occupancy
  towt_base <- nmecr::model_with_TOWT(training_data = df_towt,
                                      model_input_options = base_model_input_options)
  
  return(towt_base)
  
}