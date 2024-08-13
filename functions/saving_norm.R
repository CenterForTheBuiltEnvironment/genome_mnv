library(tidyverse)
library(lubridate)
library(dplyr)
source("./functions/model_fit.R")
source("./functions/model_pred.R")
source("./functions/resample.R")

saving_norm <- function(dataframe, tmy, method){
  
  #' Energy saving calculation defined for TOWT prediction results (contains model_pred function)
  #' @description 
  #' by calculating the difference between baseline estimation and intervention estimation,
  #' and normalized on baseline for fractional savings.
  #'
  #' @param dataframe Dataframe that contains towt prediction info
  #' @param carryover_check Specify whether estimate for carryover purpose
  #' @param baseline_results If carryover_check == TRUE, baseline_results dataframe is needed
  #'
  #' @return A numerical percentage as annual fraction saving

  # load weather file  
  weather <- tmy %>%
    mutate(eload = NA) %>% 
    select(-site)
    
  # normal energy saving calculation
  towt_base <- dataframe %>%
    resample(.) %>%
    filter(strategy == 1) %>%
    select(-strategy, -week) %>%
    model_fit()
  
  towt_s2 <- dataframe %>%
    resample(.) %>%
    filter(strategy == 2) %>%
    select(-strategy, -week) %>%
    model_fit()
  
  base_results <- model_pred(weather, towt_base) %>%
    rename("datetime" = "time")
  
  interv_results <- model_pred(weather, towt_s2) %>%
    rename("datetime" = "time")
  
  
  saving_df <- data.frame(time = weather$time,
                          temp = weather$temp,
                          S2 = base_results$towt - interv_results$towt)
  
  saving_df_week <- saving_df %>%
    mutate(week = week(time)) %>%
    group_by(week) %>%
    summarise(S2_avg = mean(S2)) %>% 
    ungroup()
  
  FS <- 100*sum(saving_df$S2, na.rm = TRUE)/sum(base_results$towt)
  
  return(FS)
}