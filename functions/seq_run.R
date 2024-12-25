library(tidyverse)
library(lubridate)
library(nmecr)
library(dplyr)
source("../functions/saving_norm.R")

seq_run <- function(param, dataframe, tmy){
  
  #' @description This is a function defined to run sequential test (including sequential probability ratio test) 
  #' @usage seq_run(param, dataframe, tmy)
  #' @param param sequential test parameters
  #' @param dataframe measurement dataframe under randomized schedule
  #' @param tmy site tmy
  #' @return a list of dataframe as sequential test results
  #' 
  
  # Run SPRT
  # prepare sequential test dataset
  sprt_hourly <- dataframe %>%
    mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor())
  
  sprt_daily <- dataframe %>%
    group_by(datetime = floor_date(datetime, unit = "day")) %>%
    summarise(strategy = unique(strategy),
              power_ave = mean(eload, na.rm = TRUE),
              t_out = mean(t_out, na.rm = TRUE)) %>%
    ungroup()
  
  # prepare dataframe for analysis
  df_sprt <- sprt_daily %>%
    mutate(strategy = as.factor(strategy),
           strategy = recode_factor(strategy, "1" = "Baseline", "2" = "Intervention")) %>%
    filter(strategy %in% c(param$baseline, param$strategy)) %>%
    pivot_longer(cols = -c(datetime, strategy), 
                 names_to = "parameter", 
                 values_to = "value") %>%
    mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor(),
           value = value) %>%
    filter(str_detect(parameter, param$parameter)) %>%
    droplevels()
  
  
  # define list to store weekly means
  df_means <- list()
  
  # calculate weekly means
  for (i in 1:param$n_weeks) {
    
    # subset data by week
    df_means[[i]] <- df_sprt %>%
      filter(week <= i) %>%
      group_by(strategy, 
               parameter) %>%
      summarise(week = i,
                value_ave = mean(value, na.rm = TRUE),
                value_sd = sd(value, na.rm = TRUE), 
                .groups = "keep") %>%
      ungroup()
    
  }
  
  # combine list into df
  df_means <- bind_rows(df_means)
  
  # define lists to store stopping criteria results
  ## SPRT results
  sprt_res <- list()
  
  ## 80% independent variable
  overlap_base <- list()
  overlap_interv <- list()
  quantile_tmy <- quantile(tmy$temp, na.rm = TRUE, c(0, 1))
  
  ## annual saving estimation
  annual_saving <- list()
  
  # do sample
  for (i in seq(2, param$n_weeks)) {
    
    # subset data by week
    df_seq <- df_sprt %>%
      filter(week <= i) %>%
      droplevels()
    
    # Calculate overlapping temperature range
    overlap <- try(tibble("n_weeks" = i, 
                                overlap_base = sprt_hourly %>%
                                  filter(week <= i & strategy == 1) %>% 
                                  ol_est(., quantile_tmy)), silent = T)
    if (inherits(overlap, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      overlap_base[[i]] <- NULL
      
    } else {
      overlap_base[[i]] <- overlap
      
    }
    
    
    overlap <- try(tibble("n_weeks" = i, 
                                  overlap_interv = sprt_hourly %>%
                                    filter(week <= i & strategy == 2) %>% 
                                    ol_est(., quantile_tmy)), silent = T)
    
    if (inherits(overlap, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      overlap_interv[[i]] <- NULL
      
    } else {
      overlap_interv[[i]] <- overlap
      
    }
    
    # start estimation after 2 months for complete ftow profile
    if (i / 4 > 2 & i %% 4 == 0) {
      
      # Update user on the week of update
      # print(paste0("updating TOWT model in week ", i))
      annual_saving[[i]] <- tibble("n_weeks" = i,
                                   annual_saving = saving_norm(sprt_hourly %>% filter(week <= i), tmy))
      
    }
    
    # do test
    results_seq <- try(seq_ttest(x = value ~ strategy, 
                             data = df_seq,
                             d = 0.5, 
                             power = 0.9, 
                             alternative = "less", # greater
                             paired = FALSE,
                             verbose = TRUE), silent = T)
    
    if (inherits(results_seq, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      skip <- T
      
    } else {
      
      skip <- F
    }
    
    # calculate effect size
    results_ci <- try(effsize::cohen.d(d = df_seq$value,
                                   f = df_seq$strategy,
                                   conf.level = 0.90,
                                   paired = FALSE,
                                   na.rm = TRUE), silent = T)
    
    if (inherits(results_ci, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      skip <- T
      
    } else {
      
      skip <- F
    }
    
    # update user
    # print(paste0("Calculating effect size for week ", i))
    
    # bootstrap effect size
    results_bs <- try(bootES::bootES(data = df_seq, R = 1000, 
                                 contrast = c(param$baseline, param$strategy),
                                 data.col = "value", group.col = "strategy"), silent = T)
    
    if (inherits(results_bs, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      skip <- T
      
    } else {
      
      skip <- F
    }
    
    # extract test statistic
    if (skip){
      
      sprt_res[[i]] <- NULL
      
    } else {
      sprt_res[[i]] <- tibble("n_weeks" = i, 
                              "threshold_lower" = exp(results_seq@B_boundary_log),
                              "threshold_upper" = exp(results_seq@A_boundary_log),
                              "statistic" = results_seq@likelihood_ratio,
                              "decision" = results_seq@decision,
                              "cohens_d" = round(results_ci$estimate, digits = 2),
                              "ci_low" = round(results_ci$conf.int[[1]], digits = 2),
                              "ci_high" = round(results_ci$conf.int[[2]], digits = 2),
                              "ns_stat" = round(results_bs$t0, digits = 2),
                              "ns_ci_low" = round(results_bs$bounds[[1]], digits = 2),
                              "ns_ci_high" = round(results_bs$bounds[[2]], digits = 2))
      
      
    }
    
  }    
  
  # combine list and vectors into df
  annual_saving <- bind_rows(annual_saving)
  
  # work out when threshold is reached
  sprt_res <- bind_rows(sprt_res)
  
  sprt_res <- sprt_res %>%
    left_join(., annual_saving, by = "n_weeks") %>%
    mutate(flag = ifelse(str_detect(decision, "accept"), 1, 0))
  
  sprt_overlap_base <- bind_rows(overlap_base) %>% 
    left_join(., annual_saving, by = "n_weeks") %>% 
    mutate(flag = ifelse(overlap_base >= 0.8, 1, 0),
           flag = ifelse(flag == lag(flag, 1), 0, flag))
  
  sprt_overlap_interv <- bind_rows(overlap_interv) %>% 
    left_join(., annual_saving, by = "n_weeks") %>% 
    mutate(flag = ifelse(overlap_interv >= 0.8, 1, 0),
           flag = ifelse(flag == lag(flag, 1), 0, flag))
  
  return(list(annual_saving = annual_saving, 
              df_means = df_means,
              sprt_res = sprt_res, 
              sprt_overlap_base = sprt_overlap_base, 
              sprt_overlap_interv = sprt_overlap_interv))
  
}