# Block design for running different sampling intervals
# Code written by Aoyu Zou




#### Setup ####
# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, here, scales, slider, patchwork, # general
               broom, ggpmisc, ggpubr, # linear models
               RMV2.0, # TOWT model
               sprtt, effsize) # sequential testing

# turn off scientific notation
options(scipen = 999)

# set default theme for ggplot
theme_set(theme_minimal())

# define base ggplot theme
theme_update(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
             plot.caption = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", colour = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             axis.text = element_text(size = 10),
             strip.text = element_text(size = 10, color = "grey20", face = "bold"),
             strip.background = element_blank())





#### READ ####
function_path <- "../functions/"
output_path <- "../readfiles/"
source(paste0(function_path, "rand_seq.R"))





#### GENERATION ####
# adding random sampling schedules
block_params <- list(start_date = "2016-03-01",
                     n_weeks = 96,
                     n_seasons = 16, 
                     block_unit = 6, 
                     pool = 20)

# schedule generation
exclude <- data.frame(date = NA)

# consecutive sampling
schedule_2_design <- data.frame()
schedule_3_design <- data.frame()
schedule_7_design <- data.frame()
i <- 1

while (i <= block_params$pool){
  
  schedule <- data.frame()
  for (sample_start in as.list(seq(as.Date(block_params$start_date), as.Date(block_params$start_date) + weeks(block_params$n_weeks), by = block_params$block_unit * 7))){

    sample_end = sample_start + weeks(block_params$block_unit) - days(1)
    block <- rand_seq(sample_start, sample_end, 2, 2, exclude)
    schedule <- rbind(schedule, block$decision_data)
  }

  df_schedule_2 <- schedule %>%
    rename(datetime = date,
           !!str_glue("strategy_{i}") := strategy)

  if (nrow(schedule_2_design) == 0) {
      schedule_2_design <- df_schedule_2
    } else {
      schedule_2_design <- left_join(schedule_2_design, df_schedule_2, by = "datetime")
    }

  # 3-day consecutive sampling
  schedule <- data.frame()
  for (sample_start in as.list(seq(as.Date(block_params$start_date), as.Date(block_params$start_date) + weeks(block_params$n_weeks), by = block_params$block_unit * 7))){

    sample_end = sample_start + weeks(block_params$block_unit) - days(1)
    block <- rand_seq(sample_start, sample_end, 3, 2, exclude)
    schedule <- rbind(schedule, block$decision_data)
  }

  df_schedule_3 <- schedule %>%
    rename(datetime = date,
           !!str_glue("strategy_{i}") := strategy)

  if (nrow(schedule_3_design) == 0) {
    schedule_3_design <- df_schedule_3
  } else {
    schedule_3_design <- left_join(schedule_3_design, df_schedule_3, by = "datetime")
  }

  # weekly switching schedules
  schedule <- data.frame()
  for (sample_start in as.list(seq(as.Date(block_params$start_date), as.Date(block_params$start_date) + weeks(block_params$n_weeks), by = block_params$block_unit * 7))){
    
    sample_end = sample_start + weeks(block_params$block_unit) - days(1)
    block <- rand_seq(sample_start, sample_end, 7, 10, exclude)
    schedule <- rbind(schedule, block$decision_data)
  }
  
  df_schedule_7 <- schedule %>%
    rename(datetime = date,
           !!str_glue("strategy_{i}") := strategy)
  
  if (nrow(schedule_7_design) == 0) {
    schedule_7_design <- df_schedule_7
  } else {
    schedule_7_design <- left_join(schedule_7_design, df_schedule_7, by = "datetime")
  }
  
  print(paste0("finished ", i))
  
  i <- i + 1
}





#### WRITE ####
write.csv(schedule_2_design, file = paste0(output_path, "schedule_2.csv"), row.names = FALSE, na = "")
write.csv(schedule_3_design, file = paste0(output_path, "schedule_3.csv"), row.names = FALSE, na = "")
write.csv(schedule_7_design, file = paste0(output_path, "schedule_7.csv"), row.names = FALSE, na = "")
