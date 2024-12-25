# Genome 2 M&V analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, # general
               broom, ggpmisc, ggpubr, # linear models
               sprtt, effsize) # sequential testing

# turn off scientific notation
options(scipen = 999)

# set default theme for ggplot
theme_set(theme_minimal())

# define base ggplot theme
theme_update(plot.title = element_text(size = 14, color = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 10, color = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
             plot.caption = element_text(size = 8, color = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", color = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             strip.text = element_text(size = 10, color = "grey20", face = "bold"),
             strip.background = element_blank())

# define strategy colors
ls_colors <- c("Baseline" = "#99d8c9",
               "Measured baseline" = "#99d8c9",
               "Adjusted baseline" = "#99d8c9",
               "Projected baseline\n(no change)" = "#1b9e77",
               "Intervention" = "#fdbb84",
               "Measured interv" = "#fdbb84",
               "True savings" = "black",
               "True normalized savings" = "black")

# load functions
function_path <- "../functions/"
source(paste0(function_path, "blocking.R"))
source(paste0(function_path, "model_fit.R"))
source(paste0(function_path, "model_pred.R"))
source(paste0(function_path, "saving_norm.R"))
source(paste0(function_path, "ol_est.R"))
source(paste0(function_path, "seq_run.R"))
source(paste0(function_path, "seq_plot.R"))
source(paste0(function_path, "err_plot.R"))
source(paste0(function_path, "prepost_plot.R"))
source(paste0(function_path, "cont_plot.R"))
source(paste0(function_path, "rand_seq.R"))

# define parameters
# section run control
run_params <- list(type = "variable", 
                   sprt = F, 
                   sprt_cont = F, 
                   interval = T, 
                   null = T)

# Adding intervention effect as advanced chiller operation
ctr_params <- list(peak_hours = 10:16,                      # accounts for peak hours
                   chwl_perc = 0.25,                        # % of chw electricity consumption (wrt building)
                   step_perc = 0.08,                        # % of chw electricity saving by raising 1 °C
                   conv_swt = 6,                            # baseline swt °C
                   weather_knots = c(15, 25),               # outdoor temperature reset steps
                   swt_knots = c(12, 7),                    # swt reset steps
                   coe_peak = 0.8,                          # coefficient adjustments for peak hours
                   coe_off = 1.2,                           # coefficient adjustments for off-peak hours
                   enable_temp = 8)                         # assume chiller operation starts


# adding random sampling schedules
block_params <- list(start_date = "2016-01-01",
                     n_weeks = 108,
                     n_seasons = 9, 
                     block_unit = 12)

# adding sprt criteria
sprt_param <- list(baseline = "Baseline",
                   strategy = "Intervention",
                   parameter = "power_ave",
                   label = "power",
                   n_weeks = 48)

cont_param <- list(baseline = "Baseline",
                   strategy = "Intervention",
                   parameter = "power_ave",
                   label = "power",
                   resamp = c(2, 8))





#### FUNCTIONS ####
# Function defined to extract legend using get_legend function
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function defined to read downloaded tmy files
get_tmy <- function(all_sites){
  
  df_tmy <- data.frame()
  
  for (site in all_sites){
    df <- read_csv(paste0(readfile_path, "tmy/", str_glue("{site}.epw")),
                   skip = 8, col_types = "ddddd-d---------------------------------",
                   col_names = c("year", "month", "day", "hour", "min", "tmy")) %>%
      mutate(year = 2017,
             time = ymd_h(paste(paste(year, month, day, sep = "-"), hour, sep = " ")),
             temp = tmy) %>%
      dplyr::select(time, temp) %>% 
      mutate(site = site)
    
    df_tmy <- rbind(df_tmy, df)
  }
  
  return(df_tmy)
}

# Function defined to adjust the plot scale
get_scale <- function(eload, range = 2){
  
  min_y <- mean(eload, na.rm = T) - range * sd(eload, na.rm = T)
  max_y <- mean(eload, na.rm = T) + range * sd(eload, na.rm = T)
  
  return(c(min_y, max_y))
}

# Function defined to add chwst reset intervention
run_reset <- function(df_baseline){
  
  mean <- mean(df_baseline$base_eload, na.rm = T) * ctr_params$chwl_perc
  
  grad <- (ctr_params$swt_knots[2] - ctr_params$swt_knots[1]) / 
    (ctr_params$weather_knots[2] - ctr_params$weather_knots[1])
  
  interc <- ctr_params$swt_knots[2] - (ctr_params$weather_knots[2] * grad)
  
  df_interv <- df_baseline %>% 
    mutate(swt = t_out * grad + interc, 
           chwl = mean,
           hour = hour(datetime)) %>% 
    mutate(swt = ifelse(swt > ctr_params$swt_knots[1], ctr_params$swt_knots[1], ifelse(swt < ctr_params$swt_knots[2], ctr_params$swt_knots[2], swt)), 
           temp_savings = ifelse(t_out >= ctr_params$enable_temp, (swt - ctr_params$conv_swt) * ctr_params$step_perc, 0), 
           time_adj = ifelse(hour %in% ctr_params$peak_hours, ctr_params$coe_peak, ctr_params$coe_off), 
           perc_savings = temp_savings * time_adj, 
           savings = chwl * perc_savings, 
           interv_eload = base_eload - savings) %>% 
    select(datetime, base_eload, interv_eload, t_out)
  
  return(df_interv)
}

# Function defined to interpolate NAs
run_interpo <- function(df_all){
  
  na_counts <- df_all %>%
    mutate(date = as.Date(timestamp)) %>%
    group_by(date) %>%
    summarize(na_hours = sum(is.na(eload)))
  
  # Filter out days with more than half of the hours having NAs
  valid_days <- na_counts %>%
    filter(na_hours <= 12) %>%
    pull(date)
  
  df_filtered <- df_all %>%
    filter(as.Date(timestamp) %in% valid_days)
  
  df_filtered <- df_filtered %>%
    mutate(across(c(eload, t_out), ~ zoo::na.approx(., na.rm = FALSE))) %>% 
    rename(datetime = timestamp, 
           base_eload = eload)
  
  return(df_filtered)
}

# Function defined to find end of blocking period
get_eob <- function(sprt_res, sprt_overlap_base, sprt_overlap_interv){
  
  sprt_check <- sprt_res %>% filter(flag == 1) %>% slice(1) %>% .$n_weeks
  bt_check <- sprt_overlap_base %>% filter(flag == 1) %>% .$n_weeks
  it_check <- sprt_overlap_interv %>% filter(flag == 1) %>% .$n_weeks
  eob <- seq(block_params$block_unit, sprt_param$n_weeks, by = block_params$block_unit)
  
  return(eob[eob >= max(c(sprt_check, bt_check, it_check))][1])
}

# Function defined to get sequential test results timeline
get_timeline <- function(sprt_res, sprt_overlap_base, sprt_overlap_interv){
  
  sprt <- sprt_res %>% filter(flag == 1) %>% slice(1) %>% .$n_weeks
  
  if (identical(sprt, integer(0))) {
    sprt <- 48
  }
  
  
  return(list(name = name,
              site = site, 
              sprt = sprt,
              base_temp = sprt_overlap_base %>% filter(flag == 1) %>% .$n_weeks,
              interv_temp = sprt_overlap_interv %>% filter(flag == 1) %>% .$n_weeks,
              eob = get_eob(sprt_res, sprt_overlap_base, sprt_overlap_interv), 
              final = sprt_param$n_weeks))
}

# Function defined to get sequential mean difference savings
get_mdsaving <- function(timeline, sprt_res){
  
  sprt_check <- timeline$sprt
  temp_check <- max(timeline$base_temp, timeline$interv_temp)
  eob <- timeline$eob
  final <- sprt_param$n_weeks
  
  return(list(name = name,
              site = site, 
              sprt = sprt_res %>% filter(n_weeks == sprt_check) %>% .$ns_stat,
              temp = sprt_res %>% filter(n_weeks == temp_check) %>% .$ns_stat,
              eob = sprt_res %>% filter(n_weeks == eob) %>% .$ns_stat,
              final = sprt_res %>% filter(n_weeks == final) %>% .$ns_stat))
}

# Function defined to get sequential normalized annual savings
get_nmsaving <- function(timeline, sprt_res){
  
  sprt_check <- timeline$sprt
  temp_check <- max(timeline$base_temp, timeline$interv_temp)
  eob <- timeline$eob
  final <- sprt_param$n_weeks
  
  return(list(name = name,
              site = site, 
              sprt = sprt_res %>% filter(n_weeks == sprt_check) %>% .$annual_saving,
              temp = sprt_res %>% filter(n_weeks == temp_check) %>% .$annual_saving,
              eob = sprt_res %>% filter(n_weeks == eob) %>% .$annual_saving,
              final = sprt_res %>% filter(n_weeks == final) %>% .$annual_saving))
}

# Function defined to get sequential fractional savings
get_frsaving <- function(timeline, dataframe){
  
  sprt_check <- timeline$sprt
  temp_check <- max(timeline$base_temp, timeline$interv_temp)
  eob <- timeline$eob
  final <- sprt_param$n_weeks
  
  return_l <- list()
  
  for (i in c(sprt_check, temp_check, eob, final)){
    
    df <- dataframe %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week <= i)
    
    # saving calculation as mean difference
    FS_rand <- (mean(df %>% filter(strategy == 1) %>% .$eload) -
                  mean(df %>% filter(strategy == 2) %>% .$eload)) /
      mean(df %>% filter(strategy == 1) %>% .$eload) * 100
    
    return_l <- c(return_l, FS_rand)
  }
  
  return(list(name = name,
              site = site, 
              sprt = return_l[[1]],
              temp = return_l[[2]],
              eob = return_l[[3]],
              final = return_l[[4]]))
  
}





#### READ DATA ####
readfile_path <- str_glue("../readfiles/{run_params$type}/")
summaryfigs_path <- str_glue("../figs/{run_params$type}/site_summary/")
combifigs_path <- str_glue("../figs/{run_params$type}/comb_analysis/")

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds"))
df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))
df_weather <- read_rds(paste0(readfile_path, "df_weather.rds"))
schedule_2_design <- read_csv(paste0(readfile_path, "../schedule_2.csv"))
schedule_3_design <- read_csv(paste0(readfile_path, "../schedule_3.csv"))
schedule_7_design <- read_csv(paste0(readfile_path, "../schedule_7.csv"))

all_sites <- df_energy %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types <- df_energy %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

all_names <- df_energy %>%
  select(name) %>%
  distinct(name)

df_tmy <- get_tmy(all_sites$site)






#### INDIVIDUAL ####
# storing results
FS_ref <- list()
MD_ref <- list()
model_acc <- list()
seq_timeline <- list()
seq_mdsaving <- list()
seq_nmsaving <- list()
seq_frsaving <- list()
seq_timeline_interval_2 <- list()
seq_timeline_interval_3 <- list()
seq_timeline_interval_7<- list()
seq_nmsaving_interval_2 <- list()
seq_frsaving_interval_2 <- list()
seq_nmsaving_interval_3 <- list()
seq_frsaving_interval_3 <- list()
seq_nmsaving_interval_7<- list()
seq_frsaving_interval_7<- list()
cont_saving <- list()
interval_saving <- list()
FS_tmy <- list()

# get site information
# for (n in 1:2){
for (n in 1:(nrow(all_names))){
  
  name <- all_names$name[n]
  
  site_info <- df_energy %>%
    filter(name == all_names$name[n]) %>%
    select(site, type) %>%
    distinct()
  
  site <- site_info$site
  
  ifelse(!dir.exists(file.path(str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}"))), dir.create(file.path(str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}"))), FALSE)
  sitefigs_path <- str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}")
  
  site_weather <- df_weather %>%
    filter(site == site_info$site) %>%
    select(timestamp, t_out) %>%
    group_by(timestamp) %>%
    summarise(t_out = mean(t_out)) %>%
    ungroup()
  
  site_tmy <- df_tmy %>% 
    filter(site == site_info$site)
  
  df_all <- df_energy %>%
    filter(name == all_names$name[n]) %>%
    select(timestamp, eload) %>%
    left_join(site_weather, by = "timestamp")
  
  # length check
  if (nrow(df_all) != (366 + 365) * 24){
    print("Incomplete/duplicate timestamp, please check")
  } else {
    print(paste0(name, " at ", site_info$site, " start"))
  }
  
  # Linear interpolation of baseline
  df_all <- df_all %>%
    run_interpo()
  
  plot_scale <- get_scale(df_all$base_eload)
  
  df_hourly_conv <- df_all %>%
    run_reset()
  
  # separate baseline and intervention
  df_base_conv <- df_hourly_conv %>%
    select(datetime,
           eload = base_eload,
           t_out) %>% 
    drop_na()
  
  df_interv_conv <- df_hourly_conv %>%
    select(datetime,
           eload = interv_eload,
           t_out) %>% 
    drop_na()
  
  # Check prediction accuracy
  towt_base <- df_base_conv %>%
    filter(datetime < as.Date("2017-01-01")) %>%
    model_fit()
  
  df_towt <- df_base_conv %>%
    filter(datetime < as.Date("2017-01-01")) %>%
    select(time = datetime,
           temp = t_out,
           eload)
  
  base_proj <- model_pred(df_towt, towt_base) %>%
    rename("datetime" = "time") %>%
    mutate(error = eload - towt)
  
  cv_rmse <- mean(sqrt(base_proj$error ^ 2)) / mean(base_proj$eload) * 100
  
  model_acc[[n]] <- tibble("name" = name,
                           "site" = site, 
                           "cvrmse" = cv_rmse)
  
  # base_proj %>%
  #   ggplot() +
  #   geom_point(aes(x = datetime, y = eload, color = "Measurement"), alpha = 0.2, size = 0.2) +
  #   geom_point(aes(x = datetime, y = towt, color = "Prediction"), alpha = 0.2, size = 0.2) +
  #   geom_smooth(aes(x = datetime, y = eload, color = "Measurement"), formula = y ~ x, method = "loess", linewidth = 0.7) +
  #   geom_smooth(aes(x = datetime, y = towt, color = "Prediction"), formula = y ~ x, method = "loess", linewidth = 0.7) +
  #   annotate(geom = "text",
  #            x = median(base_proj$datetime),
  #            y = median(base_proj$eload) + 2 * sd(base_proj$eload),
  #            label = paste0("CV(RMSE): ", round(cv_rmse, digits = 2), "%")) +
  #   scale_color_brewer(palette = "Set1") +
  #   scale_x_datetime(date_breaks = "2 months",
  #                    date_labels = "%b")  +
  #   scale_y_continuous(expand = c(0, 0),
  #                      breaks = breaks_pretty(n = 3),
  #                      labels = number_format(suffix = " kW")) +
  #   labs(x = NULL,
  #        y = NULL,
  #        color = NULL,
  #        title = "TOWT model prediction results") +
  #   theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
  #         legend.direction = "horizontal",
  #         legend.position = "bottom",
  #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  # 
  # ggsave(filename = "towt_acc.png", path = sitefigs_path, units = "in", height = 6, width = 10, dpi = 300)
  
  # TOWT baseline project for post retrofit period
  towt_base <- df_base_conv %>%
    filter(datetime < as.Date("2017-01-01")) %>%
    model_fit()
  
  df_towt <- df_interv_conv %>%
    filter(datetime >= as.Date("2017-01-01")) %>%
    select(time = datetime,
           temp = t_out,
           eload)
  
  base_proj <- model_pred(df_towt, towt_base) %>%
    rename("datetime" = "time")
  
  
  
  
  #### RANDOMIZATION ####
  schedule <- blocking(start_date = block_params$start_date,
                       n_weeks = block_params$n_weeks,
                       n_seasons = block_params$n_seasons,
                       seed = sample(1: 2^15, 1),
                       searches = 20,
                       jumps = 20,
                       treatments = 2,
                       consec = 1)
  
  # schedule summary
  schedule$weekday_summary
  
  df_schedule <- schedule$schedule %>%
    select(datetime = date,
           strategy)
  
  df_rand <- df_hourly_conv %>%
    left_join(df_schedule, by = "datetime") %>%
    fill(strategy, .direction = "down") %>%
    filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
    pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
    filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
    select(-eload_type) %>% 
    drop_na()
  
  err_plot(df_rand, df_base_conv, df_interv_conv)
  ggsave(filename = "sampling_err.png", path = sitefigs_path, units = "in", height = 9, width = 8, dpi = 300)
  
  # saving calculation as mean difference
  FS_true <- (mean(df_base_conv$eload) - mean(df_interv_conv$eload)) / mean(df_base_conv$eload) * 100
  FS_conv <- (mean(base_proj$towt) - mean(base_proj$eload)) / mean(base_proj$towt) * 100
  FS_rand <- (mean(df_rand %>% filter(strategy == 1) %>% .$eload) -
                mean(df_rand %>% filter(strategy == 2) %>% .$eload)) /
    mean(df_rand %>% filter(strategy == 1) %>% .$eload) * 100
  
  MD_true <- mean(df_base_conv$eload) - mean(df_interv_conv$eload)
  MD_conv <- mean(base_proj$towt) - mean(base_proj$eload)
  MD_rand <- mean(df_rand %>% filter(strategy == 1) %>% .$eload) -
    mean(df_rand %>% filter(strategy == 2) %>% .$eload)
  
  FS_ref[[n]] <- tibble("name" = name,
                        "site" = site, 
                        "ref_true" = FS_true,
                        "ref_conv" = FS_conv,
                        "ref_rand" = FS_rand)
  
  MD_ref[[n]] <- tibble("name" = name,
                        "site" = site, 
                        "ref_true" = MD_true,
                        "ref_conv" = MD_conv,
                        "ref_rand" = MD_rand)
  
  # p1 <- df_hourly_conv %>%
  #   mutate(savings = base_eload - interv_eload,
  #          year = as.factor(year(datetime))) %>%
  #   ggplot(aes(x = datetime, y = savings)) +
  #   geom_point(alpha = 0.2, size = 0.2) +
  #   geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.7) +
  #   facet_wrap(~year, scales = "free_x") +
  #   scale_x_datetime(date_breaks = "2 months",
  #                    date_labels = "%b")  +
  #   scale_y_continuous(expand = c(0, 0),
  #                      breaks = breaks_pretty(n = 3),
  #                      labels = number_format(suffix = " kW")) +
  #   labs(x = NULL,
  #        y = NULL,
  #        color = NULL,
  #        subtitle = "Calculated savings") +
  #   theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
  #         legend.direction = "horizontal",
  #         axis.text.x = element_blank(),
  #         legend.position = "bottom",
  #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  # 
  # p2 <- df_hourly_conv %>%
  #   mutate(year = as.factor(year(datetime))) %>%
  #   ggplot(aes(x = datetime, y = t_out)) +
  #   geom_point(alpha = 0.2, size = 0.2) +
  #   geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.7) +
  #   facet_wrap(~year, scales = "free_x") +
  #   scale_x_datetime(date_breaks = "2 months",
  #                    date_labels = "%b")  +
  #   scale_y_continuous(expand = c(0, 0),
  #                      breaks = breaks_pretty(n = 3),
  #                      labels = number_format(suffix = " °C")) +
  #   labs(x = NULL,
  #        y = NULL,
  #        color = NULL,
  #        subtitle = "Measured outdoor weather") +
  #   theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
  #         legend.direction = "horizontal",
  #         legend.position = "bottom",
  #         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  # 
  # ggarrange(p1, p2,
  #           ncol=1, nrow=2,
  #           labels = c("a)", "b)"),
  #           widths = c(1, 1),
  #           align = "h",
  #           common.legend = TRUE,
  #           legend="bottom") +
  #   plot_annotation(title = "Case study building intervention savings")
  # 
  # ggsave(filename = "savings_true.png", path = sitefigs_path, units = "in", height = 8, width = 10, dpi = 300)
  
  
  
  
  
  #### SPRT ####
  if (run_params$sprt){
    seq_res <- seq_run(sprt_param, df_rand, site_tmy)
    annual_saving <- seq_res$annual_saving
    df_means <- seq_res$df_means
    sprt_res <- seq_res$sprt_res
    sprt_overlap_base <- seq_res$sprt_overlap_base
    sprt_overlap_interv <- seq_res$sprt_overlap_interv
    
    # get true savings
    df_week <- df_base_conv %>% 
      select(datetime, 
             base_eload = eload) %>% 
      left_join(df_interv_conv, by = "datetime") %>% 
      mutate(savings = eload - base_eload) %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week <= sprt_param$n_weeks)
    
    true_saving <- list()
    
    for (i in c(2:sprt_param$n_weeks)){
      saving <- df_week %>% 
        filter(week <= i)
      
      true_saving[[i]] <- tibble("n_weeks" = i, 
                                 savings = mean(saving %>% .$savings))
    }
    
    true_saving <- bind_rows(true_saving)
    
    # get sequential test timeline
    eob <- get_eob(sprt_res, sprt_overlap_base, sprt_overlap_interv)
    seq_timeline[[n]] <- get_timeline(sprt_res, sprt_overlap_base, sprt_overlap_interv)
    
    # plot overall results
    seq_plot(df_means, sprt_res, sprt_overlap_base, sprt_overlap_interv, annual_saving, true_saving, eob)
    ggsave(filename = "overall_seq.png", path = sitefigs_path, units = "in", height = 9, width = 8, dpi = 300)
    
    # savings at timeline
    seq_mdsaving[[n]] <- get_mdsaving(seq_timeline[[n]], sprt_res)
    seq_nmsaving[[n]] <- get_nmsaving(seq_timeline[[n]], sprt_res)
    seq_frsaving[[n]] <- get_frsaving(seq_timeline[[n]], df_rand)
    
  }
  
  
  
  
  if (run_params$sprt_cont){
    
    cont_param$last_weeks <- block_params$n_weeks - eob
    cont_param$n_weeks <- block_params$n_weeks
    cont_param$start_date <- as.Date("2016-01-01") + weeks(eob)
    cont_param$n_seasons <- round(cont_param$last_weeks / block_params$block_unit)
    date_seq <- seq(from = cont_param$start_date, by = "day", length.out = cont_param$last_weeks * 7)
    strategy <- sample(c(1, 2), size = length(date_seq), replace = TRUE, prob = cont_param$resamp)
    
    df_schedule_cont <- data.frame(datetime = date_seq,
                                   strategy = strategy)
    
    df_rand_old <- df_rand %>% 
      filter(datetime < cont_param$start_date)
    
    df_rand_new <- df_hourly_conv %>%
      filter(datetime >= cont_param$start_date) %>% 
      left_join(df_schedule_cont, by = "datetime") %>%
      fill(strategy, .direction = "down") %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      select(-eload_type) %>% 
      drop_na()
    
    df_rand_cont <- rbind(df_rand_old, df_rand_new)
    
    rand_tmy <- saving_norm(df_rand_cont %>% mutate(week = NA), site_tmy)
    
    # saving calculation as mean difference
    rand_fs <- (mean(df_rand_cont %>% filter(strategy == 1) %>% .$eload) -
                  mean(df_rand_cont %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_cont %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md <- mean(df_rand_cont %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_cont %>% filter(strategy == 2) %>% .$eload)
    
    cont_saving[[n]] <- list("name" = name, 
                             "site" = site, 
                             "cont_md" = rand_md, 
                             "cont_fs" = rand_fs, 
                             "cont_tmy" = rand_tmy)
    
  }
  
  
  
  
  
  #### TMY ####
  rand_tmy <- saving_norm(df_rand %>% mutate(week = NA), site_tmy)
  
  df_conv_tmy <- df_hourly_conv %>% 
    mutate(strategy = ifelse(year(datetime) == "2016", 1, 2), 
           eload = ifelse(year(datetime) == "2016", base_eload, interv_eload), 
           week = NA) %>% 
    select(datetime, eload, strategy, t_out, week)
  
  conv_tmy <- saving_norm(df_conv_tmy %>% mutate(week = NA), site_tmy)
  
  
  FS_tmy[[n]] <- tibble("name" = name,
                        "site" = site, 
                        "rand" = rand_tmy, 
                        "conv" = conv_tmy)
  
  
  
  
  
  #### INTERVAL ####
  if (run_params$interval){
    
    seed <- sample(c(1:20), 1)
    
    # running on 2-day sampling interval 
    df_schedule_2 <- schedule_2_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch)) 
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_2, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_2 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                  mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_2 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_2 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    seq_res <- seq_run(sprt_param, df_rand_new, site_tmy)
    sprt_res <- seq_res$sprt_res
    sprt_overlap_base <- seq_res$sprt_overlap_base
    sprt_overlap_interv <- seq_res$sprt_overlap_interv
    seq_timeline_interval_2[[n]] <- get_timeline(sprt_res, sprt_overlap_base, sprt_overlap_interv)
    seq_nmsaving_interval_2[[n]] <- get_nmsaving(seq_timeline_interval_2[[n]], sprt_res)
    seq_frsaving_interval_2[[n]] <- get_frsaving(seq_timeline_interval_2[[n]], df_rand_new)
    
    # running on 3-day sampling interval 
    df_schedule_3 <- schedule_3_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch))
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_3, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_3 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                    mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_3 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_3 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    seq_res <- seq_run(sprt_param, df_rand_new, site_tmy)
    sprt_res <- seq_res$sprt_res
    sprt_overlap_base <- seq_res$sprt_overlap_base
    sprt_overlap_interv <- seq_res$sprt_overlap_interv
    seq_timeline_interval_3[[n]] <- get_timeline(sprt_res, sprt_overlap_base, sprt_overlap_interv)
    seq_nmsaving_interval_3[[n]] <- get_nmsaving(seq_timeline_interval_3[[n]], sprt_res)
    seq_frsaving_interval_3[[n]] <- get_frsaving(seq_timeline_interval_3[[n]], df_rand_new)
  
    
    # running with weekly interval
    df_schedule_7 <- schedule_7_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch))
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_7, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_7 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                    mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_7 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_7 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    seq_res <- seq_run(sprt_param, df_rand_new, site_tmy)
    sprt_res <- seq_res$sprt_res
    sprt_overlap_base <- seq_res$sprt_overlap_base
    sprt_overlap_interv <- seq_res$sprt_overlap_interv
    seq_timeline_interval_7[[n]] <- get_timeline(sprt_res, sprt_overlap_base, sprt_overlap_interv)
    seq_nmsaving_interval_7[[n]] <- get_nmsaving(seq_timeline_interval_7[[n]], sprt_res)
    seq_frsaving_interval_7[[n]] <- get_frsaving(seq_timeline_interval_7[[n]], df_rand_new)
    
    interval_saving[[n]] <- list("name" = name,
                                 "site" = site,
                                 "interval_md_2" = rand_md_2, 
                                 "interval_fs_2" = rand_fs_2, 
                                 "interval_tmy_2" = rand_tmy_2, 
                                 "interval_md_3" = rand_md_3, 
                                 "interval_fs_3" = rand_fs_3, 
                                 "interval_tmy_3" = rand_tmy_3, 
                                 "interval_md_7" = rand_md_7, 
                                 "interval_fs_7" = rand_fs_7, 
                                 "interval_tmy_7" = rand_tmy_7)
  }
  
  print(paste0("Finished: ", n, "/", nrow(all_names)))
  
} 
  




#### NULL ####
if (run_params$null) {
  
  FS_ref_null <- list()
  MD_ref_null <- list()
  FS_tmy_null <- list()
  seq_timeline_null <- list()
  seq_mdsaving_null <- list()
  seq_nmsaving_null <- list()
  seq_frsaving_null <- list()
  interval_null <- list()
  
  failed_sprt <- 0
  
  # for (n in 1:10){
  for (n in 1:(nrow(all_names))){
    
    name <- all_names$name[n]
    
    site_info <- df_energy %>%
      filter(name == all_names$name[n]) %>%
      select(site, type) %>%
      distinct()
    
    site <- site_info$site
    
    ifelse(!dir.exists(file.path(str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}"))), dir.create(file.path(str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}"))), FALSE)
    sitefigs_path <- str_glue("../figs/{run_params$type}/site_analysis/{site}/{name}")
    
    site_weather <- df_weather %>%
      filter(site == site_info$site) %>%
      select(timestamp, t_out) %>%
      group_by(timestamp) %>%
      summarise(t_out = mean(t_out)) %>%
      ungroup()
    
    site_tmy <- df_tmy %>% 
      filter(site == site_info$site)
    
    df_all <- df_energy %>%
      filter(name == all_names$name[n]) %>%
      select(timestamp, eload) %>%
      left_join(site_weather, by = "timestamp")
    
    # length check
    if (nrow(df_all) != (366 + 365) * 24){
      print("Incomplete/duplicate timestamp, please check")
    } else {
      print(paste0(name, " at ", site_info$site, " start"))
    }
    
    # Linear interpolation of baseline
    df_all <- df_all %>%
      run_interpo()
    
    plot_scale <- get_scale(df_all$base_eload)
    
    df_hourly_conv <- df_all %>%
      mutate(interv_eload = base_eload)
    
    # separate baseline and intervention
    df_base_conv <- df_hourly_conv %>%
      select(datetime,
             eload = base_eload,
             t_out) %>% 
      drop_na()
    
    df_interv_conv <- df_hourly_conv %>%
      select(datetime,
             eload = interv_eload,
             t_out) %>% 
      drop_na()
    
    # Check prediction accuracy
    towt_base <- df_base_conv %>%
      filter(datetime < as.Date("2017-01-01")) %>%
      model_fit()
    
    df_towt <- df_base_conv %>%
      filter(datetime < as.Date("2017-01-01")) %>%
      select(time = datetime,
             temp = t_out,
             eload)
    
    base_proj <- model_pred(df_towt, towt_base) %>%
      rename("datetime" = "time")
    
    
    # TOWT baseline project for post retrofit period
    towt_base <- df_base_conv %>%
      filter(datetime < as.Date("2017-01-01")) %>%
      model_fit()
    
    df_towt <- df_interv_conv %>%
      filter(datetime >= as.Date("2017-01-01")) %>%
      select(time = datetime,
             temp = t_out,
             eload)
    
    base_proj <- model_pred(df_towt, towt_base) %>%
      rename("datetime" = "time")
    
    # Repeat randomization
    schedule <- blocking(start_date = block_params$start_date,
                         n_weeks = block_params$n_weeks,
                         n_seasons = block_params$n_seasons,
                         seed = sample(1: 2^15, 1),
                         searches = 20,
                         jumps = 20,
                         treatments = 2,
                         consec = 1)
    
    # schedule summary
    schedule$weekday_summary
    
    df_schedule <- schedule$schedule %>%
      select(datetime = date,
             strategy)
    
    df_rand <- df_hourly_conv %>%
      left_join(df_schedule, by = "datetime") %>%
      fill(strategy, .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      select(-eload_type) %>% 
      drop_na()
    
    err_plot(df_rand, df_base_conv, df_interv_conv)
    ggsave(filename = "sampling_err_null.png", path = sitefigs_path, units = "in", height = 9, width = 8, dpi = 300)
    
    # saving calculation as mean difference
    FS_true <- (mean(df_base_conv$eload) - mean(df_interv_conv$eload)) / mean(df_base_conv$eload) * 100
    FS_conv <- (mean(base_proj$towt) - mean(base_proj$eload)) / mean(base_proj$towt) * 100
    FS_rand <- (mean(df_rand %>% filter(strategy == 1) %>% .$eload) -
                  mean(df_rand %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand %>% filter(strategy == 1) %>% .$eload) * 100
    
    MD_true <- mean(df_base_conv$eload) - mean(df_interv_conv$eload)
    MD_conv <- mean(base_proj$towt) - mean(base_proj$eload)
    MD_rand <- mean(df_rand %>% filter(strategy == 1) %>% .$eload) -
      mean(df_rand %>% filter(strategy == 2) %>% .$eload)
    
    FS_ref_null[[n]] <- tibble("name" = name,
                          "site" = site, 
                          "ref_true" = FS_true,
                          "ref_conv" = FS_conv,
                          "ref_rand" = FS_rand)
    
    MD_ref_null[[n]] <- tibble("name" = name,
                          "site" = site, 
                          "ref_true" = MD_true,
                          "ref_conv" = MD_conv,
                          "ref_rand" = MD_rand)
    
    # get true savings
    df_week <- df_base_conv %>% 
      select(datetime, 
             base_eload = eload) %>% 
      left_join(df_interv_conv, by = "datetime") %>% 
      mutate(savings = eload - base_eload) %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week <= sprt_param$n_weeks)
    
    true_saving_null <- list()
    
    for (i in 2:sprt_param$n_weeks){
      saving <- df_week %>% 
        filter(week <= i)
      
      true_saving_null[[i]] <- tibble("n_weeks" = i, 
                                       savings = mean(saving %>% .$savings))
    }
    
    true_saving_null <- bind_rows(true_saving_null)
    
    seq_res <- try(seq_run(sprt_param, df_rand, site_tmy), silent = TRUE)  
    
    if (inherits(seq_res, "try-error")) {
      
      message("An error occurred. Skipping this part...")
      failed_sprt <- failed_sprt + 1
      
    } else {
      
      annual_saving <- seq_res$annual_saving
      df_means <- seq_res$df_means
      sprt_res <- seq_res$sprt_res
      sprt_overlap_base <- seq_res$sprt_overlap_base
      sprt_overlap_interv <- seq_res$sprt_overlap_interv
      
      # get sequential test timeline
      eob <- get_eob(sprt_res, sprt_overlap_base, sprt_overlap_interv)
      seq_timeline_null[[n]] <- get_timeline(sprt_res, sprt_overlap_base, sprt_overlap_interv)
      
      # plot overall results
      seq_plot(df_means, sprt_res, sprt_overlap_base, sprt_overlap_interv, annual_saving, true_saving_null, eob)
      ggsave(filename = "overall_seq_null.png", path = sitefigs_path, units = "in", height = 9, width = 8, dpi = 300)
      
      # savings at timeline
      seq_mdsaving_null[[n]] <- get_mdsaving(seq_timeline_null[[n]], sprt_res)
      seq_nmsaving_null[[n]] <- get_nmsaving(seq_timeline_null[[n]], sprt_res)
      seq_frsaving_null[[n]] <- get_frsaving(seq_timeline_null[[n]], df_rand)
      
    }
    
    # TMY
    rand_tmy <- saving_norm(df_rand %>% mutate(week = NA), site_tmy)
    
    df_conv_tmy <- df_hourly_conv %>% 
      mutate(strategy = ifelse(year(datetime) == "2016", 1, 2), 
             eload = ifelse(year(datetime) == "2016", base_eload, interv_eload), 
             week = NA) %>% 
      select(datetime, eload, strategy, t_out, week)
    
    conv_tmy <- saving_norm(df_conv_tmy %>% mutate(week = NA), site_tmy)
    
    
    FS_tmy_null[[n]] <- tibble("name" = name,
                                "site" = site, 
                                "rand" = rand_tmy, 
                                "conv" = conv_tmy)
    
    # running on 2-day sampling interval 
    seed <- sample(c(1:20), 1)
    
    # running on 2-day sampling interval 
    df_schedule_2 <- schedule_2_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch))
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_2, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_2 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                    mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_2 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_2 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    # running on 3-day sampling interval 
    df_schedule_3 <- schedule_3_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch))
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_3, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_3 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                    mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_3 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_3 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    # running with weekly interval
    df_schedule_7 <- schedule_7_design %>% 
      select(datetime, 
             strategy = str_glue("strategy_{seed}")) %>% 
      mutate(switch = ifelse(strategy != lag(strategy), T, F)) %>% 
      mutate(switch = ifelse(row_number() == 1, F, switch))
    
    df_rand_new <- df_hourly_conv %>%
      left_join(df_schedule_7, by = "datetime") %>%
      fill(c(strategy, switch), .direction = "down") %>%
      filter(datetime <= as.Date("2016-01-01") + weeks(block_params$n_weeks)) %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      filter(switch == 0) %>% 
      select(-c(eload_type, switch)) %>% 
      drop_na()
    
    # saving calculation as mean difference
    rand_fs_7 <- (mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) -
                       mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) * 100
    
    # savings at timeline
    rand_md_7 <- mean(df_rand_new %>% filter(strategy == 1) %>% .$eload) - mean(df_rand_new %>% filter(strategy == 2) %>% .$eload)
    
    rand_tmy_7 <- saving_norm(df_rand_new %>% mutate(week = NA), site_tmy)
    
    interval_null[[n]] <- list("name" = name,
                               "site" = site,
                               "interval_md_2" = rand_md_2, 
                               "interval_fs_2" = rand_fs_2, 
                               "interval_tmy_2" = rand_tmy_2, 
                               "interval_md_3" = rand_md_3, 
                               "interval_fs_3" = rand_fs_3, 
                               "interval_tmy_3" = rand_tmy_3, 
                               "interval_md_7" = rand_md_7, 
                               "interval_fs_7" = rand_fs_7, 
                               "interval_tmy_7" = rand_tmy_7)

    print(paste0("Finished: ", n, "/", nrow(all_names)))
  
  }

}






#### BIND ####
# savings calculation
df_seq_FS <- bind_rows(seq_frsaving) %>%
  pivot_longer(-c(name, site), names_to = "seq", values_to = "FS")

df_MD <- bind_rows(MD_ref) %>%
  pivot_longer(-c(name, site), names_to = "type", values_to = "savings") %>%
  separate(type, into = c("scenario", "method"), sep = "_")

df_FS <- bind_rows(FS_ref) %>%
  pivot_longer(-c(name, site), names_to = "type", values_to = "savings") %>%
  separate(type, into = c("scenario", "method"), sep = "_")

df_sprt_all <- bind_rows(seq_mdsaving) %>%
  pivot_longer(-c(name, site), names_to = "seq", values_to = "sprt") %>%
  left_join(bind_rows(seq_nmsaving) %>%
              pivot_longer(-c(name, site), names_to = "seq", values_to = "annual"),
            by = c("name", "site", "seq")) %>%
  left_join(bind_rows(seq_timeline) %>%
              mutate(temp = pmax(base_temp, interv_temp)) %>%
              select(-c(base_temp, interv_temp)) %>%
              pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
            by = c("name", "site", "seq"))

df_model_acc <- bind_rows(model_acc) 

df_FS_tmy <- bind_rows(FS_tmy)

# store savings and timeline
write_rds(df_sprt_all, paste0(readfile_path, "df_sprt_all.rds"), compress = "gz")
write_rds(df_seq_FS, paste0(readfile_path, "df_seq_FS.rds"), compress = "gz")
write_rds(df_MD, paste0(readfile_path, "df_MD.rds"), compress = "gz")
write_rds(df_FS, paste0(readfile_path, "df_FS.rds"), compress = "gz")
write_rds(df_model_acc, paste0(readfile_path, "df_model_acc.rds"), compress = "gz")
write_rds(df_FS_tmy, paste0(readfile_path, "df_FS_tmy.rds"), compress = "gz")


if (run_params$sprt_cont){
  
  df_cont <- bind_rows(cont_saving)
  write_rds(df_cont, paste0(readfile_path, "df_cont.rds"), compress = "gz")
  
}

if (run_params$interval){
  
  seq_FS_interval_2 <- bind_rows(seq_frsaving_interval_2) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "FS") %>% 
    left_join(bind_rows(seq_timeline_interval_2) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 2)
  
  seq_FS_interval_3 <- bind_rows(seq_frsaving_interval_3) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "FS") %>% 
    left_join(bind_rows(seq_timeline_interval_3) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 3)
  
  seq_FS_interval_7 <- bind_rows(seq_frsaving_interval_7) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "FS") %>% 
    left_join(bind_rows(seq_timeline_interval_7) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 7)
  
  seq_nm_interval_2 <- bind_rows(seq_nmsaving_interval_2) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "annual") %>% 
    left_join(bind_rows(seq_timeline_interval_2) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 2)
  
  seq_nm_interval_3 <- bind_rows(seq_nmsaving_interval_3) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "annual") %>% 
    left_join(bind_rows(seq_timeline_interval_3) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 3)
  
  seq_nm_interval_7 <- bind_rows(seq_nmsaving_interval_7) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "annual") %>% 
    left_join(bind_rows(seq_timeline_interval_7) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq")) %>% 
    mutate(interval = 7)
  
  df_seq_interval_FS <- bind_rows(seq_FS_interval_2, seq_FS_interval_3, seq_FS_interval_7)
  df_seq_interval_nm <- bind_rows(seq_nm_interval_2, seq_nm_interval_3, seq_nm_interval_7)
  
  df_interval <- bind_rows(interval_saving)
  
  write_rds(df_interval, paste0(readfile_path, "df_interval.rds"), compress = "gz")
  write_rds(df_seq_interval_FS, paste0(readfile_path, "df_seq_interval_FS.rds"), compress = "gz")
  write_rds(df_seq_interval_nm, paste0(readfile_path, "df_seq_interval_nm.rds"), compress = "gz")
  
}

if (run_params$null){
  
  df_seq_FS_null <- bind_rows(seq_frsaving_null) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "FS")
  
  df_null_all <- bind_rows(seq_mdsaving_null) %>%
    pivot_longer(-c(name, site), names_to = "seq", values_to = "sprt") %>%
    left_join(bind_rows(seq_nmsaving_null) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "annual"),
              by = c("name", "site", "seq")) %>%
    left_join(bind_rows(seq_timeline_null) %>%
                mutate(temp = pmax(base_temp, interv_temp)) %>%
                select(-c(base_temp, interv_temp)) %>%
                pivot_longer(-c(name, site), names_to = "seq", values_to = "n_weeks"),
              by = c("name", "site", "seq"))
  
  df_MD_null <- bind_rows(MD_ref_null) %>%
    pivot_longer(-c(name, site), names_to = "type", values_to = "savings") %>%
    separate(type, into = c("scenario", "method"), sep = "_")
  
  df_FS_null <- bind_rows(FS_ref_null) %>%
    pivot_longer(-c(name, site), names_to = "type", values_to = "savings") %>%
    separate(type, into = c("scenario", "method"), sep = "_")
  
  df_interval_null <- bind_rows(interval_null)
  df_FS_tmy_null <- bind_rows(FS_tmy_null)
  
  write_rds(df_seq_FS_null, paste0(readfile_path, "df_seq_FS_null.rds"), compress = "gz")
  write_rds(df_null_all, paste0(readfile_path, "df_null_all.rds"), compress = "gz")
  write_rds(df_FS_null, paste0(readfile_path, "df_FS_null.rds"), compress = 'gz')
  write_rds(df_MD_null, paste0(readfile_path, "df_MD_null.rds"), compress = 'gz')
  write_rds(df_interval_null, paste0(readfile_path, "df_interval_null.rds"), compress = "gz")
  write_rds(df_FS_tmy_null, paste0(readfile_path, "df_FS_tmy_null.rds"), compress = "gz")
  
}

