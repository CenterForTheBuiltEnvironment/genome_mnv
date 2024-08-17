# Genome 2 M&V analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, # general
               broom, ggpmisc, ggpubr, # linear models
               RMV2.0, # TOWT model
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
function_path <- "./functions/"
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

# define parameters
# section run control
run_params <- list(type = "tidy", 
                   site = F, 
                   sprt = T, 
                   sprt_cont = T, 
                   nre_occ = T)

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
                   resamp = c(2, 8), 
                   cont_weeks = 36)

# NRE: Occupancy change
occ_params <- list(change_start = c(1, 5, 9),
                   change_end = c(4, 8, 12),
                   change = c(10, 20))





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
  
  return(list(name = name,
              site = site, 
              sprt = sprt_res %>% filter(flag == 1) %>% slice(1) %>% .$n_weeks,
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
readfile_path <- str_glue("./readfiles/{run_params$type}/")
summaryfigs_path <- str_glue("./figs/{run_params$type}/site_summary/")
combifigs_path <- str_glue("./figs/{run_params$type}/comb_analysis/")

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds"))
df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))
df_weather <- read_rds(paste0(readfile_path, "df_weather.rds"))

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





#### SITE ####
if (run_params$site){
  colors <- brewer.pal(nrow(all_types), "Set3")
  type_colors <- setNames(as.list(colors), all_types$type)
  df_tmy <- get_tmy(all_sites$site)
  
  
  # site and type number summary
  p1 <- df_energy %>%
    group_by(type) %>%
    distinct(name) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(type = fct_reorder(type, n, .desc = F)) %>%
    ggplot(aes(x = 1, y = n, fill = as.factor(type))) +
    geom_col() +
    geom_text(aes(label = n), color = "black", position = position_stack(vjust = 0.5)) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4)) +
    scale_x_discrete(expand = c(0, 0.1)) +
    scale_fill_manual(values = type_colors) +
    labs(x = NULL,
         y = NULL,
         subtitle = "Across all sites",
         fill = NULL) +
    theme(axis.text.x = element_blank(),
          panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  p2 <- df_energy %>%
    group_by(site, type) %>%
    distinct(name) %>%
    summarise(n = n()) %>%
    mutate(proportion = n / sum(n),
           total = sum(n),
           ymax = cumsum(proportion),
           ymin = c(0, head(ymax, n = -1))) %>%
    mutate(label_pos = (ymax + ymin) / 2) %>%
    ungroup() %>%
    group_by(type) %>%
    mutate(order = sum(n)) %>%
    ungroup() %>%
    mutate(type = fct_reorder(type, order, .desc = F)) %>%
    ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = type)) +
    geom_rect() +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    facet_wrap(~ site, nrow = 2) +
    labs(x = NULL,
         y = NULL,
         subtitle = "For each site",
         fill = NULL) +
    scale_fill_manual(values = type_colors) +
    geom_text(aes(x = 3.5, y = label_pos, label = n)) +
    geom_text(aes(x = 2, y = 0, label = paste0("Total\n", total)), color = "black") +
    theme(legend.direction = "horizontal",
          axis.text = element_blank(),
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggarrange(p1, p2,
            ncol = 2, nrow = 1,
            widths = c(0.25, 1),
            common.legend = TRUE,
            legend = "bottom") +
    plot_annotation(title = "Case study building type summary")
  
  ggsave(filename = "site_summary.png", path = summaryfigs_path, units = "in", height = 6, width = 8, dpi = 300)
  
  # EUI summary across all buildings
  df_eui <- df_energy %>%
    group_by(name) %>%
    summarise(total_eload = sum(eload, na.rm = T),
              type = unique(type),
              site = unique(site)) %>%
    ungroup() %>%
    left_join(df_meta, by = c("site", "type", "name")) %>%
    mutate(eui = total_eload / sqm) %>%
    group_by(site, type, name) %>%
    summarise(avg_eui = eui / 2) %>%
    ungroup() %>%
    mutate(name = fct_reorder(name, type, .desc = F))
  
  
  bar_plot <- df_eui %>%
    ggplot() +
    geom_col(aes(x = name, y = avg_eui, fill = type)) +
    scale_fill_manual(values = type_colors) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kWh/m\U00B2")) +
    labs(x = "All buildings",
         y = NULL,
         fill = NULL) +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          axis.text.x = element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 0.2, r = -0.5, unit = "cm"))
  
  hist_plot <- df_eui %>%
    ggplot(aes(x = avg_eui)) +
    geom_histogram(stat = "bin", position = "identity") +
    labs(title = NULL,
         x = NULL,
         y = "Frequency") +
    coord_flip() +
    theme(panel.grid.major.x = element_line(color = "grey80"),
          axis.text.y = element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 0.2, l = -0.5, unit = "cm"))
  
  combined_plot <- plot_grid(
    bar_plot,
    hist_plot,
    align = 'hv',
    nrow = 1,
    rel_widths = c(2, 1)
  )
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Electricity usage intensity for all buildings", fontface = 'bold', size = 18),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  print(final_plot)
  
  ggsave(filename = "eui_all.png", path = summaryfigs_path, units = "in", height = 5, width = 8, dpi = 300)
  
  df_energy %>%
    mutate(year = as.factor(year(timestamp))) %>%
    group_by(name, year) %>%
    summarise(total_eload = sum(eload, na.rm = T),
              type = unique(type),
              site = unique(site)) %>%
    ungroup() %>%
    left_join(df_meta, by = c("site", "type", "name")) %>%
    mutate(eui = total_eload / sqm) %>%
    ggplot() +
    geom_boxplot(aes(x = type, y = eui, fill = year), outlier.shape = NA) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, 450),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kWh/m\U00B2")) +
    labs(x = NULL,
         y = NULL,
         fill = NULL,
         title = "Electricity usage intensity for each building type") +
    coord_flip() +
    theme(panel.grid.major.x = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = "eui_type.png", path = summaryfigs_path, units = "in", height = 5, width = 6, dpi = 300)
  
  # Weather condition summary across all sites
  df_weather %>%
    mutate(year = as.factor(year(timestamp)),
           datetime = as.POSIXct(format(timestamp, "%m-%d %H:%M:%S"), format="%m-%d %H:%M:%OS")) %>%
    ggplot(aes(x = datetime, y = t_out, color = year)) +
    geom_point(size = 0.1, alpha = 0.05) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 1) +
    scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " °C")) +
    facet_wrap(~site, nrow = 2) +
    coord_cartesian(ylim = c(-8, 38)) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         title = "Outdoor weather conditions across all sites") +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = "weather_summary.png", path = summaryfigs_path, units = "in", height = 5, width = 10, dpi = 300)
  
  # density plot of tmy and actual tout
  df_weather %>%
    filter(year(timestamp) == 2016) %>%
    mutate(type = "observed") %>%
    ggplot() +
    geom_histogram(aes(x = t_out, fill = type), alpha = 0.6, color = NA) +
    geom_histogram(data = df_tmy %>% mutate(type = "tmy"), aes(x = temp, fill = type), alpha = 0.6, color = NA) +
    scale_fill_manual(values = c("observed" = "#5941A9", "tmy" = "#E5D4ED")) +
    facet_wrap(~site, nrow = 2) +
    labs(x = NULL, y = NULL, color = NULL, fill = NULL,
         title = "Outdoor temperature histogram\nfor observed and typical conditions") +
    scale_x_continuous(labels = number_format(suffix = " °C")) +
    coord_cartesian(clip = "off") +
    theme(axis.text.y = element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom")
  
  ggsave(filename = "weather_tmy.png", path = summaryfigs_path, units = "in", height = 5, width = 10, dpi = 300)
  
  # Energy and temperature dependence across each site
  df_all <- df_energy %>%
    left_join(df_weather, by = c("timestamp", "site"))
  
  df_daily <- df_all %>%
    group_by(site, type, date = floor_date(timestamp, unit = "day")) %>%
    summarise(t_out = mean(t_out, na.rm = T),
              eload = mean(eload, na.rm = T)) %>%
    ungroup()
  
  plot_list <- list()
  
  for (i in all_sites$site) {
    
    df <- df_daily %>%
      filter(site == i)
    
    df$type <- factor(df$type, levels = all_types$type)
    
    p <- df %>%
      ggplot(aes(x = t_out, y = eload, color = type)) +
      geom_point(data = .%>%
                   group_by(type) %>%
                   slice_sample(n = 100),
                 size = 0.5,
                 alpha = 0.5,
                 shape = 16) +
      geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.5, alpha = 0.1) +
      scale_x_continuous(expand = c(0, 0),
                         breaks = breaks_pretty(n = 3),
                         labels = number_format(suffix = " °C")) +
      scale_y_continuous(expand = c(0, 0),
                         breaks = breaks_pretty(n = 4),
                         labels = number_format(suffix = " kW")) +
      scale_color_manual(values = type_colors) +
      labs(x = NULL,
           y = NULL,
           color = NULL,
           subtitle = str_glue("{i}")) +
      theme(panel.grid.major.y = element_line(),
            legend.direction = "horizontal",
            legend.position = "none")
    
    # Store the plot in the list
    plot_list[[i]] <- p
  }
  
  temp_legend <- ggplot(df_daily, aes(x=t_out, y = eload, color = type)) +
    geom_point(size = 0.5, alpha = 0.2, shape = 16) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.5, alpha = 0.1) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " °C")) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = type_colors) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         subtitle = str_glue("{i}")) +
    theme(panel.grid.major.y = element_line(),
          legend.direction = "horizontal",
          legend.position = "bottom")
  
  legend <- get_legend(temp_legend)
  
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 2)
  combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.2))
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Daily average power usage by average outdoor weather conditions", fontface = 'bold', size = 18),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  print(final_plot)
  
  ggsave(filename = "energy_weather.png", path = summaryfigs_path, units = "in", height = 5, width = 10, dpi = 300)
  
  # energy consumption by each hour of the day
  df_hour <- df_energy %>%
    group_by(site,
             type,
             timestamp) %>%
    summarize_if(is.numeric, mean, na.rm=T) %>%
    ungroup()
  
  plot_list <- list()
  
  for (i in all_sites$site) {
    
    df <- df_hour %>%
      filter(site == i)
    
    p <- df %>%
      ggplot(aes(x = hour(timestamp), y = eload, color = type)) +
      geom_point(data = .%>%
                   group_by(type) %>%
                   slice_sample(n = 100),
                 size = 0.5,
                 alpha = 0.5,
                 shape = 16) +
      geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.5, alpha = 0.1) +
      scale_x_continuous(breaks = c(0, 6, 12, 18),
                         labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
      scale_y_continuous(expand = c(0, 0),
                         breaks = breaks_pretty(n = 4),
                         labels = number_format(suffix = " kW")) +
      scale_color_manual(values = type_colors) +
      labs(x = NULL,
           y = NULL,
           color = NULL,
           subtitle = str_glue("{i}")) +
      theme(panel.grid.major.y = element_line(),
            legend.direction = "horizontal",
            legend.position = "none")
    
    # Store the plot in the list
    plot_list[[i]] <- p
  }
  
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 2)
  combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.2))
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Average power usage by each hour of the day", fontface = 'bold', size = 18),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  print(final_plot)
  
  ggsave(filename = "energy_hour_site.png", path = summaryfigs_path, units = "in", height = 5, width = 10, dpi = 300)
}





#### INDIVIDUAL ####

# storing results
FS_ref <- list()
FS_occ <- list()
MD_ref <- list()
model_acc <- list()
seq_timeline <- list()
seq_mdsaving <- list()
seq_nmsaving <- list()
seq_frsaving <- list()
cont_mdsaving <- list()
cont_frsaving <- list()
energy <- list()

# get site information
# for (n in 1:2){
for (n in 1:(nrow(all_names))){
  
  name <- all_names$name[n]
  
  site_info <- df_energy %>%
    filter(name == all_names$name[n]) %>%
    select(site, type) %>%
    distinct()
  
  site <- site_info$site
  
  ifelse(!dir.exists(file.path(str_glue("./figs/{run_params$type}/site_analysis/{site}/{name}"))), dir.create(file.path(str_glue("./figs/{run_params$type}/site_analysis/{site}/{name}"))), FALSE)
  sitefigs_path <- str_glue("./figs/{run_params$type}/site_analysis/{site}/{name}")
  
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
  
  # power-temp plot
  p1 <- df_hourly_conv %>%
    pivot_longer(c(base_eload, interv_eload), names_to = "strategy", values_to = "eload") %>%
    mutate(strategy = as.factor(strategy),
           strategy = recode_factor(strategy, "base_eload" = "Baseline", "interv_eload" = "Intervention")) %>%
    ggplot(aes(x = t_out, y = eload, color = strategy)) +
    geom_point(data= .%>%
                 group_by(strategy) %>%
                 slice_sample(n = 1000),
               size = 0.7,
               alpha = 0.4,
               shape = 16) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 1.25, alpha = 0.15) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " °C")) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = plot_scale) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         subtitle = "by outdoor drybulb temperature") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  # hour-of-day visualization
  p2 <- df_hourly_conv %>%
    pivot_longer(c(base_eload, interv_eload), names_to = "strategy", values_to = "eload") %>%
    mutate(strategy = as.factor(strategy),
           strategy = recode_factor(strategy, "base_eload" = "Baseline", "interv_eload" = "Intervention")) %>%
    ggplot(aes(x = hour(datetime), y = eload, color = strategy)) +
    geom_point(data= .%>%
                 group_by(strategy) %>%
                 slice_sample(n = 1000),
               size = 0.7,
               alpha = 0.4,
               shape = 16) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 1.25, alpha = 0.15) +
    scale_x_continuous(breaks = c(0, 6, 12, 18),
                       labels = c("12 AM", "6 AM", "12 PM", "6 PM")) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = plot_scale) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         subtitle = "by each hour of the day") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.y = element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggarrange(p1, p2,
            ncol=2, nrow=1,
            labels = c("a)", "b)"),
            widths = c(1, 1),
            align = "v",
            common.legend = TRUE,
            legend="bottom") +
    plot_annotation(title = "Case study building power consumption (hourly average)")
  
  ggsave(filename = "temp_time_power.png", path = sitefigs_path, units = "in", height = 5, width = 10, dpi = 300)
  
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
  
  base_proj %>%
    ggplot() +
    geom_point(aes(x = datetime, y = eload, color = "Measurement"), alpha = 0.2, size = 0.2) +
    geom_point(aes(x = datetime, y = towt, color = "Prediction"), alpha = 0.2, size = 0.2) +
    geom_smooth(aes(x = datetime, y = eload, color = "Measurement"), formula = y ~ x, method = "loess", linewidth = 0.7) +
    geom_smooth(aes(x = datetime, y = towt, color = "Prediction"), formula = y ~ x, method = "loess", linewidth = 0.7) +
    annotate(geom = "text",
             x = median(base_proj$datetime),
             y = median(base_proj$eload) + 2 * sd(base_proj$eload),
             label = paste0("CV(RMSE): ", round(cv_rmse, digits = 2), "%")) +
    scale_color_brewer(palette = "Set1") +
    scale_x_datetime(date_breaks = "2 months",
                     date_labels = "%b")  +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " kW")) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         title = "TOWT model prediction results") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = "towt_acc.png", path = sitefigs_path, units = "in", height = 6, width = 10, dpi = 300)
  
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
                       seed = sample(1, 2^15, 1),
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
  
  p1 <- df_hourly_conv %>%
    mutate(savings = base_eload - interv_eload,
           year = as.factor(year(datetime))) %>%
    ggplot(aes(x = datetime, y = savings)) +
    geom_point(alpha = 0.2, size = 0.2) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.7) +
    facet_wrap(~year, scales = "free_x") +
    scale_x_datetime(date_breaks = "2 months",
                     date_labels = "%b")  +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " kW")) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         subtitle = "Calculated savings") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.direction = "horizontal",
          axis.text.x = element_blank(),
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  p2 <- df_hourly_conv %>%
    mutate(year = as.factor(year(datetime))) %>%
    ggplot(aes(x = datetime, y = t_out)) +
    geom_point(alpha = 0.2, size = 0.2) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = 0.7) +
    facet_wrap(~year, scales = "free_x") +
    scale_x_datetime(date_breaks = "2 months",
                     date_labels = "%b")  +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " °C")) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         subtitle = "Measured outdoor weather") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggarrange(p1, p2,
            ncol=1, nrow=2,
            labels = c("a)", "b)"),
            widths = c(1, 1),
            align = "h",
            common.legend = TRUE,
            legend="bottom") +
    plot_annotation(title = "Case study building intervention savings")
  
  ggsave(filename = "savings_true.png", path = sitefigs_path, units = "in", height = 8, width = 10, dpi = 300)
  
  
  
  
  
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
    
    for (i in 2:sprt_param$n_weeks){
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
    
    cont_param$n_weeks <- cont_param$cont_weeks + eob
    start_datetime <- as.Date("2016-01-01") + weeks(eob)
    end_datetime <- start_datetime + weeks(cont_param$cont_weeks)  
    
    date_seq <- seq(from = start_datetime, by = "day", length.out = cont_param$cont_weeks * 7)
    strategy <- sample(c(1, 2), size = length(date_seq), replace = TRUE, prob = cont_param$resamp)
    
    df_schedule_cont <- data.frame(datetime = date_seq,
                                   strategy = strategy)
    
    df_rand_old <- df_rand %>% 
      filter(datetime < start_datetime)
    
    df_rand_new <- df_rand %>% 
      filter(datetime >= start_datetime & datetime < end_datetime) %>% 
      select(-strategy) %>% 
      left_join(df_schedule_cont, by = "datetime") %>% 
      fill(strategy, .direction = "down")
    
    df_rand_new <- df_hourly_conv %>%
      filter(datetime >= start_datetime & datetime < end_datetime) %>% 
      left_join(df_schedule_cont, by = "datetime") %>%
      fill(strategy, .direction = "down") %>%
      pivot_longer(c(base_eload, interv_eload), names_to = "eload_type", values_to = "eload") %>%
      filter((strategy == 1 & eload_type == "base_eload") | (strategy == 2 & eload_type == "interv_eload")) %>%
      select(-eload_type) %>% 
      drop_na()
    
    df_rand_cont <- rbind(df_rand_old, df_rand_new)
    
    seq_res <- seq_run(cont_param, df_rand_cont, site_tmy)
    annual_saving <- seq_res$annual_saving
    df_means <- seq_res$df_means
    sprt_res <- seq_res$sprt_res
    
    # get true savings
    df_week <- df_base_conv %>% 
      select(datetime, 
             base_eload = eload) %>% 
      left_join(df_interv_conv, by = "datetime") %>% 
      mutate(savings = eload - base_eload) %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week <= cont_param$n_weeks)
    
    true_saving <- list()
    
    for (i in 2:cont_param$n_weeks){
      saving <- df_week %>% 
        filter(week <= i)
      
      true_saving[[i]] <- tibble("n_weeks" = i, 
                                 savings = mean(saving %>% .$savings))
    }
    
    true_saving <- bind_rows(true_saving)
    
    # plot overall results
    cont_plot(df_means, sprt_res, annual_saving, true_saving, eob)
    ggsave(filename = "cont_seq.png", path = sitefigs_path, units = "in", height = 8, width = 8, dpi = 300)
    
    # savings at timeline
    cont_mdsaving[[n]] <- list("name" = name, 
                               "site" = site, 
                               "cont" = sprt_res %>% filter(n_weeks == cont_param$n_weeks) %>% .$ns_stat)
    
    # saving calculation as mean difference
    fs_rand <- (mean(df_rand_cont %>% filter(strategy == 1) %>% .$eload) -
                  mean(df_rand_cont %>% filter(strategy == 2) %>% .$eload)) /
      mean(df_rand_cont %>% filter(strategy == 1) %>% .$eload) * 100
    
    cont_frsaving[[n]] <- list("name" = name, 
                               "site" = site, 
                               "cont" = fs_rand)
    
    energy[[n]] <- list("name" = name, 
                        "site" = site, 
                        "ref" = sum(df_base_conv %>% filter(datetime < as.Date("2017-01-01")) %>% .$eload), 
                        "rand" = sum(df_rand_cont %>% filter(datetime < as.Date("2017-01-01")) %>% .$eload), 
                        "interv" = sum(df_interv_conv %>% filter(datetime < as.Date("2017-01-01")) %>% .$eload))
    
  }
  
  
  
  
  #### NRE ####
  if (run_params$nre_occ) {
    
    # Occupancy change
    scenario <- 1
    tibble_occ <- list()
    tibble_occ[[n]] <- tibble("name" = name, 
                              "site" = site)
    
    for (time in 1:length(occ_params$change_start)){
      for (change in occ_params$change){
        
        change_start <- occ_params$change_start[time]
        change_end <- occ_params$change_end[time]
        
        base_pre_meas <- df_base_conv %>%
          filter(datetime < as.Date("2017-01-01")) %>%
          mutate(eload = ifelse(month(datetime) >= month(change_start) & month(datetime) <= month(change_end), eload * (1 - change / 100), eload))
        
        interv_pre_true <- df_interv_conv %>%
          filter(datetime < as.Date("2017-01-01")) %>%
          mutate(eload = ifelse(month(datetime) >= month(change_start) & month(datetime) <= month(change_end), eload * (1 - change / 100), eload))
        
        # Baseline projection
        towt_base <- base_pre_meas %>%
          model_fit()
        
        df_towt <- df_interv_conv %>%
          filter(datetime >= as.Date("2017-01-01")) %>%
          select(time = datetime,
                 temp = t_out,
                 eload)
        
        base_pos_proj <- model_pred(df_towt, towt_base) %>%
          rename("datetime" = "time") %>%
          select(datetime, towt, eload)
        
        base_pos_true <- df_base_conv %>%
          filter(datetime >= as.Date("2017-01-01"))
        
        interv_pos_meas <- df_interv_conv %>%
          filter(datetime >= as.Date("2017-01-01"))
        
        rand <- df_rand %>%
          filter(datetime < as.Date("2017-01-01")) %>%
          mutate(eload = ifelse(month(datetime) >= month(change_start) & month(datetime) <= month(change_end), eload * (1 - change / 100), eload)) %>%
          rbind(df_rand %>% filter(datetime >= as.Date("2017-01-01")))
        
        prepost_plot(base_pre_meas,  
                     base_pos_proj, 
                     base_pos_true, 
                     interv_pos_meas, 
                     str_glue("Building energy consumption over a 2-year period\n(2016-{change_start}-1 ~ 2016-{change_end}-30: {change}% change)"))
        
        ggsave(filename = str_glue("occ_change_conv_S{scenario}.png"), path = sitefigs_path, units = "in", height = 5, width = 10, dpi = 300)
        
        # saving calculation as mean difference
        base_true <- rbind(base_pre_meas, base_pos_true)
        interv_true <- rbind(interv_pre_true, interv_pos_meas)
        
        FS_true <- (mean(base_true$eload) - mean(interv_true$eload)) / mean(base_true$eload) * 100
        FS_conv <- (mean(base_pos_proj$towt) - mean(interv_pos_meas$eload)) / mean(base_pos_proj$towt) * 100
        FS_rand <- (mean(rand %>% filter(strategy == 1) %>% .$eload) -
                      mean(rand %>% filter(strategy == 2) %>% .$eload)) /
          mean(rand %>% filter(strategy == 1) %>% .$eload) * 100
        
        new_tibble <- as_tibble(setNames(list(FS_true, FS_conv, FS_rand),
                                         c(str_glue("S{scenario}_true"),
                                           str_glue("S{scenario}_conv"),
                                           str_glue("S{scenario}_rand"))))
        
        tibble_occ[[n]] <- bind_cols(tibble_occ[[n]], new_tibble)
        
        scenario <- scenario + 1
        
      }
      
    }
    
    FS_occ <- bind_rows(FS_occ, tibble_occ)
    
  }
  
  print(paste0("Finished: ", n, "/", nrow(all_names)))
  
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

df_NRE_occ <- df_FS %>%
  rbind(FS_occ %>% 
          pivot_longer(-c(name, site), names_to = "type", values_to = "savings") %>%
          separate(type, into = c("scenario", "method"), sep = "_"))

df_eui <- bind_rows(energy) %>% 
  left_join(df_meta, by = c("name", "site")) 

# store savings and timeline
write_rds(df_sprt_all, paste0(readfile_path, "df_sprt_all.rds"), compress = "gz")
write_rds(df_seq_FS, paste0(readfile_path, "df_seq_FS.rds"), compress = "gz")
write_rds(df_NRE_occ, paste0(readfile_path, "df_NRE_occ.rds"), compress = "gz")
write_rds(df_MD, paste0(readfile_path, "df_MD.rds"), compress = "gz")
write_rds(df_FS, paste0(readfile_path, "df_FS.rds"), compress = "gz")
write_rds(df_eui, paste0(readfile_path, "df_eui.rds"), compress = "gz")

# Normalized savings
p1 <- df_sprt_all %>% 
  filter(seq == "eob") %>% 
  left_join(df_NRE_occ %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  left_join(df_seq_FS %>% filter(seq == "eob"), by = c("name", "site", "seq")) %>% 
  mutate(Deviation = abs(savings - FS),
         n_weeks = as.factor(n_weeks)) %>% 
  rename("Normalized on TMY" = annual, 
         "Measured" = FS) %>% 
  pivot_longer(c("Measured", "Normalized on TMY", Deviation), names_to = "parameter", values_to = "value") %>% 
  ggplot(aes(x = n_weeks, y = value, fill = parameter)) +
  stat_boxplot(geom ='errorbar', width = 0.75) +
  geom_boxplot(width = 0.5, position = position_dodge(0.75), outlier.shape = NA) +
  scale_x_discrete(expand = c(0.5, 0),
                   labels = c("24\nweeks", "36\nweeks")) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     limits = c(0, 16), 
                     labels = number_format(suffix = "%")) +
  labs(x = NULL, 
       y = NULL, 
       fill = NULL, 
       subtitle = "Early stop after satisfying all criteria") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- df_sprt_all %>% 
  filter(seq == "final") %>% 
  left_join(df_NRE_occ %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  left_join(df_seq_FS %>% filter(seq == "final"), by = c("name", "site", "seq")) %>% 
  mutate(Deviation = abs(savings - FS),
         n_weeks = as.factor(n_weeks)) %>% 
  rename("Normalized on TMY" = annual, 
         "Measured" = FS) %>% 
  pivot_longer(c("Measured", "Normalized on TMY", Deviation), names_to = "parameter", values_to = "value") %>% 
  ggplot(aes(x = n_weeks, y = value, fill = parameter)) +
  stat_boxplot(geom ='errorbar', width = 0.75) +
  geom_boxplot(width = 0.5, position = position_dodge(0.75), outlier.shape = NA) +
  scale_x_discrete(expand = c(0.5, 0),
                   labels = c("48\nweeks")) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     limits = c(0, 16), 
                     labels = number_format(suffix = "%")) +
  labs(x = NULL, 
       y = NULL, 
       fill = NULL, 
       subtitle = "At the end of the full testing period") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text.y = element_blank(), 
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(2, 1),
          align = "hv",
          common.legend = T, 
          legend="bottom") +
  plot_annotation(title = "Fractional savings estimation using randomized M&V")

ggsave(filename = str_glue("FS_Dev.png"), path = combifigs_path, units = "in", height = 5, width = 8, dpi = 300)

# Sequential mean difference 
plot_data <- df_sprt_all %>% 
  filter(seq == "eob") %>% 
  left_join(df_MD %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - - sprt), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(x = abs_diff), width = 0.1) +
  geom_text(aes(x = median(abs_diff), y = 0, label = paste0(round(median(abs_diff), digits = 1), " kW")), 
            check_overlap = TRUE, 
            position = position_nudge(y = -0.1)) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = " kW")) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- plot_data %>% 
  mutate(bin = cut(abs_diff, breaks = 20, include.lowest = T)) %>%
  group_by(bin) %>%
  summarise(count = n(), 
            knot = max(abs_diff)) %>%
  arrange(knot) %>%
  ungroup() %>% 
  mutate(cumulative_count = cumsum(count)) %>% 
  ggplot() +
  geom_line(aes(x = knot, y = cumulative_count), alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(x = knot, y = cumulative_count)) +
  geom_text(aes(x = knot, 
                y = cumulative_count, 
                label = paste0("(", round(knot, digits = 1), " kW, ", cumulative_count, ")")), 
            position = position_nudge(y = 1), 
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = " kW")) +
  scale_y_continuous(limits = c(50, nrow(all_names) + 1)) +
  labs(x = "Absolute difference in measured savings", 
       y = "Number of buildings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol=1, nrow=2,
          labels = c("a)", "b)"),
          heights = c(1, 4),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in saving estimation at early stop\n(Measured savings from observed weather)"))

ggsave(filename = str_glue("md_comp.png"), path = combifigs_path, units = "in", height = 8, width = 8, dpi = 300)

# continuous sprt mean difference 
plot_data <- bind_rows(cont_mdsaving) %>% 
  left_join(df_MD %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - - cont), 
         plot_max = max(abs_diff))

p1 <-  plot_data %>% 
  ggplot() +
  geom_boxplot(aes(x = abs_diff), width = 0.1) +
  geom_text(aes(x = median(abs_diff), y = 0, label = paste0(round(median(abs_diff), digits = 1), " kW")), 
            check_overlap = TRUE, 
            position = position_nudge(y = -0.1)) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = " kW")) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- plot_data %>% 
  mutate(bin = cut(abs_diff, breaks = 20, include.lowest = T)) %>%
  group_by(bin) %>%
  summarise(count = n(), 
            knot = max(abs_diff)) %>%
  arrange(knot) %>%
  ungroup() %>% 
  mutate(cumulative_count = cumsum(count)) %>% 
  ggplot() +
  geom_line(aes(x = knot, y = cumulative_count), alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(x = knot, y = cumulative_count)) +
  geom_text(aes(x = knot, 
                y = cumulative_count, 
                label = paste0("(", round(knot, digits = 1), " kW, ", cumulative_count, ")")), 
            position = position_nudge(y = 1), 
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = " kW")) +
  scale_y_continuous(limits = c(50, nrow(all_names) + 1)) +
  labs(x = "Absolute difference in measured savings", 
       y = "Number of buildings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol=1, nrow=2,
          labels = c("a)", "b)"),
          heights = c(1, 4),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in saving estimation after sequential test\n(Measured savings from observed weather)"))

ggsave(filename = str_glue("md_comp_cont.png"), path = combifigs_path, units = "in", height = 8, width = 8, dpi = 300)

# Sequential fractional savings
plot_data <- df_seq_FS %>% 
  filter(seq == "eob") %>% 
  left_join(df_NRE_occ %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - FS), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(x = abs_diff), width = 0.1) +
  geom_text(aes(x = median(abs_diff), y = 0, label = paste0(round(median(abs_diff), digits = 1), "%")), 
            check_overlap = TRUE, 
            position = position_nudge(y = -0.1)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = "%")) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- plot_data %>% 
  mutate(bin = cut(abs_diff, breaks = 20, include.lowest = T)) %>%
  group_by(bin) %>%
  summarise(count = n(), 
            knot = max(abs_diff)) %>%
  arrange(knot) %>%
  ungroup() %>% 
  mutate(cumulative_count = cumsum(count)) %>% 
  ggplot() +
  geom_line(aes(x = knot, y = cumulative_count), alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(x = knot, y = cumulative_count)) +
  geom_text(aes(x = knot, 
                y = cumulative_count, 
                label = paste0("(", round(knot, digits = 1), "%, ", cumulative_count, ")")), 
            position = position_nudge(y = 1), 
            size = 4) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = "%")) +
  scale_y_continuous(limits = c(40, nrow(all_names) + 1)) +
  labs(x = "Absolute difference in normalized fractional savings", 
       y = "Number of buildings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol=1, nrow=2,
          labels = c("a)", "b)"),
          heights = c(1, 4),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in fractional saving estimation at early stop\n(Measured savings from observed weather)"))

ggsave(filename = str_glue("fr_comp.png"), path = combifigs_path, units = "in", height = 8, width = 8, dpi = 300)

# continuous sprt fractional savings
plot_data <- bind_rows(cont_frsaving) %>% 
  left_join(df_FS %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - cont), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(x = abs_diff), width = 0.1) +
  geom_text(aes(x = median(abs_diff), y = 0, label = paste0(round(median(abs_diff), digits = 1), "%")), 
            check_overlap = TRUE, 
            position = position_nudge(y = -0.1)) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = "%")) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p2 <- plot_data %>% 
  mutate(bin = cut(abs_diff, breaks = 20, include.lowest = T)) %>%
  group_by(bin) %>%
  summarise(count = n(), 
            knot = max(abs_diff)) %>%
  arrange(knot) %>%
  ungroup() %>% 
  mutate(cumulative_count = cumsum(count)) %>% 
  ggplot() +
  geom_line(aes(x = knot, y = cumulative_count), alpha = 0.2, linewidth = 1.2) +
  geom_point(aes(x = knot, y = cumulative_count)) +
  geom_text(aes(x = knot, 
                y = cumulative_count, 
                label = paste0("(", round(knot, digits = 1), "%, ", cumulative_count, ")")), 
            position = position_nudge(y = 1), 
            size = 4) +
  scale_x_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3),
                     labels = number_format(suffix = "%")) +
  scale_y_continuous(limits = c(40, nrow(all_names) + 1)) +
  labs(x = "Absolute difference in normalized fractional savings", 
       y = "Number of buildings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(xlim = c(0, max(plot_data$plot_max) + 1)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p1, p2,
          ncol=1, nrow=2,
          labels = c("a)", "b)"),
          heights = c(1, 4),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in fractional saving estimation after sequential test\n(Measured savings from observed weather)"))

ggsave(filename = str_glue("fr_comp_cont.png"), path = combifigs_path, units = "in", height = 8, width = 8, dpi = 300)


# Timeline
df_sprt_all %>% 
  filter(seq != "final") %>% 
  group_by(name) %>% 
  arrange(name, n_weeks) %>% 
  mutate(increment = n_weeks - lag(n_weeks)) %>% 
  mutate(start = min(n_weeks)) %>% 
  mutate(across(everything(), ~ ifelse(is.na(.), start, .))) %>% 
  mutate(seq = as.factor(seq), 
         seq = recode_factor(seq, "sprt" = "SPRT check", "temp" = "Temperature check", "eob" = "Block end check")) %>% 
  ungroup() %>% 
  ggplot(aes(group = site)) +
  geom_col(aes(x = name, y = increment, fill = seq), position = "stack") +
  facet_wrap(~site, ncol = 1, scales = "free_y") +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(df_sprt_all %>% filter(seq != "final") %>% .$n_weeks) + 0.5),
                     breaks = seq(0, sprt_param$n_weeks, by = block_params$block_unit), 
                     labels = number_format(suffix = "\nweeks")) +
  geom_hline(yintercept = seq(block_params$block_unit, max(df_sprt_all %>% filter(seq != "final") %>% .$n_weeks), by = block_params$block_unit), 
             lty = "dashed", 
             color = "grey20") +
  scale_fill_brewer() +
  coord_flip() +
  labs(x = NULL,
       fill = NULL,
       y = NULL,
       title = "Overall sequential test complete timeline") +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("seq_timeline.png"), path = combifigs_path, units = "in", height = 18, width = 10, dpi = 300)

# continue sampling visualization
df_rand_cont %>% 
  mutate(strategy = as.factor(strategy), 
         strategy = recode_factor(strategy, "1" = "Baseline", "2" = "Intervention")) %>% 
  mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor(), 
         sample = as.factor(ifelse(week <= eob, "Before", "After")), 
         sample = fct_relevel(sample, c("Before", "After"))) %>%
  group_by(sample, strategy) %>% 
  summarise(n = n()) %>% 
  mutate(n = n, 
         perc = n / sum(n)) %>%
  ungroup() %>% 
  ggplot(aes(x = sample, y = perc, fill = strategy)) +
  geom_bar(position="fill", stat="identity") +
  geom_hline(yintercept = 0.5, 
             linetype = "dashed", 
             color = "red") +
  annotate(geom = "text",
           color = "red", 
           x = 0.5, 
           y = 0.52, 
           label = "50%") +
  labs(x = NULL,
       fill = NULL,
       y = NULL,
       title = "Sampling ratio between baseline and intervention", 
       subtitle = "Before and after the sequential test") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = "%", scale = 100)) +
  scale_fill_manual(values = ls_colors) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# Continue sampling energy saving
mean_perc <- df_eui %>% 
  mutate(ref_saving = (ref - interv) / sqm, 
         rand_saving = (ref - rand) / sqm, 
         perc = rand_saving / ref_saving * 100) %>% 
  group_by(site) %>%
  summarise(perc_mean = round(mean(perc), digits = 0),
            pos = round(n() / 2) + 1)

df_eui %>% 
  mutate(ref_saving = (ref - interv) / sqm, 
         rand_saving = (ref - rand) / sqm) %>% 
  ggplot() +
  geom_col(aes(x = name, y = ref_saving, fill = "Estimated savings"), position = "identity") +
  geom_col(aes(x = name, y = rand_saving, fill = "Measured savings"), position = "identity") +
  geom_text(data = mean_perc, aes(x = pos, y = -.5, group = site, label = paste0(perc_mean, "%"))) +
  geom_text(aes(x = pos, y = conv), 
            mean_perc %>% 
              filter(site == df_eui %>% tail(1) %>% .$site) %>% 
              mutate(conv = 20), 
            color = "grey40",
            label = "Conventional M&V\nbase-year savings: 0") +
  facet_wrap(~site, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = c(0.1, 0),
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " kWh/m^2")) +
  scale_fill_brewer() +
  labs(x = NULL,
       fill = NULL,
       y = NULL,
       color = NULL, 
       title = "First-year normalized savings on floor area through randomization", 
       subtitle = "baseline and intervention sampled at 50%/50% followed by 20%/80%") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("saving_cont.png"), path = combifigs_path, units = "in", height = 5, width = 15, dpi = 300)

# NRE saving estimation: occupancy change
for (s in c("ref", "S1", "S2", "S3", "S4", "S5", "S6")){
  dev_conv <- df_NRE_occ %>%
    filter(scenario == s) %>%
    filter(method != "rand") %>%
    group_by(name, site, scenario) %>%
    summarise(conv = abs(diff(savings))) %>%
    ungroup()
  
  dev_rand <- df_NRE_occ %>%
    filter(scenario == s) %>%
    filter(method != "conv") %>%
    group_by(name, site, scenario) %>%
    summarise(rand = abs(diff(savings))) %>%
    ungroup()
  
  dev_FS <- dev_rand %>%
    left_join(dev_conv, by = c("name", "site", "scenario")) %>% 
    mutate(diff_in_diff = conv - rand)
  
  mean_diff <- dev_FS %>%
    group_by(site) %>%
    summarise(mean = round(mean(diff_in_diff), digits = 1),
              pos = round(n() / 2) + 1,
              .groups = 'keep')
  
  dev_FS %>%
    ggplot(aes(group = site)) +
    geom_col(aes(x = name, y = diff_in_diff), position = "identity", alpha = 0.5) +
    facet_wrap(~site, nrow = 1, scales = "free_x") +
    scale_y_continuous(expand = c(0.1, 0),
                       breaks = breaks_pretty(n = 4), 
                       labels = number_format(suffix = "%")) +
    geom_text(data = mean_diff,
              aes(x = pos, y = -.5, group = site, label = paste0(mean, "%"))) +
    geom_line(aes(x = name, y = rand, color = "Absolute deviation of randomized method"), alpha = 0.5) +
    geom_point(aes(x = name, y = rand, color = "Absolute deviation of randomized method"), alpha = 0.5) +
    labs(x = NULL,
         fill = NULL,
         y = NULL,
         color = NULL, 
         title = "Difference-in-difference of fractional savings calculated for each building") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = str_glue("occ_{s}_savings.png"), path = combifigs_path, units = "in", height = 5, width = 15, dpi = 300)
}
