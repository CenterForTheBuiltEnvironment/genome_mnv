#### LIBRARIES ####
require(pacman)

# load packages using pacman
pacman::p_load(tidyverse, lubridate, here, stats, zoo, scales, lvplot, ggpubr, gridExtra, patchwork, RColorBrewer)

# turn off scientific notation
options(scipen = 999, digits = 15)

# set directory
here::i_am("manuscript.rmd")

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

# colors
ls_colors <- c("Baseline" = "#c6dbef",
               "Measured baseline" = "#c6dbef",
               "Adjusted baseline" = "#c6dbef",
               "Projected baseline" = "#2171b5",
               "Intervention" = "#fdbb84",
               "Measured interv" = "#fdbb84",
               "Conventional" = "grey70", 
               "Randomized" = "#99d8c9", 
               "Randomized\n(24 months)" = "#66c2a4",
               "Randomized\n(50/50)" = "#66c2a4",
               "Daily\nsampling" = "#66c2a4",
               "2-day\nsampling" = "#41ae76",
               "3-day\nsampling" = "#238b45",
               "7-day\nsampling" = "#006d2c",
               "Randomized\n(20/80)" = "#ccece6",
               "Buildings finishing randomized M&V" =  "#99d8c9",
               "Buildings satisfying 80% TMY range" = "black", 
               "Buildings satisfying SPRT" = "black")

# parameters
ctr_params <- list(peak_hours = 10:16,
                   chwl_perc = 0.25,
                   step_perc = 0.08,
                   conv_swt = 6,
                   weather_knots = c(15, 25),
                   swt_knots = c(12, 7),
                   coe_peak = 0.8,
                   coe_off = 1.2,
                   enable_temp = 8)

occ_params <- list(change_start = 5,
                   change_end = 8,
                   change = 20)





#### FUNCTIONS ####
# Function defined to read downloaded tmy files
get_tmy <- function(all_sites, readfile_path){
  
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
    mutate(date = date(timestamp)) %>%
    group_by(date) %>%
    summarize(na_hours = sum(is.na(eload)))
  
  # Filter out days with more than half of the hours having NAs
  valid_days <- na_counts %>%
    filter(na_hours <= 12) %>%
    pull(date)
  
  df_filtered <- df_all %>%
    filter(date(timestamp) %in% valid_days)
  
  df_filtered <- df_filtered %>%
    mutate(across(c(eload, t_out), ~ zoo::na.approx(., na.rm = FALSE))) %>% 
    rename(datetime = timestamp, 
           base_eload = eload)
  
  return(df_filtered)
}

# Function defined to adjust the plot scale
get_scale <- function(eload, range = 2){
  
  min_y <- mean(eload, na.rm = T) - range * sd(eload, na.rm = T)
  max_y <- mean(eload, na.rm = T) + range * sd(eload, na.rm = T)
  
  return(c(min_y, max_y))
}





#### READ DATA ####
df_loc <- read.csv("../readfiles/loc_map.csv")

# stable dataset
readfile_stable <- str_glue("../readfiles/stable/")
fig_path = "../figs/manuscript/"

df_energy_stable <- read_rds(paste0(readfile_stable, "df_energy.rds"))
df_meta_stable <- read_rds(paste0(readfile_stable, "df_meta.rds"))
df_weather_stable <- read_rds(paste0(readfile_stable, "df_weather.rds"))
df_sprt_all_stable <- read_rds(paste0(readfile_stable, "df_sprt_all.rds"))
df_seq_FS_stable <- read_rds(paste0(readfile_stable, "df_seq_FS.rds"))
df_MD_stable <- read_rds(paste0(readfile_stable, "df_MD.rds"))
df_FS_stable <- read_rds(paste0(readfile_stable, "df_FS.rds"))
df_cont_stable <- read_rds(paste0(readfile_stable, "df_cont.rds"))
df_FS_null_stable <- read_rds(paste0(readfile_stable, "df_FS_null.rds"))
df_MD_null_stable <- read_rds(paste0(readfile_stable, "df_MD_null.rds"))
df_seq_FS_null_stable <- read_rds(paste0(readfile_stable, "df_seq_FS_null.rds"))
df_model_acc_stable <- read_rds(paste0(readfile_stable, "df_model_acc.rds"))
df_FS_tmy_stable <- read_rds(paste0(readfile_stable, "df_FS_tmy.rds"))
df_FS_tmy_null_stable <- read_rds(paste0(readfile_stable, "df_FS_tmy_null.rds"))
df_interval_stable <- read_rds(paste0(readfile_stable, "df_interval.rds"))
df_interval_null_stable <- read_rds(paste0(readfile_stable, "df_interval_null.rds"))
df_seq_interval_FS_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_FS.rds"))
df_seq_interval_nm_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_nm.rds"))


all_sites_stable <- df_energy_stable %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types_stable <- df_energy_stable %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

all_names_stable <- df_energy_stable %>%
  select(name) %>%
  distinct(name)

# variable dataset
readfile_variable <- str_glue("../readfiles/variable/")

df_energy_variable <- read_rds(paste0(readfile_variable, "df_energy.rds"))
df_meta_variable <- read_rds(paste0(readfile_variable, "df_meta.rds"))
df_weather_variable <- read_rds(paste0(readfile_variable, "df_weather.rds"))
df_sprt_all_variable <- read_rds(paste0(readfile_variable, "df_sprt_all.rds"))
df_seq_FS_variable <- read_rds(paste0(readfile_variable, "df_seq_FS.rds"))
df_MD_variable <- read_rds(paste0(readfile_variable, "df_MD.rds"))
df_FS_variable <- read_rds(paste0(readfile_variable, "df_FS.rds"))
df_cont_variable <- read_rds(paste0(readfile_variable, "df_cont.rds"))
df_FS_null_variable <- read_rds(paste0(readfile_variable, "df_FS_null.rds"))
df_MD_null_variable <- read_rds(paste0(readfile_variable, "df_MD_null.rds"))
df_seq_FS_null_variable <- read_rds(paste0(readfile_variable, "df_seq_FS_null.rds"))
df_model_acc_variable <- read_rds(paste0(readfile_variable, "df_model_acc.rds"))

df_FS_tmy_variable <- read_rds(paste0(readfile_variable, "df_FS_tmy.rds"))
df_FS_tmy_null_variable <- read_rds(paste0(readfile_variable, "df_FS_tmy_null.rds"))
df_interval_variable <- read_rds(paste0(readfile_variable, "df_interval.rds"))
df_interval_null_variable <- read_rds(paste0(readfile_variable, "df_interval_null.rds"))
df_seq_interval_FS_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_FS.rds"))
df_seq_interval_nm_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_nm.rds"))

all_sites_variable <- df_energy_variable %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types_variable <- df_energy_variable %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

all_names_variable <- df_energy_variable %>%
  select(name) %>%
  distinct(name)

# read functions
function_path <- "../functions/"
source(paste0(function_path, "model_fit.R"))
source(paste0(function_path, "model_pred.R"))
source(paste0(function_path, "prepost_plot.R"))

S_building <- nrow(df_FS_tmy_stable)
V_building <- nrow(df_FS_tmy_variable)
A_building <- S_building + V_building




#### GA ####
rand_eob_S <- df_seq_FS_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_FS_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site"))%>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY
# plot for the stable subset
rand_eob_S <- df_sprt_all_stable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_TW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

p_middle <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TMY weather")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# Timeline plot
df_time <- df_sprt_all_variable %>% 
  filter(seq != "final") %>% 
  bind_rows(df_sprt_all_stable %>% filter(seq != "final")) %>% 
  select(name, seq, n_weeks)

count <- list()
n <- 1
for (i in seq(0, 36, by = 3)){
  
  df_sprt <- df_time %>% 
    filter(seq == "sprt",
           n_weeks <= i)
  
  df_eob <- df_time %>% 
    filter(seq == "eob",
           n_weeks <= i)
  
  df_temp <- df_time %>% 
    filter(seq == "temp", 
           n_weeks <= i)
  
  count[[n]] <- tibble("n_weeks" = i, 
                       "eob" = nrow(df_eob), 
                       "temp" = nrow(df_temp), 
                       "sprt" = nrow(df_sprt))
  
  n <- n + 1
}

count <- bind_rows(count) 

p_bottom <- count %>% 
  ggplot() +
  geom_bar(data = . %>% 
             filter(n_weeks == 24 | n_weeks == 36),
           aes(x = n_weeks, y = eob, fill = "Buildings finishing randomized M&V"), 
           stat = "identity", 
           position = "stack",
           alpha = 0.8, 
           width = 3) +
  geom_line(aes(x = n_weeks, y = temp, color = "Buildings satisfying 80% TMY range"),
            alpha = 0.2) +
  geom_point(aes(x = n_weeks, y = temp, color = "Buildings satisfying 80% TMY range"), 
             size = 1.5) +
  geom_line(aes(x = n_weeks, y = sprt, color = "Buildings satisfying SPRT"),
            alpha = 0.2) +
  geom_point(aes(x = n_weeks, y = sprt, color = "Buildings satisfying SPRT"), 
             size = 1.5, 
             shape = 17) +
  geom_segment(aes(x = 36.5, y = max(eob), xend = 95.5, yend = max(eob)),
               arrow = arrow(length = unit(0.25, "in")),   
               linewidth = 1.1,  
               color = "#fb8072") + 
  annotate(geom = "text", 
           x = 66, 
           y = max(count$eob) - 30, 
           size = 5,
           label = "Excess time by conventional M&V") +
  geom_vline(xintercept = c(12, 24, 36, 48, 96), lty = "dashed", color = "grey80") +
  annotate(geom = "text", 
           x = seq(6, 48, by = 12), 
           y = 200, 
           label = paste0("12-week\nblock"), 
           alpha = 0.5) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = ls_colors) +
  scale_x_continuous(expand = c(0, 0), 
                     limits = c(0, 100),
                     breaks = c(12, 24, 36, 48, 96), 
                     labels = c("12 weeks", "24 weeks", "36 weeks", "1 year", "2 years")) +
  coord_cartesian(ylim = c(0, 650)) +
  labs(x = NULL, 
       y = "Number of buildings", 
       fill = NULL, 
       color = NULL, 
       subtitle = "timeline comparison") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_middle, p_bottom, 
          ncol = 1, nrow = 3,
          heights = c(1.2, 1.2, 1),
          labels = c("a)", "b)", "c)"),
          common.legend = T, 
          legend = "bottom") +
  plot_annotation(title = "M&V overall results comparison")

ggsave(filename = "abstract.png", path = fig_path, units = "in", height = 12, width = 10, dpi = 300)





#### ABS ####
# plot for the stable subset
p1 <- df_MW_S %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{S_building} stable buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the variable subset
p2 <- df_MW_V %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{V_building} variable buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# plot for combined dataset
p3 <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_top <- ggarrange(p3, p1, p2,
                   ncol = 3, nrow = 1,
                   labels = c("a)", "b)", "c)"),
                   align = "hv",
                   legend="none") +
  plot_annotation(subtitle = "with measured weather conditions")

# TMY
# plot for the stable subset
p1 <- df_TW_S %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{S_building} stable buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the variable subset
p2 <- df_TW_V %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{V_building} variable buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for combined dataset
p3 <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_bottom <- ggarrange(p3, p1, p2,
                      ncol = 3, nrow = 1,
                      labels = c("d)", "e)", "f)"),
                      align = "hv",
                      legend="none") +
  plot_annotation(subtitle = "with typical meteorological weather")

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "abs.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)





#### ABS-NULL ####
rand_eob_S <- df_seq_FS_null_stable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = abs(0 - FS), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_null_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_FS_null_variable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = abs(0 - FS), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_null_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_null_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1.5, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       title = "No-saving M&V accuracy comparison", 
       subtitle = str_glue("All {A_building} buildings wither measured weather")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "abs_null.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)





#### ABS-CONT ####
rand_cont_S <- df_cont_stable %>% 
  select(name, site, cont_fs) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - cont_fs),
         method = "rand_fs") %>% 
  select(name, method, diff)

conv_S <- df_FS_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_cont_S, conv_S, rand_final_S)

# plot for the variable subset
rand_cont_V <- df_cont_variable %>% 
  select(name, site, cont_fs) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - cont_fs),
         method = "rand_fs") %>% 
  select(name, method, diff)

conv_V <- df_FS_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_cont_V, conv_V, rand_final_V)

# plot for combined dataset
rand_cont_A <- bind_rows(rand_cont_V, rand_cont_S) %>% 
  select(name, rand_cont_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_cont_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_cont_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_cont_diff" = "Randomized\n(20/80)", "rand_final_diff" = "Randomized\n(50/50)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       title = "M&V accuracy comparison",
       subtitle = str_glue("All {A_building} buildings with measured weather conditions")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "abs_cont.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)






#### ABS-INT ####
interval_1_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_7), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_7), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY version
interval_1_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_2)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_3)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_7)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_7), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_2)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_3)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_7)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_tmy_7), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_TW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>%
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TMY weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "abs_interval.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)

# sprt sequence plot for different sampling intervals
interval_1_S <- df_seq_FS_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

# variable set
interval_1_V <- df_seq_FS_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - FS), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>%
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# with tmy version
interval_1_S <- df_sprt_all_stable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

# variable set
interval_1_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_TW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>%
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TMY weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "abs_interval_sprt.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)

# null results
interval_1_S <- df_FS_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_7), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_null_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_2), 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_3), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - interval_fs_7), 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings",
       title = "No-saving M&V accuracy comparison", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "abs_interval_null.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)







#### MEAN ####
# plot for the stable subset
rand_eob_S <- df_seq_FS_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

p1 <- df_MW_S %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{S_building} stable buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the variable subset
rand_eob_V <- df_seq_FS_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site"))%>% 
  mutate(diff = savings - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

p2 <- df_MW_V %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{V_building} variable buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

p3 <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_top <- ggarrange(p3, p1, p2,
                   ncol = 3, nrow = 1,
                   labels = c("a)", "b)", "c)"),
                   align = "hv",
                   legend="none") +
  plot_annotation(subtitle = "Accuracy comparison with measured weather conditions")

# TMY 
# plot for the stable subset
rand_eob_S <- df_sprt_all_stable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

p1 <- df_TW_S %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{S_building} stable buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the variable subset
rand_eob_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

p2 <- df_TW_V %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = str_glue("{V_building} variable buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_TW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

p3 <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_bottom <- ggarrange(p3, p1, p2,
                      ncol = 3, nrow = 1,
                      labels = c("d)", "e)", "f)"),
                      align = "hv",
                      legend="none") +
  plot_annotation(subtitle = "Accuracy comparison with typical meteorological weather")

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "mean.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)






#### MEAN-NULL ####
rand_eob_S <- df_seq_FS_null_stable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = 0 - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_FS_null_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_FS_null_variable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = 0 - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_FS_null_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_null_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_eob_V, conv_V, rand_final_V)

# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_V, rand_eob_S) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "conv_diff" = "Conventional", 
                                "rand_eob_diff" = "Randomized", 
                                "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings",
       title = "No-saving M&V accuracy comparison", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "mean_null.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)







#### MEAN-CONT ####
rand_cont_S <- df_cont_stable %>% 
  select(name, site, cont_fs) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - cont_fs,
         method = "rand_fs") %>% 
  select(name, method, diff)

conv_S <- df_FS_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_cont_S, conv_S, rand_final_S)

# plot for the variable subset
rand_cont_V <- df_cont_variable %>% 
  select(name, site, cont_fs) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - cont_fs,
         method = "rand_fs") %>% 
  select(name, method, diff)

conv_V <- df_FS_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_V <- rbind(rand_cont_V, conv_V, rand_final_V)

# plot for combined dataset
rand_cont_A <- bind_rows(rand_cont_V, rand_cont_S) %>% 
  select(name, rand_cont_diff = diff)

conv_A <- bind_rows(conv_V, conv_S) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rand_cont_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_cont_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_cont_diff" = "Randomized\n(20/80)", "rand_final_diff" = "Randomized\n(50/50)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       title = "M&V accuracy comparison",
       subtitle = str_glue("All {A_building} buildings with measured weather conditions")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "mean_cont.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)





#### MEAN-INT ####
interval_1_S <- df_FS_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_stable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_7, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_variable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_7, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY version
interval_1_S <- df_FS_tmy_stable %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_2)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_3)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_stable %>% 
  select(c(name, site, interval_tmy_7)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_7, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_tmy_variable %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_2)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_3)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_variable %>% 
  select(c(name, site, interval_tmy_7)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_tmy_7, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_TW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TMY weather")) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "mean_interval.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)

# sprt sequence plot for different sampling intervals
interval_1_S <- df_seq_FS_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_seq_interval_FS_stable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

# variable set
interval_1_V <- df_seq_FS_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_seq_interval_FS_variable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - FS, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# with tmy version
interval_1_S <- df_sprt_all_stable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_seq_interval_nm_stable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

# variable set
interval_1_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, sprt, n_weeks)) %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 2 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 3 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_seq_interval_nm_variable %>% 
  filter(interval == 7 & seq == "eob") %>% 
  left_join(df_FS_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_TW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_TW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TMY weather")) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          legend = "none") +
  plot_annotation(title = "M&V accuracy comparison")

ggsave(filename = "mean_interval_sprt.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)

# null results
interval_1_S <- df_FS_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_S <- df_interval_null_stable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_null_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_7, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_S <- rbind(interval_1_S, interval_2_S, interval_3_S, interval_7_S)

interval_1_V <- df_FS_null_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_interval_1") %>% 
  select(name, method, diff)

interval_2_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_2)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_2, 
         method = "rand_interval_2") %>% 
  select(name, method, diff)

interval_3_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_3)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_3, 
         method = "rand_interval_3") %>% 
  select(name, method, diff)

interval_7_V <- df_interval_null_variable %>% 
  select(c(name, site, interval_fs_7)) %>% 
  left_join(df_FS_null_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - interval_fs_7, 
         method = "rand_interval_7") %>% 
  select(name, method, diff)

df_MW_V <- rbind(interval_1_V, interval_2_V, interval_3_V, interval_7_V)

# plot for combined dataset
interval_1_A <- bind_rows(interval_1_V, interval_1_S) %>% 
  select(name, interval_1_diff = diff)

interval_2_A <- bind_rows(interval_2_V, interval_2_S) %>% 
  select(name, interval_2_diff = diff)

interval_3_A <- bind_rows(interval_3_V, interval_3_S) %>% 
  select(name, interval_3_diff = diff)

interval_7_A <- bind_rows(interval_7_V, interval_7_S) %>% 
  select(name, interval_7_diff = diff)

df_MW_A <- interval_1_A %>% 
  left_join(interval_2_A, by = "name") %>% 
  left_join(interval_3_A, by = "name") %>% 
  left_join(interval_7_A, by = "name") %>% 
  pivot_longer(c(interval_1_diff, interval_2_diff, interval_3_diff, interval_7_diff), names_to = "method", values_to = "diff")

df_MW_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1_diff" = "Daily\nsampling", 
                                "interval_2_diff" = "2-day\nsampling", 
                                "interval_3_diff" = "3-day\nsampling", 
                                "interval_7_diff" = "7-day\nsampling")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = method)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(method) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = method, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Error in fractional savings", 
       title = "No-saving M&V accuracy comparison",
       subtitle = str_glue("All {A_building} buildings with measure weather conditions")) +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "mean_interval_null.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)





#### TIME-INT ####
# timeline to finish sequential test
interval_1_count <- bind_rows(df_sprt_all_stable, df_sprt_all_variable) %>% 
  select(c(seq, n_weeks)) %>% 
  filter(seq != "final") %>% 
  mutate(seq = as.factor(seq)) %>% 
  mutate(method = "interval_1")

interval_2_count <- bind_rows(df_seq_interval_FS_stable, df_seq_interval_FS_variable) %>% 
  drop_na() %>% 
  filter(seq != "final" & interval == 2) %>% 
  select(c(seq, n_weeks)) %>% 
  mutate(seq = as.factor(seq)) %>% 
  mutate(method = "interval_2")

interval_3_count <- bind_rows(df_seq_interval_FS_stable, df_seq_interval_FS_variable) %>% 
  drop_na() %>% 
  filter(seq != "final" & interval == 3) %>% 
  select(c(seq, n_weeks)) %>% 
  mutate(seq = as.factor(seq)) %>% 
  mutate(method = "interval_3")

interval_7_count <- bind_rows(df_seq_interval_FS_stable, df_seq_interval_FS_variable) %>% 
  drop_na() %>% 
  filter(seq != "final" & interval == 7) %>% 
  select(c(seq, n_weeks)) %>% 
  mutate(seq = as.factor(seq)) %>% 
  mutate(method = "interval_7")

df_count <- bind_rows(interval_1_count, interval_2_count, interval_3_count, interval_7_count)

df_count %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, 
                                "interval_1" = "Daily\nsampling", 
                                "interval_2" = "2-day\nsampling", 
                                "interval_3" = "3-day\nsampling", 
                                "interval_7" = "7-day\nsampling"),
         seq = recode_factor(seq, "sprt" = "SPRT", "temp" = "80% TMY range", "eob" = "M&V complete")) %>% 
  ggplot() +
  geom_boxplot(aes(x = seq, y = n_weeks, fill = method)) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = seq(0, 50, by = 12)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Satisfied stopping criterion", 
       y = "Number of weeks", 
       fill = NULL, 
       title = "Overall M&V timeline by different sampling intervals") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "timeline_interval_sprt.png", path = fig_path, units = "in", height = 6, width = 9, dpi = 300)


