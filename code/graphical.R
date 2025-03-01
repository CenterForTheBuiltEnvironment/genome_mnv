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
               "2-day\nsampling" = "#66c2a4",
               "3-day\nsampling" = "#66c2a4",
               "7-day\nsampling" = "#66c2a4",
               "Dropped" = "#41ae76",
               "Kept" = "#006d2c",
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
df_weather_stable <- read_rds(paste0(readfile_stable, "df_weather.rds"))
df_sprt_all_stable <- read_rds(paste0(readfile_stable, "df_sprt_all.rds"))
df_seq_fs_stable <- read_rds(paste0(readfile_stable, "df_seq_fs.rds"))
df_fs_stable <- read_rds(paste0(readfile_stable, "df_fs.rds"))
df_cont_stable <- read_rds(paste0(readfile_stable, "df_cont.rds"))
df_fs_null_stable <- read_rds(paste0(readfile_stable, "df_fs_null.rds"))
df_seq_fs_null_stable <- read_rds(paste0(readfile_stable, "df_seq_fs_null.rds"))
df_model_acc_stable <- read_rds(paste0(readfile_stable, "df_model_acc.rds"))
df_fs_tmy_stable <- read_rds(paste0(readfile_stable, "df_fs_tmy.rds"))
df_fs_tmy_null_stable <- read_rds(paste0(readfile_stable, "df_fs_tmy_null.rds"))
df_interval_drop_stable <- read_rds(paste0(readfile_stable, "df_interval_drop.rds"))
df_interval_keep_stable <- read_rds(paste0(readfile_stable, "df_interval_keep.rds"))
df_seq_interval_fs_drop_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_fs_drop.rds"))
df_seq_interval_fs_keep_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_fs_keep.rds"))
df_seq_interval_nm_drop_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_nm_drop.rds"))
df_seq_interval_nm_keep_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_nm_keep.rds"))
df_seq_interval_tl_drop_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_tl_drop.rds"))
df_seq_interval_tl_keep_stable <- read_rds(paste0(readfile_stable, "df_seq_interval_tl_keep.rds"))


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
df_weather_variable <- read_rds(paste0(readfile_variable, "df_weather.rds"))
df_sprt_all_variable <- read_rds(paste0(readfile_variable, "df_sprt_all.rds"))
df_seq_fs_variable <- read_rds(paste0(readfile_variable, "df_seq_fs.rds"))
df_fs_variable <- read_rds(paste0(readfile_variable, "df_fs.rds"))
df_cont_variable <- read_rds(paste0(readfile_variable, "df_cont.rds"))
df_fs_null_variable <- read_rds(paste0(readfile_variable, "df_fs_null.rds"))
df_seq_fs_null_variable <- read_rds(paste0(readfile_variable, "df_seq_fs_null.rds"))
df_model_acc_variable <- read_rds(paste0(readfile_variable, "df_model_acc.rds"))
df_fs_tmy_variable <- read_rds(paste0(readfile_variable, "df_fs_tmy.rds"))
df_fs_tmy_null_variable <- read_rds(paste0(readfile_variable, "df_fs_tmy_null.rds"))
df_interval_drop_variable <- read_rds(paste0(readfile_variable, "df_interval_drop.rds"))
df_interval_keep_variable <- read_rds(paste0(readfile_variable, "df_interval_keep.rds"))
df_seq_interval_fs_drop_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_fs_drop.rds"))
df_seq_interval_fs_keep_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_fs_keep.rds"))
df_seq_interval_nm_drop_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_nm_drop.rds"))
df_seq_interval_nm_keep_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_nm_keep.rds"))
df_seq_interval_tl_drop_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_tl_drop.rds"))
df_seq_interval_tl_keep_variable <- read_rds(paste0(readfile_variable, "df_seq_interval_tl_keep.rds"))

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

S_building <- nrow(df_fs_tmy_stable)
V_building <- nrow(df_fs_tmy_variable)
A_building <- S_building + V_building




#### GA ####
rand_eob_S <- df_seq_fs_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - fs), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_fs_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site"))%>% 
  mutate(diff = abs(savings - fs), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_variable %>% 
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
  select(-c(seq, n_weeks)) %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_tmy_stable %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_tmy_stable %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, n_weeks)) %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - annual), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_tmy_variable %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_tmy_variable %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
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
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_accuracy <- ggarrange(p_top, p_middle, 
                        ncol = 1, nrow = 2,
                        labels = c("a)", "b)"),
                        common.legend = T, 
                        legend = "bottom")
  
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

p_timeline <- count %>% 
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

ggarrange(p_accuracy, p_timeline, 
          ncol = 1, nrow = 2,
          heights = c(2.2, 1),
          labels = c("", "c)"),
          common.legend = F, 
          legend = "bottom") +
  plot_annotation(title = "Overall comparison of conventional and randomized M&V")

ggsave(filename = "abstract.png", path = fig_path, units = "in", height = 10, width = 12, dpi = 300)





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
  plot_annotation(title = "Absolute error in savings estimation by different M&V methods")

ggsave(filename = "abs.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)





#### ABS-NULL ####
rand_eob_S <- df_seq_fs_null_stable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = abs(0 - fs), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_null_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_fs_null_variable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = abs(0 - fs), 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_null_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - conv), 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_null_variable %>% 
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
       y = "Absolute effect in fractional savings", 
       title = "Absolute effect of non-routine events when using different M&V methods", 
       subtitle = str_glue("All {A_building} buildings wither measured weather")) +
  coord_cartesian(ylim = c(0, 23)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "abs_null.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)





#### ABS-CONT ####
rand_cont_S <- df_cont_stable %>% 
  select(-cont_tmy) %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - cont_fs)) %>% 
  select(name, diff, ratio)

rand_final_S <- df_fs_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         ratio = 0) %>% 
  select(name, diff, ratio)

df_MW_S <- rbind(rand_cont_S, rand_final_S)

# plot for the variable subset
rand_cont_V <- df_cont_variable %>% 
  select(-cont_tmy) %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = abs(savings - cont_fs)) %>% 
  select(name, diff, ratio)

rand_final_V <- df_fs_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = abs(true - rand), 
         ratio = 0) %>% 
  select(name, diff, ratio)

df_MW_V <- rbind(rand_cont_V, rand_final_V)

# plot for combined dataset
rand_cont_A <- bind_rows(rand_cont_V, rand_cont_S) %>% 
  select(name, rand_cont_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rbind(df_MW_S, df_MW_V)

df_MW_A %>% 
  mutate(ratio = as.factor(ratio), 
         ratio = recode_factor(ratio, "0" = "Randomized\n(50/50)", "1" = "Randomized\n(20/80)", "2" = "Randomized\n(10/90)")) %>% 
  ggplot(aes(x = ratio, y = diff, fill = "Randomized")) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = ratio)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% 
              group_by(ratio) %>% 
              summarise(mean = mean(diff)) %>% 
              ungroup(), 
            aes(x = ratio, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       title = "Absolute error in savings estimation by different sampling ratio",
       subtitle = str_glue("All {A_building} buildings with measured weather conditions")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "abs_cont.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)






#### ABS-INT ####
# dropped version
df_interval_keep_stable %<>% 
  left_join(df_fs_stable %>% filter(method == "rand"), by = c("name", "site")) %>% 
  select(-c("scenario", "method")) %>% 
  rename("interval_fs_1" = "savings") %>% 
  left_join(df_fs_tmy_stable, by = c("name", "site")) %>% 
  select(-conv) %>% 
  rename("interval_tmy_1" = "rand")

df_interval_keep_variable %<>% 
  left_join(df_fs_variable %>% filter(method == "rand"), by = c("name", "site")) %>% 
  select(-c("scenario", "method")) %>% 
  rename("interval_fs_1" = "savings") %>% 
  left_join(df_fs_tmy_variable, by = c("name", "site")) %>% 
  select(-conv) %>% 
  rename("interval_tmy_1" = "rand")

df_drop <- df_interval_drop_stable %>% 
  bind_rows(df_interval_drop_variable) %>% 
  select(name, contains("fs")) %>% 
  rename("drop_1" = "interval_fs_1", 
         "drop_2" = "interval_fs_2", 
         "drop_3" = "interval_fs_3", 
         "drop_7" = "interval_fs_7")  

df_keep <- df_interval_keep_stable %>% 
  bind_rows(df_interval_keep_variable) %>% 
  select(name, contains("fs")) %>% 
  rename("keep_1" = "interval_fs_1",
         "keep_2" = "interval_fs_2", 
         "keep_3" = "interval_fs_3", 
         "keep_7" = "interval_fs_7")

df_MW_A <- df_drop %>% 
  left_join(df_keep, by = c("name")) %>%
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(drop_1 = abs(savings - drop_1),
         drop_2 = abs(savings - drop_2), 
         drop_3 = abs(savings - drop_3), 
         drop_7 = abs(savings - drop_7), 
         keep_1 = abs(savings - keep_1),
         keep_2 = abs(savings - keep_2), 
         keep_3 = abs(savings - keep_3), 
         keep_7 = abs(savings - keep_7)) %>% 
  pivot_longer(cols = contains("drop") | contains("keep"),  
               names_to = c("group", "interval"),
               names_pattern = "(drop|keep)_(\\d+)",
               values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75), color = "white") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY version
df_drop <- df_interval_drop_stable %>% 
  bind_rows(df_interval_drop_variable) %>% 
  select(name, contains("tmy")) %>% 
  rename("drop_1" = "interval_tmy_1",
         "drop_2" = "interval_tmy_2", 
         "drop_3" = "interval_tmy_3", 
         "drop_7" = "interval_tmy_7")  

df_keep <- df_interval_keep_stable %>% 
  bind_rows(df_interval_keep_variable) %>% 
  select(name, contains("tmy")) %>% 
  rename("keep_1" = "interval_tmy_1",
         "keep_2" = "interval_tmy_2", 
         "keep_3" = "interval_tmy_3", 
         "keep_7" = "interval_tmy_7")

df_TW_A <- df_drop %>% 
  left_join(df_keep, by = c("name")) %>%
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(drop_1 = abs(savings - drop_1),
         drop_2 = abs(savings - drop_2), 
         drop_3 = abs(savings - drop_3), 
         drop_7 = abs(savings - drop_7), 
         keep_1 = abs(savings - keep_1),
         keep_2 = abs(savings - keep_2), 
         keep_3 = abs(savings - keep_3), 
         keep_7 = abs(savings - keep_7)) %>% 
  pivot_longer(cols = contains("drop") | contains("keep"),  
               names_to = c("group", "interval"),
               names_pattern = "(drop|keep)_(\\d+)",
               values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75), color = "white") +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(0, 10)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          common.legend = T, 
          legend = "bottom") +
  plot_annotation(title = "Absolute error in savings estimation of ranodmized M&V by different sampling intervals")

ggsave(filename = "abs_interval.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)

# sprt sequence plot for different sampling intervals
df_seq_interval_fs_keep_stable %<>% 
  select(-c(sprt, temp, final)) %>% 
  bind_rows(df_seq_fs_stable %>% 
              filter(seq == "eob") %>% 
              select(-seq) %>% 
              mutate(interval = 1) %>% 
              rename(eob = "fs"))  

df_seq_interval_fs_keep_variable %<>% 
  select(-c(sprt, temp, final)) %>% 
  bind_rows(df_seq_fs_variable %>% 
              filter(seq == "eob") %>% 
              select(-seq) %>% 
              mutate(interval = 1) %>% 
              rename(eob = "fs")) 

df_drop <- df_seq_interval_fs_drop_stable %>% 
  bind_rows(df_seq_interval_fs_drop_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "drop")

df_keep <- df_seq_interval_fs_keep_stable %>% 
  bind_rows(df_seq_interval_fs_keep_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "keep")

df_MW_A <- df_drop %>% 
  bind_rows(df_keep) %>% 
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(diff = abs(savings - eob)) %>% 
  filter(is.finite(diff))

p_top <- df_MW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, inf.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# with tmy version
df_seq_interval_nm_keep_stable %<>% 
  select(-c(sprt, temp, final)) %>% 
  bind_rows(df_sprt_all_stable %>% 
              filter(seq == "eob") %>% 
              select(-c(seq, n_weeks)) %>% 
              mutate(interval = 1) %>% 
              rename(eob = "annual")) 

df_seq_interval_nm_keep_variable %<>% 
  select(-c(sprt, temp, final)) %>% 
  bind_rows(df_sprt_all_variable %>% 
              filter(seq == "eob") %>% 
              select(-c(seq, n_weeks)) %>% 
              mutate(interval = 1) %>% 
              rename(eob = "annual")) 

df_drop <- df_seq_interval_nm_drop_stable %>% 
  bind_rows(df_seq_interval_nm_drop_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "drop")

df_keep <- df_seq_interval_nm_keep_stable %>% 
  bind_rows(df_seq_interval_nm_keep_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "keep")

df_TW_A <- df_drop %>% 
  bind_rows(df_keep) %>% 
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(diff = abs(savings - eob)) %>% 
  filter(is.finite(diff))

p_bottom <- df_TW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(0, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          common.legend = T, 
          legend = "bottom") +
  plot_annotation(title = "Absolute error in savings estimation of randomized M&V after satisfying all stopping criteria")

ggsave(filename = "abs_interval_sprt.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)






#### MEAN ####
# plot for the stable subset
rand_eob_S <- df_seq_fs_stable %>% 
  filter(seq == "eob") %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - fs, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_fs_variable %>% 
  filter(seq == "eob") %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site"))%>% 
  mutate(diff = savings - fs, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_variable %>% 
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
       y = "Mean error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-20, 20)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY
# plot for the stable subset
rand_eob_S <- df_sprt_all_stable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, n_weeks)) %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_tmy_stable %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_tmy_stable %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_TW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_sprt_all_variable %>% 
  filter(seq == "eob") %>% 
  select(-c(seq, n_weeks)) %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - annual, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_tmy_variable %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_tmy_variable %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - rand, 
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
       y = "Mean error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(-20, 20)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_middle, 
          ncol = 1, nrow = 2,
          labels = c("a)", "b)")) +
  plot_annotation(title = "Mean error in savings estimation by different M&V methods")

ggsave(filename = "mean.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)






#### MEAN-NULL ####
rand_eob_S <- df_seq_fs_null_stable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = 0 - fs, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_S <- df_fs_null_stable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_S <- df_fs_null_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_MW_S <- rbind(rand_eob_S, conv_S, rand_final_S)

# plot for the variable subset
rand_eob_V <- df_seq_fs_null_variable %>% 
  filter(seq == "eob") %>% 
  mutate(diff = 0 - fs, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_V <- df_fs_null_variable %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_V <- df_fs_null_variable %>% 
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
                     breaks = breaks_pretty(n = 6), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Mean effect in fractional savings",
       title = "Mean effect of non-routine events when using different M&V methods", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "mean_null.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)







#### MEAN-CONT ####
rand_cont_S <- df_cont_stable %>% 
  select(-cont_tmy) %>% 
  left_join(df_fs_stable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - cont_fs) %>% 
  select(name, diff, ratio)

rand_final_S <- df_fs_stable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         ratio = 0) %>% 
  select(name, diff, ratio)

df_MW_S <- rbind(rand_cont_S, rand_final_S)

# plot for the variable subset
rand_cont_V <- df_cont_variable %>% 
  select(-cont_tmy) %>% 
  left_join(df_fs_variable %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(diff = savings - cont_fs) %>% 
  select(name, diff, ratio)

rand_final_V <- df_fs_variable %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         ratio = 0) %>% 
  select(name, diff, ratio)

df_MW_V <- rbind(rand_cont_V, rand_final_V)

# plot for combined dataset
rand_cont_A <- bind_rows(rand_cont_V, rand_cont_S) %>% 
  select(name, rand_cont_diff = diff)

rand_final_A <- bind_rows(rand_final_V, rand_final_S) %>% 
  select(name, rand_final_diff = diff)

df_MW_A <- rbind(df_MW_S, df_MW_V)

df_MW_A %>% 
  mutate(ratio = as.factor(ratio), 
         ratio = recode_factor(ratio, "0" = "Randomized\n(50/50)", "1" = "Randomized\n(20/80)", "2" = "Randomized\n(10/90)")) %>% 
  ggplot(aes(x = ratio, y = diff, fill = "Randomized")) +
  geom_jitter(width = 0.2, alpha = 0.8, size = 0.5) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = ratio)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(ratio) %>% summarise(mean = mean(diff)) %>% ungroup(), 
            aes(x = ratio, y = mean, label = paste0(round(mean, digits = 1), " %"))) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       x = NULL, 
       y = "Mean error in fractional savings", 
       title = "Mean error in savings estimation by different sampling ratio",
       subtitle = str_glue("All {A_building} buildings with measured weather conditions")) +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = "mean_cont.png", path = fig_path, units = "in", height = 6, width = 12, dpi = 300)






#### MEAN-INT ####
# dropped version
df_drop <- df_interval_drop_stable %>% 
  bind_rows(df_interval_drop_variable) %>% 
  select(name, contains("fs")) %>% 
  rename("drop_1" = "interval_fs_1", 
         "drop_2" = "interval_fs_2", 
         "drop_3" = "interval_fs_3", 
         "drop_7" = "interval_fs_7")  

df_keep <- df_interval_keep_stable %>% 
  bind_rows(df_interval_keep_variable) %>% 
  select(name, contains("fs")) %>% 
  rename("keep_1" = "interval_fs_1",
         "keep_2" = "interval_fs_2", 
         "keep_3" = "interval_fs_3", 
         "keep_7" = "interval_fs_7")

df_MW_A <- df_drop %>% 
  left_join(df_keep, by = c("name")) %>%
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(drop_1 = savings - drop_1,
         drop_2 = savings - drop_2, 
         drop_3 = savings - drop_3, 
         drop_7 = savings - drop_7, 
         keep_1 = savings - keep_1,
         keep_2 = savings - keep_2, 
         keep_3 = savings - keep_3, 
         keep_7 = savings - keep_7) %>% 
  pivot_longer(cols = contains("drop") | contains("keep"),  
               names_to = c("group", "interval"),
               names_pattern = "(drop|keep)_(\\d+)",
               values_to = "diff")

p_top <- df_MW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group)) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Mean error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# TMY version
df_drop <- df_interval_drop_stable %>% 
  bind_rows(df_interval_drop_variable) %>% 
  select(name, contains("tmy")) %>% 
  rename("drop_1" = "interval_tmy_1",
         "drop_2" = "interval_tmy_2", 
         "drop_3" = "interval_tmy_3", 
         "drop_7" = "interval_tmy_7")  

df_keep <- df_interval_keep_stable %>% 
  bind_rows(df_interval_keep_variable) %>% 
  select(name, contains("tmy")) %>% 
  rename("keep_1" = "interval_tmy_1",
         "keep_2" = "interval_tmy_2", 
         "keep_3" = "interval_tmy_3", 
         "keep_7" = "interval_tmy_7")

df_TW_A <- df_drop %>% 
  left_join(df_keep, by = c("name")) %>%
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(drop_1 = savings - drop_1,
         drop_2 = savings - drop_2, 
         drop_3 = savings - drop_3, 
         drop_7 = savings - drop_7, 
         keep_1 = savings - keep_1,
         keep_2 = savings - keep_2, 
         keep_3 = savings - keep_3, 
         keep_7 = savings - keep_7) %>% 
  pivot_longer(cols = contains("drop") | contains("keep"),  
               names_to = c("group", "interval"),
               names_pattern = "(drop|keep)_(\\d+)",
               values_to = "diff")

p_bottom <- df_TW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Mean error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          common.legend = T, 
          legend = "bottom") +
  plot_annotation(title = "Mean error in savings estimation of randomized M&V by different sampling intervals")

ggsave(filename = "mean_interval.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)

# sprt sequence plot for different sampling intervals
df_drop <- df_seq_interval_fs_drop_stable %>% 
  bind_rows(df_seq_interval_fs_drop_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "drop")

df_keep <- df_seq_interval_fs_keep_stable %>% 
  bind_rows(df_seq_interval_fs_keep_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "keep")

df_MW_A <- df_drop %>% 
  bind_rows(df_keep) %>% 
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(diff = (savings - eob)) %>% 
  filter(is.finite(diff))

p_top <- df_MW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, inf.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with measured weather")) +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# with tmy version
df_drop <- df_seq_interval_nm_drop_stable %>% 
  bind_rows(df_seq_interval_nm_drop_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "drop")

df_keep <- df_seq_interval_nm_keep_stable %>% 
  bind_rows(df_seq_interval_nm_keep_variable) %>% 
  select(name, eob, interval) %>% 
  mutate(group = "keep")

df_TW_A <- df_drop %>% 
  bind_rows(df_keep) %>% 
  left_join(rbind(df_fs_stable %>% 
                    filter(scenario == "ref" & method == "true"), 
                  df_fs_variable %>% 
                    filter(scenario == "ref" & method == "true")), by = "name") %>% 
  mutate(diff = (savings - eob)) %>% 
  filter(is.finite(diff))

p_bottom <- df_TW_A %>% 
  mutate(interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = diff, fill = group)) +
  geom_jitter(alpha = 0.8, 
              size = 0.5, 
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.alpha = 0, coef = 0, fill = "#00000000", aes(color = group), width = 0.8) +
  geom_hline(yintercept = 0, color = "#fb8072", linewidth = 1, lty = "dashed") +
  geom_text(data = . %>% group_by(group, interval) %>% summarise(mean = mean(diff, na.rm = T)) %>% ungroup(), 
            aes(x = interval, y = mean, group = group,
                label = paste0(round(mean, digits = 1), " %")), 
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = c("grey80", "grey80")) + 
  labs(fill = NULL, 
       color = NULL, 
       x = NULL, 
       y = "Absolute error in fractional savings", 
       subtitle = str_glue("All {A_building} buildings with TOWT model and TMY weather")) +
  coord_cartesian(ylim = c(-15, 15)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2, 
          align = "hv", 
          common.legend = T, 
          legend = "bottom") +
  plot_annotation(title = "Absolute error in savings estimation of randomized M&V after satisfying all stopping criteria")

ggsave(filename = "mean_interval_sprt.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)





#### TIME-INT ####
# timeline to finish sequential test
df_drop <- bind_rows(df_seq_interval_tl_drop_stable, df_seq_interval_tl_drop_variable) %>%
  mutate(group = "drop", 
         temp = ifelse(base_temp >= interv_temp, base_temp, interv_temp)) %>% 
  pivot_longer(c(sprt, temp, eob, final), names_to = "seq", values_to = "n_weeks") %>% 
  select(name, seq, n_weeks, interval, group)

interval_1_keep <- bind_rows(df_sprt_all_stable, df_sprt_all_variable) %>%
  select(c(name, seq, n_weeks)) %>%
  filter(seq != "final") %>%
  mutate(seq = as.factor(seq)) %>%
  mutate(interval = 1, 
         group = "keep")
  
df_keep <- bind_rows(df_seq_interval_tl_keep_stable, df_seq_interval_tl_keep_variable) %>%
  mutate(group = "keep", 
         temp = ifelse(base_temp >= interv_temp, base_temp, interv_temp)) %>% 
  pivot_longer(c(sprt, temp, eob, final), names_to = "seq", values_to = "n_weeks") %>% 
  select(name, seq, n_weeks, interval, group) %>% 
  bind_rows(interval_1_keep)

df_MW_A <- bind_rows(df_keep, df_drop) %>% 
  filter(seq != "final")

p_top <- df_MW_A %>% 
  filter(seq != "eob") %>% 
  mutate(seq = as.factor(seq), 
         seq = recode_factor(seq, 
                             "sprt" = "Buildings satisfying SPRT", 
                             "temp" = "Buildings satisfying temperature range", 
                             "eob" = "Buildings finishing all criteria"), 
         interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  ggplot(aes(x = interval, y = n_weeks, fill = group)) +
  geom_jitter(alpha = 0.8,
              size = 0.5,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  geom_lv(k = 4, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_boxplot(outlier.shape = NA, coef = 0, color = "grey50", width = 0.8) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = ls_colors) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 50, by = 12), 
                     labels = number_format(suffix = " weeks")) +
  coord_cartesian(ylim = c(0, 50)) +
  facet_wrap(~seq, nrow = 1) +
  labs(x = NULL,
       y = NULL,
       fill = NULL,
       subtitle = "Satisfying EACH stopping criterion") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_bottom <- df_MW_A %>% 
  filter(seq == "eob") %>% 
  mutate(seq = as.factor(seq), 
         seq = recode_factor(seq, 
                             "eob" = "Buildings finishing all criteria"), 
         interval = as.factor(interval), 
         interval = recode_factor(interval, 
                                  "1" = "1-day\nsampling",
                                  "2" = "2-day\nsampling", 
                                  "3" = "3-day\nsampling", 
                                  "7" = "7-day\nsampling"), 
         group = as.factor(group), 
         group = recode_factor(group, "drop" = "Dropped", "keep" = "Kept")) %>% 
  group_by(n_weeks, group, interval) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(group, interval) %>% 
  mutate(n = n, 
         perc = n / sum(n) * 100) %>% 
  ungroup() %>% 
  ggplot(aes(x = n_weeks, y = perc, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~interval, nrow = 1) +
  scale_fill_manual(values = ls_colors) +
  scale_color_manual(values = ls_colors) +
  scale_x_continuous(breaks = c(12, 24, 36, 48), 
                     labels = number_format(suffix = "\nweeks")) +
  scale_y_continuous(breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = " %")) +
  labs(x = NULL,
       y = "Percentage of buildings finishing randomized M&V",
       fill = NULL,
       subtitle = "Satisfying ALL stopping criteria") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p_top, p_bottom, 
          ncol = 1, 
          common.legend = T, 
          legend = "bottom") +
  plot_annotation("Overall M&V timeline by different sampling intervals")

ggsave(filename = "timeline_interval_sprt.png", path = fig_path, units = "in", height = 6, width = 9, dpi = 300)
