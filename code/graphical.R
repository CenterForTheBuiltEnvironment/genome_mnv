#### LIBRARIES ####
require(pacman)

# load packages using pacman
pacman::p_load(tidyverse, lubridate, here, stats, zoo, scales, ggpubr, patchwork, RColorBrewer)

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
               "Conventional M&V" = "grey70", 
               "Randomized M&V" = "#66c2a4", 
               "Buildings satisfying all criteria" = "#66c2a4", 
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

# Tidy dataset
readfile_tidy <- str_glue("../readfiles/tidy/")
fig_path = "../figs/manuscript/"

df_energy_tidy <- read_rds(paste0(readfile_tidy, "df_energy.rds"))
df_meta_tidy <- read_rds(paste0(readfile_tidy, "df_meta.rds"))
df_weather_tidy <- read_rds(paste0(readfile_tidy, "df_weather.rds"))
df_sprt_all_tidy <- read_rds(paste0(readfile_tidy, "df_sprt_all.rds"))
df_seq_FS_tidy <- read_rds(paste0(readfile_tidy, "df_seq_FS.rds"))
df_NRE_occ <- read_rds(paste0(readfile_tidy, "df_NRE_occ.rds"))
df_MD_tidy <- read_rds(paste0(readfile_tidy, "df_MD.rds"))
df_FS_tidy <- read_rds(paste0(readfile_tidy, "df_FS.rds"))
df_eui_tidy <- read_rds(paste0(readfile_tidy, "df_eui.rds"))
df_cont_MD_tidy <- read_rds(paste0(readfile_tidy, "df_cont_MD.rds"))
df_cont_FS_tidy <- read_rds(paste0(readfile_tidy, "df_cont_FS.rds"))
df_model_acc_tidy <- read_rds(paste0(readfile_tidy, "df_model_acc.rds"))

all_sites_tidy <- df_energy_tidy %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types_tidy <- df_energy_tidy %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

all_names_tidy <- df_energy_tidy %>%
  select(name) %>%
  distinct(name)

# Messy dataset
readfile_messy <- str_glue("../readfiles/messy/")

df_energy_messy <- read_rds(paste0(readfile_messy, "df_energy.rds"))
df_meta_messy <- read_rds(paste0(readfile_messy, "df_meta.rds"))
df_weather_messy <- read_rds(paste0(readfile_messy, "df_weather.rds"))
df_sprt_all_messy <- read_rds(paste0(readfile_messy, "df_sprt_all.rds"))
df_seq_FS_messy <- read_rds(paste0(readfile_messy, "df_seq_FS.rds"))
df_MD_messy <- read_rds(paste0(readfile_messy, "df_MD.rds"))
df_FS_messy <- read_rds(paste0(readfile_messy, "df_FS.rds"))
df_eui_messy <- read_rds(paste0(readfile_messy, "df_eui.rds"))
df_cont_MD_messy <- read_rds(paste0(readfile_messy, "df_cont_MD.rds"))
df_cont_FS_messy <- read_rds(paste0(readfile_messy, "df_cont_FS.rds"))
df_FS_nsprt <- read_rds(paste0(readfile_messy, "df_FS_nsprt.rds"))
df_MD_nsprt <- read_rds(paste0(readfile_messy, "df_MD_nsprt.rds"))
df_seq_FS_nsprt <- read_rds(paste0(readfile_messy, "df_seq_FS_nsprt.rds"))
df_model_acc_messy <- read_rds(paste0(readfile_messy, "df_model_acc.rds"))

all_sites_messy <- df_energy_messy %>%
  select(site) %>%
  distinct() %>%
  arrange(site)

all_types_messy <- df_energy_messy %>%
  select(type) %>%
  mutate(type = as.factor(type)) %>%
  distinct()

all_names_messy <- df_energy_messy %>%
  select(name) %>%
  distinct(name)

# read functions
function_path <- "../functions/"
source(paste0(function_path, "model_fit.R"))
source(paste0(function_path, "model_pred.R"))
source(paste0(function_path, "prepost_plot.R"))




#### ABSTRACT ####
# plot for the tidy subset
rand_eob_T <- df_seq_FS_tidy %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_tidy %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(rand_diff_eob = savings - FS) %>% 
  select(name, site, rand_diff_eob)

conv_T <- df_FS_tidy %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(conv_diff = true - conv) %>% 
  select(name, site, conv_diff)

rand_final_T <- df_seq_FS_tidy %>% 
  filter(seq == "final") %>% 
  left_join(df_FS_tidy %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(rand_diff_final = savings - FS) %>% 
  select(name, site, rand_diff_final)

p1 <- ggplot() +
  geom_boxplot(data = conv_T, aes(x = "Conventional", y = conv_diff, fill = "Conventional M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_eob_T, aes(x = "Randomized", y = rand_diff_eob, fill = "Randomized M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_final_T, aes(x = "Randomized\n(24 months)", y = rand_diff_final, fill = "Randomized M&V (24 months)"), outlier.shape = NA) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = "Tidy set") +
  coord_cartesian(ylim = c(-8, 8)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the messy subset
rand_eob_M <- df_seq_FS_messy %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_messy %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(rand_diff_eob = savings - FS) %>% 
  select(name, site, rand_diff_eob)

conv_M <- df_FS_messy %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(conv_diff = true - conv) %>% 
  select(name, site, conv_diff)

rand_final_M <- df_seq_FS_messy %>% 
  filter(seq == "final") %>% 
  left_join(df_FS_messy %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(rand_diff_final = savings - FS) %>% 
  select(name, site, rand_diff_final)

p2 <- ggplot() +
  geom_boxplot(data = conv_M, aes(x = "Conventional", y = conv_diff, fill = "Conventional M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_eob_M, aes(x = "Randomized", y = rand_diff_eob, fill = "Randomized M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_final_M, aes(x = "Randomized\n(24 months)", y = rand_diff_final, fill = "Randomized M&V (24 months)"), outlier.shape = NA) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  labs(fill = NULL, 
       x = NULL, 
       y = NULL, 
       subtitle = "Messy set") +
  coord_cartesian(ylim = c(-8, 8)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_M, rand_eob_T)
conv_A <- bind_rows(conv_M, conv_T)
rand_final_A <- bind_rows(rand_final_M, rand_final_T)

p3 <- ggplot() +
  geom_boxplot(data = conv_A, aes(x = "Conventional", y = conv_diff, fill = "Conventional M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_eob_A, aes(x = "Randomized", y = rand_diff_eob, fill = "Randomized M&V"), outlier.shape = NA) +
  geom_boxplot(data = rand_final_A, aes(x = "Randomized\n(24 months)", y = rand_diff_final, fill = "Randomized M&V (24 months)"), outlier.shape = NA) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 3), 
                     labels = number_format(suffix = " %")) +
  scale_fill_manual(values = ls_colors) +
  labs(fill = NULL, 
       x = NULL, 
       y = "Difference in fractional savings", 
       subtitle = "All buildings") +
  coord_cartesian(ylim = c(-8, 8)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p3, p1, p2,
          ncol = 3, nrow = 1,
          labels = c("a)", "b)", "c)"),
          align = "hv",
          legend="none") +
  plot_annotation(title = "Savings estimation accuracy comparison", 
                  subtitle = "with measured weather conditions")

