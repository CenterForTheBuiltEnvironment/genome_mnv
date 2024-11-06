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
               "Randomized" = "#66c2a4", 
               "Randomized\n(24 months)" = "#238b45",
               "Buildings (tidy) finishing randomized M&V" = "#bcbddc",
               "Buildings (messy) finishing randomized M&V" = "#756bb1",
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
fig_path = "../figs/graphical_abs/"

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
  mutate(diff = savings - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_T <- df_FS_tidy %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_T <- df_FS_tidy %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_acc_T <- rbind(rand_eob_T, conv_T, rand_final_T)

p1 <- df_acc_T %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_lv(alpha = 0.4, k = 4, outlier.size = 0.5) +
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
       subtitle = str_glue("Tidy set ({nrow(df_acc_T)/3} buildings)")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

# plot for the messy subset
rand_eob_M <- df_seq_FS_messy %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS_messy %>% filter(scenario == "ref" & method == "true"), by = c("name", "site"))%>% 
  mutate(diff = savings - FS, 
         method = "rand_eob") %>% 
  select(name, method, diff)

conv_M <- df_FS_messy %>% 
  filter(method != "rand") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - conv, 
         method = "conv") %>% 
  select(name, method, diff)

rand_final_M <- df_FS_messy %>% 
  filter(method != "conv") %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  mutate(diff = true - rand, 
         method = "rand_final") %>% 
  select(name, method, diff)

df_acc_M <- rbind(rand_eob_M, conv_M, rand_final_M)

p2 <- df_acc_M %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv" = "Conventional", "rand_eob" = "Randomized", "rand_final" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_lv(alpha = 0.4, k = 4, outlier.size = 0.5) +
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
       subtitle = str_glue("Messy set ({nrow(df_acc_M)/3} buildings)")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        axis.text.y = element_blank(), 
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


# plot for combined dataset
rand_eob_A <- bind_rows(rand_eob_M, rand_eob_T) %>% 
  select(name, rand_eob_diff = diff)

conv_A <- bind_rows(conv_M, conv_T) %>% 
  select(name, conv_diff = diff)

rand_final_A <- bind_rows(rand_final_M, rand_final_T) %>% 
  select(name, rand_final_diff = diff)

df_acc_A <- rand_eob_A %>% 
  left_join(conv_A, by = "name") %>% 
  left_join(rand_final_A, by = "name") %>% 
  pivot_longer(c(rand_eob_diff, conv_diff, rand_final_diff), names_to = "method", values_to = "diff")

p3 <- df_acc_A %>% 
  mutate(method = as.factor(method), 
         method = recode_factor(method, "conv_diff" = "Conventional", "rand_eob_diff" = "Randomized", "rand_final_diff" = "Randomized\n(24 months)")) %>% 
  ggplot(aes(x = method, y = diff, fill = method)) +
  geom_lv(alpha = 0.4, k = 4, outlier.size = 0.5) +
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
       subtitle = str_glue("All ({nrow(df_acc_A)/3} buildings)")) +
  coord_cartesian(ylim = c(-18, 18)) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.position = "none",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

p_top <- ggarrange(p3, p1, p2,
          ncol = 3, nrow = 1,
          labels = c("a)", "b)", "c)"),
          align = "hv",
          legend="none") +
  plot_annotation(title = "M&V accuracy comparison", 
                  subtitle = "with measured weather conditions")

# Timeline plot
df_time <- df_sprt_all_messy %>% 
  filter(seq != "final") %>% 
  mutate(set = 'messy') %>% 
  bind_rows(df_sprt_all_tidy %>% filter(seq != "final") %>% mutate(set = 'tidy')) %>% 
  select(name, seq, n_weeks, set)

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
  
  n_tidy <- df_temp %>% 
    filter(set == "tidy")
  
  n_messy <- df_temp %>% 
    filter(set == "messy")
  
  count[[n]] <- tibble("n_weeks" = i, 
                       "eob" = nrow(df_eob), 
                       "temp" = nrow(df_temp), 
                       "sprt" = nrow(df_sprt), 
                       "tidy" = nrow(n_tidy), 
                       "messy" = nrow(n_messy))
  
  n <- n + 1
}

count <- bind_rows(count) 

p_bottom <- count %>% 
  ggplot() +
  geom_bar(data = . %>% 
             filter(n_weeks == 24 | n_weeks == 36) %>% 
             pivot_longer(c(tidy, messy), names_to = "set", values_to = "value") %>% 
             mutate(set = as.factor(set), 
                    set = recode_factor(set, "tidy" = "Buildings (tidy) finishing randomized M&V", "messy" = "Buildings (messy) finishing randomized M&V")), 
           aes(x = n_weeks, y = value, fill = set), 
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
               arrow = arrow(length = unit(0.3, "in")),   
               linewidth = 2,  
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
       title = "M&V timeline comparison") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


ggarrange(p_top, p_bottom, 
          ncol = 1, nrow = 2,
          heights = c(1.5, 1),
          labels = c(" ", "d)"),
          legend = "bottom")

ggsave(filename = "graphical_abs.png", path = fig_path, units = "in", height = 9, width = 12, dpi = 300)
