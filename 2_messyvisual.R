# Genome 2 M&V analysis
# written by Aoyu



#### SETUP ####

# use pacman
require(pacman)

# load libraries
pacman::p_load(tidyverse, lubridate, scales, slider, cowplot, patchwork, RColorBrewer, # general
               broom, ggpmisc, ggpubr) # plotting

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

# define parameters
run_params <- list(type = "messy", 
                   site = T)

sprt_param <- list(baseline = "Baseline",
                   strategy = "Intervention",
                   parameter = "power_ave",
                   label = "power",
                   n_weeks = 48)

block_params <- list(start_date = "2016-01-01",
                     n_weeks = 108,
                     n_seasons = 9, 
                     block_unit = 12)



#### FUNCTIONS ####
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

#### READ DATA ####
readfile_path <- str_glue("./readfiles/{run_params$type}/")
summaryfigs_path <- str_glue("./figs/{run_params$type}/site_summary/")
combifigs_path <- str_glue("./figs/{run_params$type}/comb_analysis/")

df_energy <- read_rds(paste0(readfile_path, "df_energy.rds"))
df_meta <- read_rds(paste0(readfile_path, "df_meta.rds"))
df_weather <- read_rds(paste0(readfile_path, "df_weather.rds"))
df_sprt_all <- read_rds(paste0(readfile_path, "df_sprt_all.rds"))
df_seq_FS <- read_rds(paste0(readfile_path, "df_seq_FS.rds"))

if (run_params$type == "tidy"){
  
  df_NRE_occ <- read_rds(paste0(readfile_path, "df_NRE_occ.rds"))
  
}

df_MD <- read_rds(paste0(readfile_path, "df_MD.rds"))
df_FS <- read_rds(paste0(readfile_path, "df_FS.rds"))
df_eui <- read_rds(paste0(readfile_path, "df_eui.rds"))
df_cont_MD <- read_rds(paste0(readfile_path, "df_cont_MD.rds"))
df_cont_FS <- read_rds(paste0(readfile_path, "df_cont_FS.rds"))

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




#### SITE ####
if (run_params$site){
  set3 <- colorRampPalette(brewer.pal('Set3',n=12))
  type_colors <- setNames(set3(13), all_types$type)
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
    geom_text(aes(label = ifelse(n > 10, as.character(n), "")), color = "black", position = position_stack(vjust = 0.5)) +
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
    facet_wrap(~ site, nrow = 3) +
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
  
  ggsave(filename = "site_summary.png", path = summaryfigs_path, units = "in", height = 8, width = 8, dpi = 300)
  
  # EUI summary across all buildings
  eui <- df_energy %>%
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
  
  
  bar_plot <- eui %>%
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
  
  hist_plot <- eui %>%
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
    facet_wrap(~site, nrow = 3) +
    coord_cartesian(ylim = c(-8, 38)) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         title = "Outdoor weather conditions across all sites") +
    theme(panel.grid.major.y = element_line(color = "grey80"),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggsave(filename = "weather_summary.png", path = summaryfigs_path, units = "in", height = 8, width = 10, dpi = 300)
  
  # density plot of tmy and actual tout
  df_weather %>%
    filter(year(timestamp) == 2016) %>%
    mutate(type = "observed") %>%
    ggplot() +
    geom_histogram(aes(x = t_out, fill = type), alpha = 0.6, color = NA) +
    geom_histogram(data = df_tmy %>% mutate(type = "tmy"), aes(x = temp, fill = type), alpha = 0.6, color = NA) +
    scale_fill_manual(values = c("observed" = "#5941A9", "tmy" = "#E5D4ED")) +
    facet_wrap(~site, nrow = 3) +
    labs(x = NULL, y = NULL, color = NULL, fill = NULL,
         title = "Outdoor temperature histogram\nfor observed and typical conditions") +
    scale_x_continuous(labels = number_format(suffix = " °C")) +
    coord_cartesian(clip = "off") +
    theme(axis.text.y = element_blank(),
          legend.direction = "horizontal",
          legend.position = "bottom")
  
  ggsave(filename = "weather_tmy.png", path = summaryfigs_path, units = "in", height = 6, width = 10, dpi = 300)
  
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
  
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 3)
  combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.2))
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Daily average power usage by average outdoor weather conditions", fontface = 'bold', size = 18),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  print(final_plot)
  
  ggsave(filename = "energy_weather.png", path = summaryfigs_path, units = "in", height = 8, width = 12, dpi = 300)
  
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
  
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 3)
  combined_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.2))
  
  final_plot <- plot_grid(
    ggdraw() + draw_label("Average power usage by each hour of the day", fontface = 'bold', size = 18),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
  
  print(final_plot)
  
  ggsave(filename = "energy_hour_site.png", path = summaryfigs_path, units = "in", height = 8, width = 12, dpi = 300)
}




#### INDIVIDUAL ####
# General accuracy checking for both method
df_FS %>% 
  pivot_wider(names_from = method, values_from = savings) %>% 
  left_join(df_seq_FS %>% filter(seq == "eob"), by = c("name", "site")) %>% 
  pivot_longer(c(true, conv, rand, FS), names_to = "parameter", values_to = "value") %>% 
  mutate(parameter = as.factor(parameter), 
         parameter = recode_factor(parameter, 
                                   "true" = "True savings", 
                                   "conv" = "Conventional M&V", 
                                   "rand" = "Randomized M&V", 
                                   "FS" = "Randomized M&V at early stop")) %>% 
  ggplot(aes(x = parameter, y = value, fill = parameter)) +
  stat_boxplot(geom ='errorbar', width = 0.5) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = breaks_pretty(n = 5), 
                     limits = c(0, 16), 
                     labels = number_format(suffix = "%")) +
  labs(x = NULL, 
       y = NULL, 
       fill = NULL, 
       title = "Savings estimation comparison", 
       subtitle = "between convention and randomized M&V (with early stop)") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("FS_Dev.png"), path = combifigs_path, units = "in", height = 5, width = 8, dpi = 300)

# Sequential mean difference 
plot_data <- df_sprt_all %>% 
  filter(seq == "eob") %>% 
  left_join(df_MD %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - - sprt), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 10)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 1), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 120, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " kW)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 120, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " kW)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " kW"),
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in randomized M&V saving estimation at early stop"))

ggsave(filename = str_glue("md_comp_rand.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)

# conventional mean difference
plot_data <- df_MD %>% 
  filter(method != "rand") %>% 
  group_by(name) %>% 
  summarise(abs_diff = abs(sum(savings))) %>% 
  ungroup() %>% 
  mutate(plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 5)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 5)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 150)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 5), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 120, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " kW)")), 
            position = position_nudge(y = 5), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 120, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " kW)")), 
            position = position_nudge(y = 5), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " kW"),
                     limits = c(0, max(plot_data$plot_max) + 5)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 5)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in conventional M&V saving estimation"))

ggsave(filename = str_glue("md_comp_conv.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)

# continuous sprt mean difference 
plot_data <- df_cont_MD %>% 
  rename(cont = sprt) %>% 
  left_join(df_MD %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - - cont), 
         plot_max = max(abs_diff))


p1 <-  plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 25)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 1), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 120, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " kW)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 120, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " kW)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " kW"),
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in randomized M&V saving estimation after sequential test\n(continue with 20%/80% sampling)"))

ggsave(filename = str_glue("md_comp_cont.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)

# Sequential fractional savings
plot_data <- df_seq_FS %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - FS), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 0.5)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 0.5)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 5)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 0.5), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 150, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " %)")), 
            position = position_nudge(y = 0.5), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 150, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " %)")), 
            position = position_nudge(y = 0.5), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " %"),
                     limits = c(0, max(plot_data$plot_max) + 0.5)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 0.5)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in randomized M&V fractional saving estimation at early stop"))

ggsave(filename = str_glue("fr_comp_rand.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)

# conventional fractional savings
plot_data <- df_FS %>% 
  filter(method != "rand") %>% 
  group_by(name) %>% 
  summarise(abs_diff = abs(diff(savings))) %>% 
  ungroup() %>% 
  mutate(plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 0.5)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 0.5)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 10)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 0.5), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 120, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " %)")), 
            position = position_nudge(y = 0.5), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 120, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " %)")), 
            position = position_nudge(y = 0.5), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " %"),
                     limits = c(0, max(plot_data$plot_max) + 0.5)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 0.5)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in conventional M&V fractional saving estimation"))

ggsave(filename = str_glue("fr_comp_conv.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)

# continuous sprt fractional savings
plot_data <- df_cont_FS %>% 
  rename(cont = FS) %>% 
  left_join(df_FS %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
  mutate(abs_diff = abs(savings - cont), 
         plot_max = max(abs_diff))

p1 <- plot_data %>% 
  ggplot() +
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 0, to = max(plot_data$abs_diff), by = 25)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
  filter(knot > 0) %>% 
  mutate(cumulative_count = sum(plot_data$abs_diff <= knot)) %>%
  ungroup()

p2 <- seq_data %>% 
  ggplot() +
  geom_hline(aes(yintercept = knot), 
             linetype = "dashed", 
             color = "grey80") +
  geom_bar(data = plot_data %>% 
             arrange(abs_diff), 
           aes(x = seq(1, nrow(all_names), by = 1), y = abs_diff), 
           stat = "identity", alpha = 0.5) +
  geom_text(aes(x = cumulative_count, 
                y = knot, 
                label = paste0("(n = ", cumulative_count, ")")), 
            position = position_nudge(x = -25, y = 1), 
            size = 4) +
  geom_hline(yintercept = median(plot_data$abs_diff), 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = median, 
                x = 150, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " %)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 150, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " %)")), 
            position = position_nudge(y = 1), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " %"),
                     limits = c(0, max(plot_data$plot_max) + 1)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 1)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in randomized M&V fractional saving estimation after sequential test\n(continue with 20%/80% sampling)"))

ggsave(filename = str_glue("fr_comp_cont.png"), path = combifigs_path, units = "in", height = 7, width = 8, dpi = 300)


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
                     limits = c(-0.5, max(df_sprt_all %>% filter(seq != "final") %>% .$n_weeks) + 0.5),
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
        axis.text.y = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("seq_timeline.png"), path = combifigs_path, units = "in", height = 18, width = 10, dpi = 300)

# continue sampling visualization
# df_rand_cont %>% 
#   mutate(strategy = as.factor(strategy), 
#          strategy = recode_factor(strategy, "1" = "Baseline", "2" = "Intervention")) %>% 
#   mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor(), 
#          sample = as.factor(ifelse(week <= eob, "Before", "After")), 
#          sample = fct_relevel(sample, c("Before", "After"))) %>%
#   group_by(sample, strategy) %>% 
#   summarise(n = n()) %>% 
#   mutate(n = n, 
#          perc = n / sum(n)) %>%
#   ungroup() %>% 
#   ggplot(aes(x = sample, y = perc, fill = strategy)) +
#   geom_bar(position="fill", stat="identity") +
#   geom_hline(yintercept = 0.5, 
#              linetype = "dashed", 
#              color = "red") +
#   annotate(geom = "text",
#            color = "red", 
#            x = 0.5, 
#            y = 0.52, 
#            label = "50%") +
#   labs(x = NULL,
#        fill = NULL,
#        y = NULL,
#        title = "Sampling ratio between baseline and intervention", 
#        subtitle = "Before and after the sequential test") +
#   scale_y_continuous(expand = c(0, 0), 
#                      breaks = breaks_pretty(n = 4), 
#                      labels = number_format(suffix = "%", scale = 100)) +
#   scale_fill_manual(values = ls_colors) +
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom",
#         plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))


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
        axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("saving_cont.png"), path = combifigs_path, units = "in", height = 5, width = 15, dpi = 300)

# NRE saving estimation: occupancy change
s = "ref"
dev_conv <- df_FS %>%
  filter(scenario == s) %>%
  filter(method != "rand") %>%
  group_by(name, site, scenario) %>%
  summarise(conv = abs(diff(savings))) %>%
  ungroup()

dev_rand <- df_FS %>%
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
  geom_col(aes(x = name, y = diff_in_diff), position = "identity", alpha = 0.3) +
  facet_wrap(~site, nrow = 1, scales = "free_x") +
  scale_y_continuous(expand = c(0.1, 0),
                     breaks = breaks_pretty(n = 4), 
                     labels = number_format(suffix = "%")) +
  geom_text(data = mean_diff,
            aes(x = pos, y = -.5, group = site, label = paste0(mean, "%"))) +
  geom_line(aes(x = name, y = rand, color = "Absolute deviation of randomized method"), alpha = 0.5) +
  geom_point(aes(x = name, y = rand, color = "Absolute deviation of randomized method"), alpha = 0.5, size = 0.1) +
  labs(x = NULL,
       fill = NULL,
       y = NULL,
       color = NULL, 
       title = "Difference-in-difference of fractional savings calculated for each building") +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
        axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggsave(filename = str_glue("comp_savings.png"), path = combifigs_path, units = "in", height = 5, width = 15, dpi = 300)

