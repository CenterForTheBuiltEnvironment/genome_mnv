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
run_params <- list(type = "messy")

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




#### PLOT ####
# Normalized savings
p1 <- df_sprt_all %>% 
  filter(seq == "eob") %>% 
  left_join(df_FS %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
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
  left_join(df_FS %>% filter(scenario == "ref" & method == "true"), by = c("name", "site")) %>% 
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
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, 0.95)

sequence <- seq(from = 5, to = max(plot_data$abs_diff), by = 10)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
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
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
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
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 50, to = max(plot_data$abs_diff), by = 50)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
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
                x = 150, 
                label = paste0("(", "50% buildings: < ", round(median, digits = 1), " kW)")), 
            position = position_nudge(y = 5), 
            check_overlap = T,
            size = 4) +
  geom_hline(yintercept = upper, 
             linetype = "dashed", 
             color = "red") +
  geom_text(aes(y = upper, 
                x = 150, 
                label = paste0("(", "95% buildings: < ", round(upper, digits = 1), " kW)")), 
            position = position_nudge(y = 5), 
            check_overlap = T,
            size = 4) +
  scale_x_continuous(expand = c(0.02, 0), 
                     breaks = breaks_pretty(n = 3)) +
  scale_y_continuous(expand = c(0, 0), 
                     breaks = sequence,
                     labels = number_format(suffix = " kW"),
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
          align = "hv",
          legend="bottom") +
  plot_annotation(title = str_glue("Deviation in saving estimation after sequential test\n(Measured savings from observed weather)"))

ggsave(filename = str_glue("md_comp_cont.png"), path = combifigs_path, units = "in", height = 8, width = 8, dpi = 300)

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
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 5, to = max(plot_data$abs_diff), by = 10)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
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
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
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
  geom_boxplot(aes(y = abs_diff)) +
  scale_y_continuous(expand = c(0, 0), 
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = NULL, 
       y = NULL, 
       subtitle = "Aggregated") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

median <- median(plot_data$abs_diff)
upper <- quantile(plot_data$abs_diff, probs = 0.95)

sequence <- seq(from = 5, to = max(plot_data$abs_diff), by = 15)
seq_df <- data.frame(knot = sequence)
seq_data <- seq_df %>%
  rowwise() %>%
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
                     limits = c(0, max(plot_data$plot_max) + 2)) +
  labs(x = "Number of buildings", 
       y = "Absolute difference in measured savings", 
       subtitle = "Accumulated breakdown by deviation") +
  coord_cartesian(ylim = c(0, max(plot_data$plot_max) + 2)) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))

ggarrange(p2, p1,
          ncol=2, nrow=1,
          labels = c("a)", "b)"),
          widths = c(4, 1),
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
  geom_col(aes(x = name, y = diff_in_diff), position = "identity", alpha = 0.5) +
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

