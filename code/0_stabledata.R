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
theme_update(plot.title = element_text(size = 14, colour = "grey20", face = "bold", hjust = 0.5),
             plot.subtitle = element_text(size = 10, colour = "grey20", face = "italic", hjust = 0.5, margin = margin(b = 10)),
             plot.caption = element_text(size = 8, colour = "grey20", face = "italic", hjust = 0.5),
             plot.background = element_rect(fill = "white", colour = NA),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             strip.text = element_text(size = 10, color = "grey20", face = "bold"),
             strip.background = element_blank())

# T test
perform_t_test <- function(df_pre, df_post, building) {
  t.test(df_pre[[building]], df_post[[building]], alternative = "two.sided")
}




#### READ DATA ####
data_path <- "../../buds-lab-building-data-genome-project-2/data/meters/cleaned/"
meta_path <- "../../buds-lab-building-data-genome-project-2/data/metadata/"
weather_path <- "../../buds-lab-building-data-genome-project-2/data/weather/"
output_path <- "../readfiles/stable/"

df_elec <- read_csv(paste0(data_path, "electricity_cleaned.csv"))
df_meta <- read_csv(paste0(meta_path, "metadata.csv"))
df_weather <- read_csv(paste0(weather_path, "weather.csv"))




#### FILTER ####
# Focus
remove_type <- c("parking", "warehouse", "utility")
remove_site <- c("Eagle", "Bobcat", "Swan", "Hog")

# NAs
na_counts <- sapply(df_elec[-1], function(x) sum(is.na(x)))
cols_to_keep <- names(na_counts[na_counts <= 1000])
cols_to_keep <- c("timestamp", cols_to_keep)

df_elec <- df_elec[, cols_to_keep]

# 0 means
means <- sapply(df_elec[-1], mean, na.rm = TRUE)
cols_to_keep <- names(means[means > 1])
cols_to_keep <- c("timestamp", cols_to_keep)

df_elec <- df_elec[, cols_to_keep]

# yearly difference
df_elec$timestamp = as.POSIXct(df_elec$timestamp, format="%Y-%m-%d %H:%M:%OS")

df_pre <- df_elec %>%
  filter(year(timestamp) == 2016)

df_post <- df_elec %>%
  filter(year(timestamp) == 2017)

buildings <- names(df_elec)[names(df_elec) != "timestamp"]

p_values <- sapply(buildings, function(building) {
  t_test_result <- perform_t_test(df_pre, df_post, building)
  t_test_result$p.value
})

results <- data.frame(Building = buildings, P_Value = p_values) %>% 
  filter(P_Value >= 0.05)

cols_to_keep <- c("timestamp", results %>% .$Building)
df_elec <- df_elec[, cols_to_keep]

# process corresponding metadata
df_meta <- df_meta %>% 
  filter(building_id %in% results$Building) %>% 
  select(building_id, sqm, sqft, eui) %>% 
  separate(building_id, into = c("site", "type", "name"), sep = "_") %>% 
  filter(!type %in% remove_type) %>% 
  filter(!site %in% remove_site)

# process corresponding weather data
df_weather <- df_weather %>% 
  select(timestamp, 
         site = site_id, 
         t_out = airTemperature) %>% 
  mutate(timestamp = as.POSIXct(timestamp, format="%Y-%m-%d %H:%M:%OS")) %>% 
  filter(site %in% df_meta$site)





#### OUTPUT ####
df_energy <- df_elec %>% 
  pivot_longer(c(-timestamp), names_to = "buildings", values_to = "eload") %>% 
  separate(buildings, into = c("site", "type", "name"), sep = "_") %>% 
  filter(!type %in% remove_type) %>% 
  filter(!site %in% remove_site)

write_rds(df_energy, paste0(output_path, "df_energy.rds"), compress = "gz")
write_rds(df_meta, paste0(output_path, "df_meta.rds"), compress = "gz")
write_rds(df_weather, paste0(output_path, "df_weather.rds"), compress = "gz")
