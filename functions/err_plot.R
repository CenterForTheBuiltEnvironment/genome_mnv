err_plot <- function(df_rand, df_base, df_interv){
  
  rand_mean <- list()
  base_mean <- list()
  interv_mean <- list()
  
  # calculate weekly means
  for (i in 2:sprt_param$n_weeks) {
    
    # subset data by week
    rand_mean[[i]] <- df_rand %>%
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week == i) %>%
      group_by(strategy) %>%
      summarise(week = i,
                rand_ave = mean(eload, na.rm = TRUE)) %>%
      ungroup()
    
    base_mean[[i]] <- df_base %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week == i) %>%
      summarise(week = i,
                base_ave = mean(eload, na.rm = TRUE)) %>%
      ungroup()
    
    interv_mean[[i]] <- df_interv %>% 
      mutate(week = interval(min(datetime), datetime) %>% as.numeric('weeks') %>% floor()) %>% 
      filter(week == i) %>%
      summarise(week = i,
                interv_ave = mean(eload, na.rm = TRUE)) %>%
      ungroup()
  }
  
  rand_mean <- bind_rows(rand_mean) %>% 
    mutate(strategy = as.factor(strategy), 
           strategy = recode(strategy, "1" = "Baseline", "2" = "Intervention"))
  
  base_mean <- bind_rows(base_mean)
  interv_mean <- bind_rows(interv_mean)
  
  p1 <- ggplot() +
    geom_line(data = rand_mean %>% 
                filter(strategy == "Baseline"), 
              aes(x = week, y = rand_ave, color = strategy)) +
    geom_line(data = base_mean, aes(x = week, y = base_ave, color = "True")) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    labs(x = NULL, 
         y = NULL, 
         color = NULL, 
         subtitle = "Baseline mean estimation") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.x = element_blank(),
          legend.position = "none",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  
  p2 <- ggplot() +
    geom_line(data = rand_mean %>% 
                filter(strategy == "Intervention"), 
              aes(x = week, y = rand_ave, color = strategy)) +
    geom_line(data = interv_mean, aes(x = week, y = interv_ave, color = "True")) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    labs(x = NULL, 
         y = NULL, 
         color = NULL, 
         subtitle = "Intervention mean estimation") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.x = element_blank(), 
          legend.position = "none",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  base_error <- rand_mean %>% 
    filter(strategy == "Baseline") %>% 
    left_join(base_mean, by = c("week")) %>% 
    mutate(base_diff = rand_ave - base_ave) %>% 
    select(week, base_diff)
  
  interv_error <- rand_mean %>% 
    filter(strategy == "Intervention") %>% 
    left_join(interv_mean, by = c("week")) %>% 
    mutate(interv_diff = rand_ave - interv_ave) %>% 
    select(week, interv_diff)
  
  saving_error <- base_error %>% 
    left_join(interv_error, by = "week") %>% 
    mutate(saving_error = interv_diff - base_diff) %>% 
    mutate(acc_error = cumsum(ifelse(is.na(saving_error), 0, saving_error)))
  
  
  p3 <- saving_error %>% 
    ggplot() +
    geom_line(aes(x = week, y = saving_error)) +
    geom_hline(yintercept = 0, color = "red", lty = "dashed") +
    annotate(geom = "text", 
             x = 3,
             y = 2, 
             size = 3, 
             color = "grey50", 
             label = "Underestimate") +
    annotate(geom = "text", 
             x = 3,
             y = -2, 
             size = 3, 
             color = "grey50", 
             label = "Overestimate") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    labs(x = NULL, 
         y = NULL, 
         color = NULL, 
         subtitle = "mean saving error estimation") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.x = element_blank(), 
          legend.position = "none",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  p4 <- saving_error %>% 
    ggplot() +
    geom_line(aes(x = week, y = acc_error)) +
    geom_hline(yintercept = 0, color = "red", lty = "dashed") +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kWh")) +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(0, sprt_param$n_weeks + 0.5),
                       breaks = seq(0, sprt_param$n_weeks, by = block_params$block_unit), 
                       labels = number_format(suffix = "\nweeks")) +
    labs(x = NULL, 
         y = NULL, 
         color = NULL, 
         subtitle = "accumulated saving error estimation") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.position = "none",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  ggarrange(p1, p2, p3, p4,
            ncol=1, nrow=4,
            labels = c("a)", "b)", "c)", "d)"),
            align = "hv") +
    plot_annotation(title = str_glue("Deviation in sequnetial saving estimation"))
  
}