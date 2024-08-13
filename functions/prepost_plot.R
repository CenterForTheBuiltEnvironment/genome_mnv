
prepost_plot <- function(base_pre_meas,  base_pos_proj, base_pos_true, interv_pos_meas, title){
  
  p1 <- ggplot() +
    geom_point(data = base_pre_meas %>%
                 slice_sample(n = 500),
               aes(x = datetime, y = eload, color = "Measured baseline"),
               size = 0.5,
               alpha = 0.3) +
    geom_smooth(data = base_pre_meas,
                aes(x = datetime, y = eload, color = "Measured baseline"), 
                formula = y ~ x, method = "loess") +
    scale_x_datetime(date_breaks = "2 months",
                     date_labels = "%b")  +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = plot_scale) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         title = NULL,
         subtitle = "Pre-retrofit period") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.direction = "horizontal",
          legend.position = "bottom",
          plot.margin = margin(t = 2, r = 7, b = 2, l = 2, unit = "mm"))
  
  p2 <- ggplot() +
    geom_point(data = base_pos_proj %>%
                 slice_sample(n = 500),
               aes(x = datetime, y = towt, color = "Projected baseline\n(no change)"),
               size = 0.5,
               alpha = 0.3) +
    geom_smooth(data = base_pos_proj,
                aes(x = datetime, y = towt, color = "Projected baseline\n(no change)"), 
                formula = y ~ x, method = "loess") +
    geom_point(data = base_pos_true %>%
                 slice_sample(n = 500),
               aes(x = datetime, y = eload, color = "Adjusted baseline"),
               size = 0.5,
               alpha = 0.3) +
    geom_smooth(data = base_pos_true,
                aes(x = datetime, y = eload, color = "Adjusted baseline"), 
                formula = y ~ x, method = "loess") +
    geom_point(data = interv_pos_meas %>%
                 slice_sample(n = 500),
               aes(x = datetime, y = eload, color = "Measured interv"),
               size = 0.5,
               alpha = 0.3) +
    geom_smooth(data = interv_pos_meas,
                aes(x = datetime, y = eload, color = "Measured interv"), 
                formula = y ~ x, method = "loess") +
    scale_x_datetime(date_breaks = "2 months",
                     date_labels = "%b")  +
    scale_y_continuous(expand = c(0, 0),
                       breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    coord_cartesian(ylim = plot_scale) +
    labs(x = NULL,
         y = NULL,
         color = NULL,
         title = NULL,
         subtitle = "post-retrofit period") +
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
            legend="bottom") +
    plot_annotation(title = title)
}