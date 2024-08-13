
cont_plot <- function(df_means, sprt_res, annual_saving, true_saving, eob){
  
  #' @description This is a function defined to plot the overal results from sequential analysis. 
  #' @usage seq_plot(df_means, sprt_res, annual_saving, true_saving, eob)
  #' @param df_means running weekly mean
  #' @param sprt_res sequential probability ratio test results
  #' @param annual_saving annual saving estimated from towt model prediction weekly
  #' @param true_saving true saving calculated from original dataset
  #' @param eob n_weeks at the end of blocking period (end of the test)
  #' @return a stacked ggplots of different results
  #' 
  
  p2 <- ggplot(sprt_res, aes(x = n_weeks, y = ns_stat)) +
    geom_ribbon(aes(ymin = ns_ci_low, ymax = ns_ci_high), alpha = 0.5, fill = "#D0E3F1") +
    geom_line(size = 1.25, color = "grey20", alpha = 0.4) +
    geom_line(data = true_saving, aes(x = n_weeks, y = savings, color = "True savings"), 
              linetype = "dashed") +
    geom_text(data = sprt_res[5, ], 
              aes(x = n_weeks, y = ns_ci_high, label = "Upper 90% CI"),
              position = position_nudge(x = 1), 
              color = "grey20", 
              size = 2.5, 
              fontface = "italic") +
    geom_text(data = sprt_res[5, ], 
              aes(x = n_weeks, y = ns_ci_low, label = "Lower 90% CI"),
              position = position_nudge(x = 1), 
              color = "grey20", 
              size = 2.5, 
              fontface = "italic") +
    geom_point(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
               size = 8, 
               shape = 16, 
               color = "#347EB3", 
               alpha = 0.5) +
    geom_point(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
               size = 4, 
               shape = 16, 
               color = "#347EB3", 
               alpha = 0.9) +
    geom_text(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
              aes(x = n_weeks, y = ns_stat, label = paste0(round(ns_stat, 1L), " kW")), 
              position = position_nudge(x = 2.5), 
              color = "grey20", 
              size = 3.0, 
              fontface = "italic") +
    geom_errorbar(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
                  aes(x = n_weeks, ymin = ns_ci_low, ymax = ns_ci_high), 
                  width = 1, 
                  color = "#347EB3", 
                  alpha = 0.5, 
                  size = 0.8) +
    geom_text(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
              aes(x = n_weeks, y = ns_ci_high, label = paste0(round(ns_ci_high, 1L), " kW")), 
              position = position_nudge(x = 2.5), 
              color = "grey20", 
              size = 2.5, 
              fontface = "italic") +
    geom_text(data = sprt_res %>% filter(n_weeks%%block_params$block_unit == 0), 
              aes(x = n_weeks, y = ns_ci_low, label = paste0(round(ns_ci_low, 1L), " kW")), 
              position = position_nudge(x = 2.5), 
              color = "grey20", 
              size = 2.5, 
              fontface = "italic") +
    geom_text(data = slice_tail(sprt_res, n = 1), 
              aes(x = n_weeks, y = ns_stat, label = paste0(round(ns_stat, 1L), " kW")), 
              position = position_nudge(x = 0.5), 
              color = "grey20", 
              size = 3.0, 
              check_overlap = TRUE, 
              hjust = 0) +
    geom_text(data = slice_tail(sprt_res, n = 1), 
              aes(x = ifelse(is.na(first(flag)), 0, n_weeks), 
                  y = ns_ci_high, 
                  label = ifelse(is.na(first(flag)), NULL, paste0(round(ns_ci_high, 1L), " kW"))), 
              position = position_nudge(x = 0.5), 
              color = "grey20", 
              size = 2.5, 
              check_overlap = TRUE, 
              hjust = 0) +
    geom_text(data = slice_tail(sprt_res, n = 1), 
              aes(x = ifelse(is.na(first(flag)), 0, n_weeks), 
                  y = ns_ci_low, 
                  label = ifelse(is.na(first(flag)), NULL, paste0(round(ns_ci_low, 1L), " kW"))), 
              position = position_nudge(x = 0.5), 
              color = "grey20", 
              size = 2.5, 
              check_overlap = TRUE, 
              hjust = 0) +
    geom_vline(xintercept = seq(0, sprt_param$n_weeks + 1, by = block_params$block_unit), 
               linetype = "dashed", 
               color = "grey20", 
               alpha = 0.3, 
               size = 0.5) +
    annotate(geom = "text", 
             x = seq(5, sprt_param$n_weeks + 1, by = block_params$block_unit), 
             y = sprt_res[2, ] %>% .$ns_ci_high, 
             label = paste0(str_glue("{block_params$block_unit}-week block")), 
             alpha = 0.5, 
             size = 3) +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(eob - 0.5, cont_param$n_weeks + 0.5),
                       labels = number_format(accuracy = 1L, suffix = "\nweeks")) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = pretty_breaks(n = 3),
                       labels = number_format(suffix = " kW")) +
    labs(title = NULL, 
         subtitle = "SPRT results and estimated difference in power consumption without weather normalization",
         x = NULL, 
         y = NULL, 
         color = NULL) +
    scale_color_manual(values = ls_colors) +
    guides(alpha = "none") +
    coord_cartesian(clip = "off") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          plot.margin = margin(3, 15, 3, 3, unit = "mm"))
  
  
  # do running mean plot
  p1 <- ggplot(df_means, aes(x = week, y = value_ave)) +
    geom_line(aes(color = strategy), linewidth = 1.25) +
    geom_text(data = filter(df_means, 
                            week == cont_param$n_weeks, 
                            parameter == cont_param$parameter), 
              aes(x = cont_param$n_weeks, 
                  y = value_ave, 
                  color = strategy, 
                  label = strategy), 
              position = position_nudge(x = 0.5), 
              size = 3.0, 
              hjust = 0) +
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(eob - 0.5, cont_param$n_weeks + 0.5)) +
    scale_y_continuous(expand = c(0, 0), 
                       breaks = breaks_pretty(n = 4),
                       labels = number_format(suffix = " kW")) +
    scale_color_manual(values = ls_colors) +
    labs(title = NULL, 
         subtitle = "Running average power consumption of the case study building",
         x = NULL, 
         y = NULL) +
    guides(alpha = "none", color = "none") +
    coord_cartesian(clip = "off") +
    theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.25),
          axis.text.x = element_blank(),
          plot.margin = margin(3, 15, 3, 3, unit = "mm"))
  
  # Annotation
  # ES: Estimation start
  # BTC: Baseline temperature check
  # ITC: Interventional temperature check
  # MinC: Minimum test check
  # SPRT: SPRT check
  # EE: End of the block period
  
  p3 <- ggplot(annual_saving) +
    geom_line(aes(x = n_weeks, y = annual_saving), 
              linewidth = 1.25, color="darkgrey") +
    geom_point(data = annual_saving %>% filter(n_weeks == eob), 
               aes(x = n_weeks, y = annual_saving), 
               size = 4,
               alpha = 0.8, 
               shape = 16, 
               color = "#fc4e2a") +
    geom_text(data = annual_saving %>% filter(n_weeks == eob),
              aes(x = n_weeks, y = annual_saving, label = paste0("EE: ", n_weeks)),
              position = position_nudge(y = -1),
              color = "grey20", 
              size = 3, 
              fontface = "bold") +
    geom_vline(xintercept = seq(0, sprt_param$n_weeks + 1, by = block_params$block_unit), 
               linetype = "dashed", 
               color = "grey20", 
               alpha = 0.3, 
               size = 0.5) +
    annotate(geom = "text", 
             x = seq(5, cont_param$n_weeks + 1, by = block_params$block_unit), 
             y = max(annual_saving$annual_saving), 
             vjust = -2, 
             label = paste0(str_glue("{block_params$block_unit}-week block")), 
             alpha = 0.5, 
             size = 3) +
    geom_text(data = annual_saving[nrow(annual_saving), ], 
              aes(x = n_weeks, y = annual_saving, label = paste0(round(annual_saving, 0), "%")), 
              position = position_nudge(x = 0.5), 
              color = "grey20", 
              size = 3.0, 
              check_overlap = TRUE, 
              hjust = 0) +
    scale_x_continuous(expand = c(0, 0), 
                       breaks = seq(0, cont_param$n_weeks, by = block_params$block_unit), 
                       limits = c(eob - 0.4, cont_param$n_weeks + 0.5),
                       labels = number_format(accuracy = 1L, suffix = "\nweeks")) +
    scale_y_continuous(breaks = breaks_pretty(n = 3),
                       labels = number_format(suffix = "%")) +
    coord_cartesian(clip = "off") +
    labs(title = NULL, 
         subtitle = "Normalized annual energy savings estimate using TOWT model predictions for a typical meteorological year",
         x = NULL, 
         y = NULL, 
         color = NULL)  +
    theme(panel.grid.major.y = element_line(),
          legend.position = "bottom",
          plot.margin = margin(3, 15, 3, 3, unit = "mm"))
  
  # Overall result
  p1 / p2 / p3 +
    plot_annotation(title = str_glue("Overall sequential evaluation results of {name} at {site}")) +
    plot_annotation(tag_levels = c('a'), tag_suffix = ')') &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(color="black"))
}