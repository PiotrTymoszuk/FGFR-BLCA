# Plots of colflunce in time of cell cultures treated with Erdafitinib

  insert_head()

# container ------

  iv_plots <- list()

# plotting -------

  insert_msg("Plotting")

  iv_plots$line_plot <- iv$data %>%
    mutate(html_label = reorder(html_label, as.integer(genotype))) %>%
    ggplot(aes(x = concentration,
               y = confluence,
               color = FGFR3_mutation)) +
    facet_wrap(facets = "html_label",
               scales = "free") +
    geom_line(aes(group = condition)) +
    geom_point(shape = 16,
               size = 2) +
    scale_color_manual(values = c(WT = "gray30",
                                  mutated = globals$mut_status_color[["mutated"]]),
                       name = html_italic("FGFR3")) +
    scale_x_continuous(breaks = unique(iv$data$concentration)) +
    scale_y_continuous(limits = c(0, 100)) +
    globals$common_theme +
    theme(strip.text.x = element_markdown(),
          legend.title = element_markdown()) +
    labs(title = "In vitro erdafitinib treatment",
         x = "erdafitinib, µM",
         y = "confluence, %")

# END ------

  insert_tail()
