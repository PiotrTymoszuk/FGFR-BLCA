# Plots of SHAP variable importance

  insert_head()

# container ------

  mibc_splots <- list()

# SHAP and SHAPviz objects -------

  insert_msg('SHAP and SHAPviz objects')

  ## FGFR-related genes

  mibc_splots$fgfr_genes <- mibc_globals$lexicon %>%
    filter(protein_type != 'other') %>%
    .$variable

  ## class names and algorithm labels

  mibc_splots$class_levels <- mibc_globals$class_levels

  mibc_splots$algo_labs <-
    c('elnet' = 'Elastic Net',
      'ranger' = 'Random Forest')

  ## the SHAP/SHAPviz objects

  mibc_splots$shap_obj <- mibc_shap$shap_obj[c("elnet", "ranger")]

  mibc_splots$shapviz <- mibc_splots$shap_obj %>%
    map(shapviz) %>%
    map(set_names, mibc_splots$class_levels)

# Mean SHAP values for all and FGFR-related genes --------

  insert_msg('Mean SHAP values')

  mibc_splots$stats <- mibc_splots$shapviz %>%
    map(get_shap_values) %>%
    #map(map, abs) %>%
    map(map, colMeans) %>%
    map(map,
        compress,
        names_to = 'variable',
        values_to = 'importance') %>%
    map(compress, names_to = 'class') %>%
    map(mutate,
        fgfr_gene = ifelse(variable %in% mibc_splots$fgfr_genes,
                           'yes', 'no'),
        fgfr_label = ifelse(fgfr_gene == 'yes',
                            variable, NA),
        fgfr_color_var = ifelse(fgfr_gene == 'no', 'other',
                                ifelse(importance > 0,
                                       'positive', 'negative')),
        fgfr_color_var = factor(fgfr_color_var,
                                c('positive', 'negative', 'other')),
        class = factor(class, mibc_splots$class_levels))

# Bee-Swarms of SHAPs and bar plots of mean absolute SHAPS ------

  insert_msg('Plots')

  mibc_splots$plots <-
    list(shapviz_obj = mibc_splots$shapviz,
         title_suffix = mibc_splots$algo_labs) %>%
    pmap(my_shap_bee,
         panels = TRUE,
         add_legend = FALSE,
         top_variables = 15,
         point_size = 1)

# Violin plots for the mean SHAP values, FGFR-related genes highlighted -------

  insert_msg('Violin plots of the mean SHAP values')

  mibc_splots$mean_plots <-
    list(x = mibc_splots$stats,
         y = mibc_splots$algo_labs) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = importance,
                      y = factor(class, rev(levels(class))))) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_violin(fill = 'cornsilk',
                       alpha = 0.25) +
           geom_point(data = x %>%
                        filter(fgfr_gene == 'no'),
                      aes(color = fgfr_color_var,
                          size = fgfr_color_var,
                          shape = fgfr_color_var,
                          alpha = fgfr_color_var),
                      position = position_jitter(height = 0.05,
                                                 width = 0)) +
           geom_point(data = x %>%
                        filter(fgfr_gene == 'yes'),
                      aes(color = fgfr_color_var,
                          size = fgfr_color_var,
                          shape = fgfr_color_var,
                          alpha = fgfr_color_var)) +
           geom_text_repel(aes(label = fgfr_label,
                               color = fgfr_color_var),
                           size = 2.5,
                           fontface = 'italic',
                           box.padding = 0.6,
                           force = 2,
                           show.legend = FALSE) +
           scale_color_manual(values = c(other = 'gray40',
                                         positive = 'firebrick',
                                         negative = 'steelblue'),
                              labels = c(other = 'other',
                                         positive = 'FGFR/FGF/FGFBP\npositive predictor',
                                         negative = 'FGFR/FGF/FGFBP\nnegative predictor'),
                              name = 'gene classification') +
           scale_size_manual(values = c(other = 1,
                                        positive = 3,
                                        negative = 3),
                             labels = c(other = 'other',
                                        positive = 'FGFR/FGF/FGFBP\npositive predictor',
                                        negative = 'FGFR/FGF/FGFBP\nnegative predictor'),
                             name = 'gene classification') +
           scale_shape_manual(values = c(other = 16,
                                         positive = 18,
                                         negative = 18),
                              labels = c(other = 'other',
                                         positive = 'FGFR/FGF/FGFBP\npositive predictor',
                                         negative = 'FGFR/FGF/FGFBP\nnegative predictor'),
                              name = 'gene classification') +
           scale_alpha_manual(values = c(other = 0.25,
                                         positive = 1,
                                         negative = 1),
                              labels = c(other = 'other',
                                         positive = 'FGFR/FGF/FGFBP\npositive predictor',
                                         negative = 'FGFR/FGF/FGFBP\nnegative predictor'),
                              name = 'gene classification') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = paste('full predictor set by',
                                 'Kamoun et al.\ngenes: n =',
                                 length(unique(x$variable))),
                x = 'mean SHAP'))

# END --------

  mibc_splots <-
    mibc_splots[c("stats", "plots", "mean_plots")]

  insert_tail()
