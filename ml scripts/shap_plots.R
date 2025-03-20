# Plots of SHAP variable importance

  insert_head()

# container ------

  ml_splots <- list()

# SHAP and SHAPviz objects -------

  insert_msg('SHAP and SHAPviz objects')

  ## class names and algorithm labels

  ml_splots$class_levels <- ml_globals$class_levels

  ## the SHAP/SHAPviz objects

  ml_splots$shap_obj <- ml_shap$shap_obj[c("elnet", "ranger")]

  ml_splots$shapviz <- ml_splots$shap_obj %>%
    map(shapviz) %>%
    map(set_names, ml_splots$class_levels)

# Mean absolute SHAP values --------

  insert_msg('Mean absolute SHAP values')

  ml_splots$stats <- ml_splots$shapviz %>%
    map(get_shap_values) %>%
    map(map, abs) %>%
    map(map, colMeans) %>%
    map(map,
        compress,
        names_to = 'variable',
        values_to = 'importance') %>%
    map(compress, names_to = 'class')

# Beeswarms of SHAPs and bar plots of mean absolute SHAPS ------

  insert_msg('Plots')

  ## version for the report

  ml_splots$plots <-
    list(shapviz_obj = ml_splots$shapviz,
         title_suffix = globals$algo_labs[names(ml_splots$shapviz)]) %>%
    pmap(my_shap_bee,
         panels = TRUE,
         add_legend = FALSE,
         top_variables = 15,
         point_size = 1,
         cust_theme = globals$common_theme)

  ## version for the paper with larger fonts

  ml_splots$paper_plots <-
    list(shapviz_obj = ml_splots$shapviz) %>%
    pmap(my_shap_bee,
         panels = TRUE,
         add_legend = FALSE,
         top_variables = 10,
         point_size = 2,
         n_breaks = 3,
         cust_theme = globals$figure_theme +
           theme(plot.title = element_text(size = 12, face = 'bold'),
                 axis.text.x = element_text(size = 9),
                 axis.text.y = element_text(size = 11, face = 'italic')))

# END --------

  ml_splots <- ml_splots[c("stats", "plots", "paper_plots")]

  insert_tail()
