# Plots of SHAP variable importance

  insert_head()

# container ------

  ml_splots <- list()

# SHAP and SHAPviz objects -------

  insert_msg('SHAP and SHAPviz objects')

  ## class names and algorithm labels

  ml_splots$class_levels <- ml_globals$class_levels

  ml_splots$algo_labs <-
    c('elnet' = 'Elastic Net',
      'ranger' = 'Random Forest')

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

  ml_splots$plots <-
    list(shapviz_obj = ml_splots$shapviz,
         title_suffix = ml_splots$algo_labs) %>%
    pmap(my_shap_bee,
         panels = TRUE,
         add_legend = FALSE,
         top_variables = 15,
         point_size = 1)

# END --------

  ml_splots <- ml_splots[c("stats", "plots")]

  insert_tail()
