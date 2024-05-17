# Predictions of the genetic clusters defined by latent class modeling in the
# GENIE cohort

  insert_head()

# container ------

  lca_pred <- list()

# analysis data and variables -------

  insert_msg('Analysis data and variables')

  ## data used for modeling

  lca_pred$factor_data <- lca_globals$factor_data %>%
    map(column_to_rownames, 'sample_id')

  lca_pred$modeling_variables <- lca_globals$modeling_variables

  ## the whole data sets of genetic features

  lca_pred$data <- lca_globals$data %>%
    map(column_to_rownames, 'sample_id')

  lca_pred$variables <- lca_globals$variables

# The model of interest and predictions -------

  insert_msg('6-class model and predictions')

  lca_pred$model <- lca_dev$tune_object$models$class_6

  lca_pred$model$class_names <-
    c('mutRB1',
      'hyperMut',
      'mutFGFR3',
      'oligoMut',
      'ampMDM2',
      'del9p21')

  ## predictions

  lca_pred$posteriors$genie <- lca_pred$model %>%
    predict

  lca_pred$posteriors[c('msk', 'tcga')] <-
    lca_pred$factor_data[c("msk", "tcga")] %>%
    map(predict,
        object = lca_pred$model)

# Diagnostic plots: class-assignment posterior p -------

  insert_msg('Plots of class assignment probabilities')

  lca_pred$posterior_plots <- lca_pred$posteriors %>%
    map(plot,
        type = 'class_posterior',
        cust_theme = globals$common_theme,
        flip = TRUE,
        point_size = 1,
        point_alpha = 0.5)

  lca_pred$posterior_plots <-
    list(x = lca_pred$posterior_plots,
         y = globals$cohort_labs[names(lca_pred$posterior_plots)]) %>%
    pmap(function(x, y) x  +
           scale_fill_manual(values = globals$genet_colors) +
           labs(title = paste('Class assignment probailities,', y),
                x = 'Cancer sample') +
           theme(panel.grid.major.x = element_blank()))

# Clustering object for all genetic features of interest ------

  insert_msg('Clustering object')

  lca_pred$genet_objects <-
    map2(lca_pred$data,
         lca_pred$posteriors,
         ~list(data = quo(!!.x),
               dist_mtx = calculate_dist(.x, 'smc'),
               dist_method = 'smc',
               clust_obj = NULL,
               clust_fun = 'prediction',
               clust_assignment = set_names(.y[c('observation', 'class')],
                                            c('observation', 'clust_id')),
               dots = list())) %>%
    map(clust_analysis)

# Ready-to-use assignment data frames -------

  insert_msg('Ready-to-use assignment tables')

  lca_pred$assignment <- lca_pred$genet_objects %>%
    map(extract, 'assignment') %>%
    map(transmute,
        sample_id = observation,
        clust_id = clust_id)

# END -----

  #lca_pred$data <- NULL
  lca_pred$factor_data <- NULL

  lca_pred <- compact(lca_pred)

  insert_tail()
