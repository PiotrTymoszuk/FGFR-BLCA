# Importance of the explanatory variables in the Elastic Net and Random forest
# models of consensus classes.
# The models use solely FGF/FGFR/FGFBP expression as explanatory factors.
#
# The importance measures are simply abs(beta) of the model coefficients
# for the Elastic Net model.
# For the Random Forest model, we're calculating the perumtation importance.

  insert_head()

# container ------

  ml_imp <- list()

# Ranger model with permutation importance -------

  insert_msg('Random forest model with permutation importance')

  ## model

  ml_imp$ranger_model <- ml_globals$train_data %>%
    select(-consensusClass) %>%
    as.matrix %>%
    call2(.fn = 'ranger',
          x = .,
          y = ml_globals$train_data$consensusClass,
          !!!ml_train$tune_obj$ranger$best_tune,
          probability = TRUE,
          importance = 'permutation') %>%
    eval

  ## permutation importance

  ml_imp$ranger_importance <- ml_imp$ranger_model %>%
    importance %>%
    compress(names_to = 'variable',
             values_to = 'importance')

# Random forest importance in a bar plot --------

  insert_msg('Permutation importance plot, Random Forest')

  ml_imp$ranger_plot <- ml_imp$ranger_importance %>%
    slice_max(importance, n = 15) %>%
    ggplot(aes(x = importance,
               y = reorder(variable, importance))) +
    geom_bar(stat = 'identity',
             color = 'white',
             fill = 'aquamarine4') +
    globals$common_theme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = 'italic')) +
    labs(title = 'Random Forest',
         x = 'permutation importance, \u0394 error')

# Elastic Net model coefficients --------

  insert_msg('Elastic Net model coefficients')

  ml_imp$coefs <- ml_train$tune_obj$elnet$model %>%
    coef %>%
    map(as.matrix) %>%
    map(as.data.frame) %>%
    map(set_names, 'beta') %>%
    map(rownames_to_column, 'variable') %>%
    map(filter,
        !variable %in% c('', '(Intercept)')) %>%
    map(mutate, abs_beta = abs(beta)) %>%
    map(as_tibble)

# Bar plots of the Elastic Net beta coefficients -------

  insert_msg('Bar plots of the Elastic Net beta coefficients')

  ## the top mos important factors are presented

  ml_imp$coef_plots <- ml_imp$coefs %>%
    map(filter, variable != 0) %>%
    map(slice_max, abs_beta, n = 15) %>%
    list(x = .,
         y = paste(names(ml_imp$coefs),
                   ', Elastic Net')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = beta,
                      y = reorder(variable, beta),
                      fill = factor(sign(beta)))) +
           geom_bar(stat = 'identity',
                    color = 'white') +
           geom_vline(xintercept = 0,
                      linetype = 'solid',
                      linewidth = 0.25) +
           scale_fill_manual(values = c('-1' = 'steelblue',
                                        '1' = 'firebrick'),
                             labels = c('-1' = 'negative',
                                        '1' = 'positive'),
                             name = 'association\nsign') +
           globals$common_theme +
           theme(axis.text.y = element_text(face = 'italic'),
                 axis.title.y = element_blank()) +
           labs(title = y,
                x = expression('importance, ' * beta[ElasticNet])))

# END ------

  ml_imp$ranger_model <- NULL

  ml_imp <- compact(ml_imp)

  insert_tail()
