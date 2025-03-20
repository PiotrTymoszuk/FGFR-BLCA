# Importance of the explanatory variables in the Elastic Net and Random forest
# models of consensus classes.
# The models use the full set of gene proposed by Kamoun et al.
# The importance measures are simply abs(beta) of the model coefficients
# for the Elastic Net model.
# For the Random Forest model, we're calculating the permutation importance.

  insert_head()

# container ------

  mibc_imp <- list()

# Ranger model with permutation importance -------

  insert_msg('Random forest model with permutation importance')

  ## model

  mibc_imp$ranger_model <- mibc_globals$train_data %>%
    select(-consensusClass) %>%
    as.matrix %>%
    call2(.fn = 'ranger',
          x = .,
          y = mibc_globals$train_data$consensusClass,
          !!!mibc_train$tune_obj$ranger$best_tune,
          probability = TRUE,
          importance = 'permutation') %>%
    eval

  ## permutation importance

  mibc_imp$ranger_importance <- mibc_imp$ranger_model %>%
    importance %>%
    compress(names_to = 'variable',
             values_to = 'importance')

# Random forest importance in a bar plot --------

  insert_msg('Permutation importance plot, Random Forest')

  mibc_imp$ranger_plot <- mibc_imp$ranger_importance %>%
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

  mibc_imp$coefs <- mibc_train$tune_obj$elnet$model %>%
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

  mibc_imp$coef_plots <- mibc_imp$coefs %>%
    map(filter, variable != 0) %>%
    map(slice_max, abs_beta, n = 15) %>%
    list(x = .,
         y = paste(names(mibc_imp$coefs),
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

  mibc_imp$ranger_model <- NULL

  mibc_imp <- compact(mibc_imp)

  insert_tail()
