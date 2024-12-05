# Variable importance for the Elastic Net Cox model of overall survival

  insert_head()

# container -------

  surv_imp <- list()

# model coefficients: they translate to log HR --------

  insert_msg('Model coefficients')

  surv_imp$coefs <- surv_train$tune_obj$model %>%
    coef %>%
    as.matrix %>%
    as.data.frame %>%
    set_names('beta') %>%
    rownames_to_column('variable') %>%
    as_tibble

  surv_imp$coefs <- surv_imp$coefs %>%
    filter(beta != 0) %>%
    mutate(hr = exp(beta))

# Bar plots of the HR values ---------

  insert_msg('Bar plots of the HR values')

  surv_imp$coef_plot <- surv_imp$coefs %>%
    ggplot(aes(x = beta,
               y = reorder(html_italic(variable), beta),
               fill = factor(sign(beta)),
               color = factor(sign(beta)))) +
    geom_bar(stat = 'identity',
             color = 'white') +
    geom_vline(xintercept = 0,
               linetype = 'solid',
               linewidth = 0.25) +
    geom_text(aes(label = signif(hr, 3),
                  x = 1.05 * beta,
                  hjust = ifelse(beta < 0, 1, 0)),
              size = 2.5,
              show.legend = FALSE) +
    scale_y_discrete(labels = function(x) stri_replace(x,
                                                       fixed = '_sec',
                                                       replacement = html_sup('2'))) +
    scale_fill_manual(values = c('-1' = 'steelblue',
                                 '1' = 'firebrick'),
                      labels = c('-1' = 'favorable',
                                 '1' = 'unfavorable'),
                      name = 'risk\nassociation') +
    scale_color_manual(values = c('-1' = 'steelblue',
                                 '1' = 'firebrick'),
                      labels = c('-1' = 'favorable',
                                 '1' = 'unfavorable'),
                      name = 'risk\nassociation') +
    globals$common_theme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_markdown()) +
    labs(title = 'Overall survival, variable importance',
         x = 'log HR, per SD')

# END -------

  insert_tail()
