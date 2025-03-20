# Distribution tests for predictions of drug resistance in form of log IC50
# (GDSC1, GDSC2, CTRP2) or area under the dose - response curve (AUC).
# Normality is tested by Shapiro-Wilk test.
#
# Inter-correlations between the drug predictions are investigated by principal
# component analysis.
# Correlations between predictions of resistance to anti-FGFR drugs are assessed
# by Spearman's permutation test.
# We're also extracting non-zero coefficients beta from the regularized linear
# models used for predictions of resistance to the anti-FGFR drugs.
# The top predictors with the largest abs(bet) are visualized.

  insert_head()

# container -------

  expl_drugs <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals -------

  insert_msg('Analysis globals')

  # variables, lexicons, and resistance metrics

  expl_drugs$lexicon <- drugs$lexicon

  expl_drugs$variables <- expl_drugs$lexicon$variable

  expl_drugs$data <- drugs$predictions %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(expl_drugs$variables))

  ## drugs that target FGFR

  expl_drugs$fgfr_drugs <- drugs$fgfr_drugs

# normality testing ------

  insert_msg('Normality testing')

  ## for multiple drug predictions in the BCAN cohort,
  ## there's clearly no normality. Those include Mitomycin-C,
  ## Gemcitabine, Bleomycin, and Metothrexate

  expl_drugs$norm_test <- expl_drugs$data %>%
    map(~f_shapiro_test(.x[expl_drugs$variables],
                        as_data_frame = TRUE)) %>%
    map(mutate,
        violated = w < 0.9) %>%
    map(arrange, p_value) %>%
    map(as_tibble)

# PCA of the drug predictions ---------

  insert_msg('PCA for the drug predictions')

  ## 3 - 4 PC explain most of the variance
  ## setting proper variable names

  expl_drugs$pca_obj <- expl_drugs$data %>%
    map(zScores) %>%
    map(reduce_data,
        kdim = 6,
        red_fun = 'pca')

  for(i in names(expl_drugs$pca_obj)) {

    expl_drugs$pca_obj[[i]]$loadings <-
      expl_drugs$pca_obj[[i]]$loadings %>%
      mutate(variable = expl_drugs$variables)

  }

  ## scree plots

  expl_drugs$scree_plots <- expl_drugs$pca_obj %>%
    map(plot, type = 'scree')

  expl_drugs$scree_plots <-
    list(x = expl_drugs$scree_plots,
         y = globals$cohort_labs[names(expl_drugs$scree_plots)]) %>%
    pmap(function(x, y) x  +
           globals$common_theme +
           labs(title = y)) %>%
    map(tag2subtitle)

  ## loadings and top eigenvector variables

  expl_drugs$loadings <- expl_drugs$pca_obj %>%
    map(~.x$loadings)

  for(i in names(expl_drugs$loadings)) {

    expl_drugs$top_loadings[[i]] <- paste0('comp_', 1:6) %>%
      set_names(paste0('comp_', 1:6)) %>%
      map(~slice_max(expl_drugs$loadings[[i]],
                     .data[[.x]],
                     n = 10)) %>%
      map(~.x$variable)

  }

  ## plots of loadings, first two dimensions:
  ## point color codes for the experiment

  expl_drugs$loading_plots <-
    map2(expl_drugs$loadings,
         expl_drugs$top_loadings %>%
           map(~.x[1:2]) %>%
           map(reduce, union),
         ~filter(.x, variable %in% .y)) %>%
    map(mutate,
        variable_label = exchange(variable,
                                  expl_drugs$lexicon),
        experiment = exchange(variable,
                              expl_drugs$lexicon,
                              value = 'experiment'))

  expl_drugs$loading_plots <-
    list(x = expl_drugs$loading_plots,
         y = globals$cohort_labs[names(expl_drugs$loading_plots)]) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = comp_1,
                      y = comp_2,
                      fill = experiment)) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = 0,
                      linetype = 'dashed') +
           geom_segment(aes(x = 0,
                            y = 0,
                            xend = comp_1,
                            yend = comp_2),
                        color = 'gray60') +
           geom_point(size = 2,
                      shape = 21) +
           geom_text_repel(aes(label = variable_label,
                               color = experiment),
                           size = 2.5,
                           show.legend = FALSE) +
           scale_fill_manual(values = globals$screen_colors,
                             name = 'screening\nexperiment') +
           scale_color_manual(values = globals$screen_colors,
                              name = 'screening\nexperiment') +
           globals$common_theme +
           labs(title = y,
                x = 'PC1',
                y = 'PC2'))

# Correlations of drugs that target FGFR --------

  insert_msg('Correlations of drugs that target FGFR')

  ## Spearman's permutation tests;
  ## FDR correction and identification of significant correlations
  ## with moderate-to-large effect size

  expl_drugs$fgfr_cor_test <- expl_drugs$data %>%
    map(~.x[expl_drugs$fgfr_drugs]) %>%
    future_map(f_cor_test,
               method = 'spearman',
               type = 'permutation',
               n_iter = 1000,
               as_data_frame = TRUE,
               .options = furrr_options(seed = TRUE)) %>%
    map(as_tibble)

  expl_drugs$fgfr_cor_test <- expl_drugs$fgfr_cor_test %>%
    map(re_adjust) %>%
    map(mutate,
        variable_label1 = exchange(variable1, expl_drugs$lexicon),
        variable_label2 = exchange(variable2, expl_drugs$lexicon),
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(rho >= 0.3, 'positive',
                                   ifelse(rho <= -0.3, 'negative', 'ns'))),
        regulation = factor(regulation, c('positive', 'negative', 'ns')),
        variable1 = factor(variable1, expl_drugs$fgfr_drugs),
        variable2 = factor(variable2, expl_drugs$fgfr_drugs))

  ## bubble plots

  expl_drugs$fgfr_bubble_plots <-
    list(x = expl_drugs$fgfr_cor_test %>%
           map(filter, variable1 != variable2),
         y = globals$cohort_labs[names(expl_drugs$fgfr_cor_test)]) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = variable1,
                  y = variable2,
                  color = regulation,
                  size = abs(rho))) +
           geom_point(shape = 16) +
           geom_text(aes(label = ifelse(regulation %in% c('positive', 'negative'),
                                        signif(rho, 2), NA)),
                     size = 2.5,
                     hjust = 0.5,
                     vjust = -1.1,
                     show.legend = FALSE) +
           scale_size_area(max_size = 5,
                           limits = c(0, 1)) +
           scale_color_manual(values = c(positive = 'firebrick',
                                         negative = 'steelblue',
                                         ns = 'gray70'),
                              name = 'correlation',
                              drop = FALSE) +
           scale_x_discrete(labels = function(x) exchange(x,
                                                          expl_drugs$lexicon,
                                                          value = 'hm_label')) +
           scale_y_discrete(labels = function(x) exchange(x,
                                                          expl_drugs$lexicon,
                                                          value = 'hm_label')) +
           guides(x = guide_axis(angle = 90)) +
           globals$common_theme +
           theme(axis.text.x = element_markdown(),
                 axis.text.y = element_markdown(),
                 axis.title = element_blank())  +
           labs(title = y,
                subtitle = paste('n =', x$n[1])))

# Key genes for predictions of the FGFR drug resistance ------

  insert_msg('Key genes for predictions of resistance of FGFR inhibitors')

  ## non-null coefficients of linear models for prediction
  ## of resistance towards FGFR-targeting drugs

  expl_drugs$coefs <- list(gdsc1 = gdsct$train_objects$gdsc1,
                           gdsc2 = gdsct$train_objects$gdsc2,
                           ctrp2 = ctrp2$train_object,
                           prism = prism$train_object) %>%
    map(~.x$models) %>%
    map(~.x[names(.x) %in% expl_drugs$fgfr_drugs]) %>%
    map(map, ~coef(.x, .x$lambda.min))

  expl_drugs$coefs <- expl_drugs$coefs %>%
    reduce(c) %>%
    map(as.matrix) %>%
    map(as.data.frame) %>%
    map(set_names, 'beta') %>%
    map(rownames_to_column, 'variable') %>%
    map(filter,
        beta != 0,
        variable != '(Intercept)') %>%
    map(function(x) if(nrow(x) == 0) NULL else x) %>%
    compact %>%
    map(arrange, -beta) %>%
    map(as_tibble)

  expl_drugs$coefs <- expl_drugs$coefs %>%
    map(mutate,
        regulation = ifelse(beta > 0, 'positive', 'negative'),
        regulation = factor(regulation, c('positive', 'negative')))

  ## the key genes shared by at least two cohorts

  expl_drugs$coef_counts <- expl_drugs$coefs %>%
    map(blast, regulation) %>%
    map(map, ~.x$variable) %>%
    transpose %>%
    map(count_features)

  expl_drugs$common_coefs <- expl_drugs$coefs %>%
    map(blast, regulation) %>%
    map(map, ~.x$variable) %>%
    transpose %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Plots of betas for the key genes --------

  insert_msg('Plots of betas for the key genes')

  expl_drugs$top_coef_plots <-
    list(x = expl_drugs$coefs %>%
           map(slice_max, abs(beta), n = 20),
         y = exchange(names(expl_drugs$coefs),
                      expl_drugs$lexicon,
                      value = 'plot_title')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = beta,
                      y = reorder(variable, beta),
                      fill = regulation)) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_bar(stat = 'identity',
                    color = 'white') +
           scale_fill_manual(values = c(positive = 'firebrick',
                                        negative = 'steelblue'),
                             name = 'association\nwith resistance') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_text(face = 'italic')) +
           labs(title = y,
                subtitle = 'key predictor genes, top 20 abs(\u03B2)',
                x = 'model coefficient estimate, (\u03B2)'))

# Word cloud plots for the predictor genes shared by at least two models --------

  insert_msg('Word clouds for the top predictor genes')

  expl_drugs$common_coef_wordcloud <- expl_drugs$coef_counts %>%
    map(filter, n >= 2) %>%
    compress(names_to = 'regulation') %>%
    mutate(regulation = ifelse(regulation == 'negative',
                               -1, 1),
           element = as.character(element)) %>%
    arrange(n) %>%
    plot_wordcloud(label_variable = 'element',
                   size_variable = 'n',
                   color_variable = 'regulation',
                   size_range = c(1.5, 4),
                   size_lab = '# models',
                   color_lab = 'association\nresistance',
                   shape = 'square',
                   fraction_rotated = 0.3,
                   fontface = 'italic',
                   show.legend = TRUE) +
    theme(legend.text = globals$common_text,
          legend.title = globals$common_text,
          plot.title = element_text(size = 8,
                                    face = 'bold')) +
    labs(title = 'Common predictors of anti-FGFR drug resistance')

# END ------

  expl_drugs$data <- NULL

  expl_drugs <- compact(expl_drugs)

  insert_tail()
