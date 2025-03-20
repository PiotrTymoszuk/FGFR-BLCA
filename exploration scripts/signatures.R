# Normality of ssGSEA scores of FGFR-related gene signatures in the TCGA,
# IMvigor, and BCAN cohorts. Investigated by Shapiro-Wilk tests

  insert_head()

# container -----

  expl_sig <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  # variables, lexicons, and resistance metrics

  expl_sig$lexicon <- sig$lexicon

  expl_sig$variables <- expl_sig$lexicon$variable

  expl_sig$data <- sig$scores %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(expl_sig$variables))

# normality testing ------

  insert_msg('Normality testing')

  ## the normality assumption is more or less met (as expected!)

  expl_sig$norm_test <- expl_sig$data %>%
    map(~f_shapiro_test(.x[expl_sig$variables],
                        as_data_frame = TRUE)) %>%
    map(mutate,
        violated = w < 0.9) %>%
    map(arrange, p_value) %>%
    map(as_tibble)

# PCA of the signature scores ---------

  insert_msg('PCA for the signature scores')

  ## 2 - 3 PC explain most of the variance
  ## setting proper variable names

  expl_sig$pca_obj <- expl_sig$data %>%
    map(zScores) %>%
    map(reduce_data,
        kdim = 6,
        red_fun = 'pca')

  for(i in names(expl_sig$pca_obj)) {

    expl_sig$pca_obj[[i]]$loadings <-
      expl_sig$pca_obj[[i]]$loadings %>%
      mutate(variable = expl_sig$variables)

  }

  ## scree plots

  expl_sig$scree_plots <- expl_sig$pca_obj %>%
    map(plot, type = 'scree')

  expl_sig$scree_plots <-
    list(x = expl_sig$scree_plots,
         y = globals$cohort_labs[names(expl_sig$scree_plots)]) %>%
    pmap(function(x, y) x  +
           globals$common_theme +
           labs(title = y)) %>%
    map(tag2subtitle)

  ## loadings and top eigenvector variables

  expl_sig$loadings <- expl_sig$pca_obj %>%
    map(~.x$loadings)

  for(i in names(expl_sig$loadings)) {

    expl_sig$top_loadings[[i]] <- paste0('comp_', 1:6) %>%
      set_names(paste0('comp_', 1:6)) %>%
      map(~slice_max(expl_sig$loadings[[i]],
                     .data[[.x]],
                     n = 10)) %>%
      map(~.x$variable)

  }

  ## plots of the top loadings, first two dimensions:
  ## point color codes for the experiment

  expl_sig$top_loading_plots <-
    map2(expl_sig$loadings,
         expl_sig$top_loadings %>%
           map(~.x[1:2]) %>%
           map(reduce, union),
         ~filter(.x, variable %in% .y)) %>%
    map(mutate,
        variable_label = exchange(variable,
                                  expl_sig$lexicon))

  expl_sig$top_loading_plots <-
    list(x = expl_sig$loading_plots,
         y = globals$cohort_labs[names(expl_sig$loading_plots)]) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = comp_1,
                      y = comp_2)) +
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
                      shape = 21,
                      fill = 'steelblue') +
           geom_text_repel(aes(label = variable_label),
                           size = 2.5,
                           show.legend = FALSE) +
           globals$common_theme +
           labs(title = y,
                x = 'PC1',
                y = 'PC2'))

# END -----

  expl_sig$lexicon <- NULL
  expl_sig$data <- NULL

  expl_sig <- compact(expl_sig)

  insert_tail()
