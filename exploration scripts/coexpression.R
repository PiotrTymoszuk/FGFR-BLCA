# Correlation of log2-transformed mRNA levels of FGF and FGFR genes
# in the cancer tissue and DepMap cell lines.
#
# Relevant effects are defined as Spearman's rho >= 0.3
# and pFDR < 0.05 in permutation test.

  insert_head()

# container -------

  expl_corr <- list()

# Parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data -------

  insert_msg('Analysis data')

  ## analysis variables

  expl_corr$variables <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    reduce(union)

  ## analysis data, appending with the DepMap cell lines

  expl_corr$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    c(list(depmap = depmap$expression)) %>%
    map(select,
        any_of(expl_corr$variables)) %>%
    map(~filter(.x, complete.cases(.x)))

  expl_corr$variables <- expl_corr$data %>%
    map(names) %>%
    reduce(intersect)

  expl_corr$data <- expl_corr$data %>%
    map(~.x[expl_corr$variables])

# Spearman's correlation -------

  insert_msg("Spearman's correlation")

  expl_corr$test <- expl_corr$data %>%
    future_map(f_cor_test,
               type = 'permutation',
               method = 'spearman',
               n_iter = 1000,
               as_data_frame = TRUE,
               .options = furrr_options(seed = TRUE))

  ## formatting the results, excluding self- and repeated correlations

  expl_corr$test <- expl_corr$test %>%
    map(mutate,
        corr_id = map2(variable1, variable2, c) %>%
          map(sort) %>%
          map_chr(paste, collapse = '_')) %>%
    map(filter,
        variable1 != variable2,
        !duplicated(corr_id)) %>%
    map(re_adjust) %>%
    map(p_formatter, text = TRUE) %>%
    map(mutate,
        variable1 = factor(variable1, expl_corr$variables),
        variable2 = factor(variable2, expl_corr$variables),
        correlation = ifelse(p_adjusted >= 0.05, 'ns',
                             ifelse(rho >= 0.3, 'positive',
                                    ifelse(rho <= -0.3, 'negative', 'ns'))),
        correlation = factor(correlation, c('positive', 'negative', 'ns'))) %>%
    map(as_tibble)

# Significant effects -------

  insert_msg('significant effects')

  expl_corr$significant <- expl_corr$test %>%
    map(filter, correlation %in% c('positive', 'negative')) %>%
    map(blast, correlation) %>%
    transpose %>%
    map(map, ~ .x$corr_id)

  ## effects shared by all data sets

  expl_corr$common_significant <- expl_corr$significant %>%
    map(reduce, intersect)

  ## plot metadata

  expl_corr$plot_data <- expl_corr$test %>%
    map(filter,
        corr_id %in% reduce(expl_corr$common_significant, union)) %>%
    map(mutate,
        plot_cap = paste0('n = ', n,
                          ', \u03C1 = ', signif(rho, 2),
                          ', ', significance))

  expl_corr$common_limits <-
    list(expl_corr$plot_data[[1]]$variable1,
         expl_corr$plot_data[[2]]$variable2) %>%
    map(droplevels) %>%
    map(as.character)

# Bubble plots for the shared significant effects -------

  insert_msg('Bubble plots with the shared significant effects')

  expl_corr$bubble_plots <-
    list(x = expl_corr$plot_data,
         y = globals$cohort_labs[names(expl_corr$test)],
         z = paste("Spearman's correlation, n =",
                   map_dbl(expl_corr$data, nrow))) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = variable1,
                      y = variable2,
                      fill = factor(sign(rho)),
                      size = abs(rho))) +
           geom_point(shape = 21) +
           geom_text(aes(label = signif(rho, 2),
                         color = factor(sign(rho))),
                     size = 2.5,
                     vjust = -1.4,
                     hjust = 0.5,
                     show.legend = FALSE) +
           scale_fill_manual(values = c('-1' = 'steelblue',
                                        '0' = 'gray60',
                                        '1' = 'firebrick'),
                             labels = c('-1' = 'negative',
                                        '0' = 'ns',
                                        '1' = 'positive'),
                             name = 'Correlation') +
           scale_color_manual(values = c('-1' = 'steelblue',
                                         '0' = 'gray60',
                                         '1' = 'firebrick'),
                              labels = c('-1' = 'negative',
                                         '0' = 'ns',
                                         '1' = 'positive'),
                              name = 'Correlation') +
           scale_size_area(max_size = 5,
                           limits = c(0, 1),
                           name = 'abs(\u03C1)') +
           globals$common_theme +
           theme(axis.title = element_blank(),
                 axis.text = element_text(face = 'italic')) +
           labs(title = y,
                subtitle = z))

# Scatter plots for the shared significant correlations -------

  insert_msg('Scatter plots for the shared significant correlations')

  for(i in names(expl_corr$plot_data)) {

    expl_corr$scatter_plots[[i]] <-
      list(variables = map2(expl_corr$plot_data[[i]]$variable1,
                            expl_corr$plot_data[[i]]$variable2,
                            c),
           plot_title = map2(expl_corr$plot_data[[i]]$variable1,
                             expl_corr$plot_data[[i]]$variable2,
                             ~paste(html_italic(.x), html_italic(.x),
                                    sep = ' and ')) %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           plot_subtitle = expl_corr$plot_data[[i]]$plot_cap,
           x_lab = expl_corr$plot_data[[i]]$variable1 %>%
             html_italic %>%
             paste0(', log<sub>2</sub> expression'),
           y_lab = expl_corr$plot_data[[i]]$variable2 %>%
             html_italic %>%
             paste0(', log<sub>2</sub> expression')) %>%
      pmap(plot_correlation,
           data = expl_corr$data[[i]],
           type = 'correlation',
           point_hjitter = 0.01,
           point_wjitter = 0.01,
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   axis.title.x = element_markdown(),
                   axis.title.y = element_markdown(),
                   plot.tag = element_blank()),
           show_trend = FALSE) %>%
      map(~.x +
            geom_smooth(method = 'loess',
                        se = FALSE,
                        span = 0.5)) %>%
      set_names(expl_corr$plot_data[[i]]$corr_id)

  }

# Result table -------

  insert_msg('Result table')

  ## for the shared significant effects

  expl_corr$result_tbl <- expl_corr$test %>%
    map(filter, corr_id %in% reduce(expl_corr$common_significant, union)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           rho = signif(rho, 2)) %>%
    select(cohort,
           variable1,
           variable2,
           rho,
           significance) %>%
    set_names(c('Cohort',
                'Gene 1, symbol',
                'Gene 2, symbol',
                'Correlation coefficient, rho',
                'Significance'))

# END ----

  expl_corr$plot_data <- NULL
  expl_corr$common_limits <- NULL
  expl_corr$data <- NULL

  expl_corr <- compact(expl_corr)

  plan('sequential')

  insert_tail()
