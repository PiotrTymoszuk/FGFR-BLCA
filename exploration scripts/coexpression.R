# Correlation of log2-transformed mRNA levels of FGF and FGFR genes
# in the cancer tissue.
#
# Relevant effects are defined as r >= 0.3 and pFDR < 0.05.

  insert_head()

# container -------

  expl_corr <- list()

# Parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data -------

  insert_msg('Analysis data')

  ## analysis variables

  expl_corr$variables <- globals[c("receptors", "ligands")]

  ## analysis data

  expl_corr$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(select,
        sample_id, any_of(reduce(expl_corr$variables, union))) %>%
    map(~filter(.x, complete.cases(.x)))

  expl_corr$variables <- expl_corr$variables %>%
    map(~.x[.x %in% names(expl_corr$data[[1]])])

  expl_corr$variable_pairs <- expl_corr$variables %>%
    reduce(union) %>%
    combn(m = 2, simplify = FALSE)

# Pearson's correlation -------

  insert_msg("Pearson's correlation")

  for(i in names(expl_corr$data)) {

    expl_corr$test[[i]] <- expl_corr$variable_pairs %>%
      future_map_dfr(~correlate_variables(expl_corr$data[[i]],
                                          variables = .x,
                                          what = 'correlation',
                                          type = 'pearson',
                                          ci = TRUE),
                     .options = furrr_options(seed = TRUE)) %>%
      re_adjust(method = 'BH') %>%
      mutate(corr_id = paste(variable1, variable2, sep = '_'),
             variable1 = factor(variable1, reduce(expl_corr$variables, c)),
             variable2 = factor(variable2, reduce(expl_corr$variables, c)))

  }

# Significant effects -------

  insert_msg('significant effects')

  expl_corr$significant <- expl_corr$test %>%
    map(filter,
        p_adjusted < 0.05,
        abs(estimate) >= 0.3) %>%
    map(~.x$corr_id)

  ## effects shared by both cohorts

  expl_corr$common_significant <- expl_corr$significant %>%
    reduce(intersect)

# Bubble plots for the shared significant effects -------

  insert_msg('Bubble plot with the shared significant effects')

  ## plot metadata

  expl_corr$plot_data <- expl_corr$test %>%
    map(filter, corr_id %in% expl_corr$common_significant)

  expl_corr$common_limits <-
    union(expl_corr$plot_data[[1]]$variable1,
          expl_corr$plot_data[[2]]$variable2) %>%
    factor(reduce(expl_corr$variables, c)) %>%
    droplevels()

  ## bubble plots

  expl_corr$bubble_plots <-
    list(x = expl_corr$plot_data,
         y = globals$cohort_labs[names(expl_corr$test)],
         z = paste("Pearson's correlation, n =",
                   map_dbl(expl_corr$data, nrow))) %>%
    pmap(function(x, y, z) x %>%
           ggplot(aes(x = variable1,
                      y = variable2,
                      fill = factor(sign(estimate)),
                      size = abs(estimate))) +
           geom_point(shape = 21) +
           geom_text(aes(label = signif(estimate, 2),
                         color = factor(sign(estimate))),
                     size = 2.5,
                     vjust = -1.2,
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
                           name = "r") +
           scale_x_discrete(limits = expl_corr$common_limits) +
           scale_y_discrete(limits = expl_corr$common_limits) +
           globals$common_theme +
           theme(axis.title = element_blank(),
                 axis.text = element_text(face = 'italic')) +
           labs(title = y,
                subtitle = z))

# Result table -------

  insert_msg('Result table')

  ## for the shared significant effects

  expl_corr$result_tbl <- expl_corr$test %>%
    map(filter, corr_id %in% expl_corr$common_significant) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           est_lab = paste0(signif(estimate, 2),
                            ' [95% CI: ', signif(lower_ci, 2),
                            ' to ', signif(upper_ci, 2), ']')) %>%
    select(cohort,
           variable1,
           variable2,
           est_lab,
           significance) %>%
    set_names(c('Cohort',
                'Gene 1, symbol',
                'Gene 2, symbol',
                'Correlation coefficient with 95% confidence interval',
                'Significcance'))

# END ----

  expl_corr$plot_data <- NULL
  expl_corr$common_limits <- NULL
  expl_corr$data <- NULL

  expl_corr <- compact(expl_corr)

  plan('sequential')

  insert_tail()
