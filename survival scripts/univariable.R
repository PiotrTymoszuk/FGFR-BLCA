# Univariable KM analysis if survival of tertile fo gene expression.
# Statistical significance is determined by Peto-Peto test corrected for
# multiple testing with the FDR method. The optimal cutoffs for stratification
# of the expression values are found with `surv_cutpoint`.
#
# Significant effects shared by at least two cohorts:
# unfavorable: TNFAIP6, GPC1, FIBP, FGF11, FGFBP1.

  insert_head()

# container --------

  surv_km <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals ---------

  insert_msg('Analysis globals')

  surv_km$variables <-
    surv_globals$variables[!stri_detect(surv_globals$variables, regex = 'sec$')]

  surv_km$data <- surv_globals$data

  surv_km$cutpoints <- surv_km$data %>%
    map(surv_cutpoint,
        time = 'os_days',
        event = 'death',
        variables = surv_km$variables,
        minprop = 0.1)

  surv_km$data <- surv_km$cutpoints %>%
    map(surv_categorize)

  ## setting factor levels according to the
  ## survminer's convention

  for(i in names(surv_km$data)) {

    surv_km$data[[i]][surv_km$variables] <-
      surv_km$data[[i]][surv_km$variables] %>%
      map_dfc(factor, c('high', 'low'))


  }

# Survfit objects ---------

  insert_msg('Survfit objects')

  for(i in names(surv_km$data)) {

    surv_km$survfit_obj[[i]] <- surv_km$variables %>%
      map(~paste('Surv(os_days, death) ~', .x)) %>%
      map(as.formula) %>%
      map(surv_fit,
          data = surv_km$data[[i]]) %>%
      set_names(surv_km$variables)

  }

# Numbers of patients and events in the strata --------

  insert_msg('Numbers of patients in the strata')

  for(i in names(surv_km$data)) {

    surv_km$n_numbers[[i]] <- surv_km$variables %>%
      set_names(surv_km$variables) %>%
      map(~split(surv_km$data[[i]], surv_km$data[[i]][[.x]])) %>%
      map(map,
          ~tibble(n_total = nrow(.x),
                  n_events = sum(.x$death))) %>%
      map(compress, names_to = 'strata') %>%
      compress(names_to = 'variable')

    ## ready-to-use legend labels for the color scale
    ## in the KM plots

    surv_km$axis_labs[[i]] <- surv_km$n_numbers[[i]] %>%
      mutate(axis_lab = paste(strata,
                              '\ntotal: n = ', n_total,
                              '\nevents: n = ', n_events)) %>%
      blast(variable) %>%
      map(~.x$axis_lab)

  }

# Median survival and tests --------

  insert_msg('Median survival and tests')

  for(i in names(surv_km$data)) {

    surv_km$stats[[i]] <- surv_km$survfit_obj[[i]] %>%
      map_dfr(surv_median) %>%
      mutate(variable = stri_split_fixed(strata,
                                         pattern = '=',
                                         simplify = TRUE)[, 1],
             strata = stri_split_fixed(strata,
                                       pattern = '=',
                                       simplify = TRUE)[, 2]) %>%
      as_tibble

    surv_km$test[[i]] <- surv_km$survfit_obj[[i]] %>%
      future_map_dfr(surv_pvalue,
                 method = 'S1',
                 .options = furrr_options(seed = TRUE)) %>%
      mutate(p_value = pval) %>%
      re_adjust(method = 'none') %>%
      mutate(raw_significance = significance) %>%
      re_adjust %>%
      mutate(plot_lab = paste('raw: ', raw_significance,
                              '\nFDR: ', significance)) %>%
      as_tibble

  }

# Significant effects --------

  insert_msg('Significant effects')

  ## single cohorts: significant upon FDR correction,
  ## common: significant effects shared by at least two cohorts

  surv_km$significant <- surv_km$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

  surv_km$common_significant <- surv_km$significant %>%
    safely(shared_features)(m = 2) %>%
    .$result %>%
    as.character

# Result table for the common significant effects -------

  insert_msg('Result table for the common significant effects')

  surv_km$result_tbl <-
    map2(surv_km$n_numbers,
         surv_km$stats,
         left_join,
         by = c('strata', 'variable')) %>%
    map(mutate,
        median = ifelse(is.na(median), '',
                        paste0(signif(median, 2), ' [95% CI: ',
                               signif(lower, 2), ' to ',
                               signif(upper, 2), ']')))

  surv_km$result_tbl <-
    map2(surv_km$result_tbl,
         surv_km$test,
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    arrange(variable, cohort, strata) %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    filter(variable %in% surv_km$common_significant)

  surv_km$result_tbl <- surv_km$result_tbl %>%
    select(variable, cohort, strata,
           n_total, n_events,
           median,
           raw_significance, significance) %>%
    set_names(c('Gene symbol', 'Cohort', 'Expression strata',
               'Patients, N', 'Events, N',
               'Median overall survival, days',
               'Raw significance', 'FDR significance'))

# Meta-data to be presented in the KM plots -------

  insert_msg('Meta-data for KM plots')

  ## testing results to be presented in the KM plots

  surv_km$plot_test <- surv_km$test %>%
    map(filter, variable %in% surv_km$common_significant) %>%
    map(mutate, variable = factor(variable, surv_km$common_significant)) %>%
    map(arrange, variable) %>%
    map(mutate, variable = as.character(variable))

  ## plot captions with the optimal cut points,
  ## and numbers of patients and events

  surv_km$plot_cap$cutpoints <- surv_km$cutpoints %>%
    map(~.x$cutpoint[surv_km$common_significant, 'cutpoint'])

  surv_km$plot_cap$n_numbers <- surv_km$data %>%
    map(~paste('total: n = ', nrow(.x), ', events: n = ', sum(.x$death)))

  surv_km$plot_cap <-
    map2(surv_km$plot_cap$cutpoints,
         surv_km$plot_cap$n_numbers,
         ~paste0('cutoff = ', signif(.x, 2), ', ', .y))

  ## labels for the color scale with numbers of patients
  ## and events in the strata

  surv_km$axis_labs <- surv_km$axis_labs %>%
    map(~.x[surv_km$common_significant])

# Kaplan-Meier plots for the common raw significant effects --------

  insert_msg('KM plots')

  for(i in names(surv_km$plot_test)) {

    surv_km$plots[[i]] <-
      list(fit = surv_km$survfit_obj[[i]][surv_km$plot_test[[i]]$variable],
           title = surv_km$plot_test[[i]]$variable %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           pval = surv_km$plot_test[[i]]$plot_lab,
           legend.title = surv_km$plot_test[[i]]$variable %>%
             html_italic) %>%
      pmap(ggsurvplot,
           pval.size = 2.5,
           pval.coord = c(0, 0.07),
           xlab = 'overall survival, days') %>%
      map(~.x$plot +
            globals$common_theme +
            theme(plot.title = element_markdown(),
                  legend.title = element_markdown()))

    ## setting the subtitles and color scales

    surv_km$plots[[i]] <-
      list(x = surv_km$plots[[i]],
           y = surv_km$plot_cap[[i]],
           z =  surv_km$axis_labs[[i]]) %>%
      pmap(function(x, y, z) x +
             labs(subtitle = y) +
             scale_color_manual(values = c('firebrick', 'steelblue'),
                                labels = z))

  }

# END ------

  surv_km <-
    surv_km[c("cutpoints",
              "n_numbers", "axis_labs",
              "stats", "test", "result_tbl",
              "significant", "common_significant",
              "plots")]

  plan('sequential')

  insert_tail()
