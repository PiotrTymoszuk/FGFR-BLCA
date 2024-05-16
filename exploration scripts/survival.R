# Comparison of overall survival between the cohorts.

  insert_head()

# container ------

  expl_os <- list()

# analysis data ------

  insert_msg('Analysis data')

  expl_os$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic[c('sample_id', 'os_days', 'death')]) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = unname(globals$cohort_labs[cohort]),
           cohort = factor(cohort, unname(globals$cohort_labs))) %>%
    filter(complete.cases(.))

# Numbers of cases and events ------

  insert_msg('Numbers of observations and events')

  expl_os$n_numbers <- expl_os$data %>%
    blast(cohort) %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x$death))) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, levels(expl_os$data$cohort)),
           legend_lab = paste0(cohort, '\ntotal: n = ',
                               n_total, '\nevents: n = ',
                               n_events))

# Surv-fit object ------

  insert_msg('Survfit objects')

  expl_os$survfit_obj <-
    surv_fit(Surv(os_days, death) ~ cohort,
             data = expl_os$data)

# Median survival times ------

  insert_msg('Median survival times')

  expl_os$stats <- expl_os$survfit_obj %>%
    surv_median %>%
    mutate(cohort = stri_replace(strata, regex = '.*=', replacement = '')) %>%
    as_tibble

# Test for differences between the cohorts ------

  insert_msg('Peto-Peto test')

  expl_os$test <- expl_os$survfit_obj %>%
    surv_pvalue(method = 'S1') %>%
    re_adjust('pval', method = 'none') %>%
    as_tibble

# Kaplan-Meier plot --------

  insert_msg('Kaplan-Meier plot')

  expl_os$plot <- expl_os$survfit_obj %>%
    ggsurvplot(pval = expl_os$test$significance,
               pval.size = 2.75)

  expl_os$plot$plot <- expl_os$plot$plot +
    scale_color_manual(values = unname(globals$cohort_colors),
                       labels = expl_os$n_numbers$legend_lab,
                       name = '') +
    globals$common_theme +
    labs(title = 'Overall survival',
         x = 'overall survival, days')

# Result table -------

  insert_msg('Result table')

  expl_os$result_tbl <- expl_os$stats %>%
    mutate(median = paste0(signif(median, 3),
                           ' [95% CI: ',
                           signif(lower, 3),
                           ' to ', signif(upper, 3), ']')) %>%
    select(cohort, median) %>%
    column_to_rownames('cohort') %>%
    t %>%
    as_tibble

  expl_os$result_tbl <- expl_os$result_tbl %>%
    mutate(variable = 'Median overall survival',
           p_value = expl_os$test$significance) %>%
    relocate(variable)

# END ------

  expl_os$data <- NULL
  expl_os$survfit_obj <- NULL

  expl_os <- compact(expl_os)

  insert_tail()
