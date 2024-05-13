# Differences in overall survival in patients with and without FGFR3 mutations.

  insert_head()

# container --------

  fgfr_os <- list()

# analysis data --------

  insert_msg('Analysis data')

  fgfr_os$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic) %>%
    map(select, sample_id, os_days, death, any_of('pt_stage'))

  fgfr_os$data <- map2(fgfr_globals$genetic_groups[names(fgfr_os$data)],
                       fgfr_os$data,
                       inner_join, by = 'sample_id')

  ## TCGA: variant of the data set without pT4 cancers

  fgfr_os$data$tcga_wo_t4 <- fgfr_os$data$tcga %>%
    filter(pt_stage != 'T4')

  fgfr_os$data <- fgfr_os$data %>%
    map(~.x[names(.x) != 'pt_stage']) %>%
    map(~filter(.x, complete.cases(.x)))

# N numbers -------

  insert_msg('N numbers')

  ## complete observations and events

  fgfr_os$n_numbers <- fgfr_os$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x$death)))

  fgfr_os$n_captions <- fgfr_os$n_numbers %>%
    map_chr(~paste0('total: n = ', .x[[1]], ', events: n = ', .x[[2]]))

  ## numbers of cases in the genetic groups

  fgfr_os$n_groups <- fgfr_os$data %>%
    map(blast, FGFR3_mutation) %>%
    map(map,
        ~tibble(n_total = nrow(.x),
                n_events = sum(.x$death))) %>%
    map(compress,
        names_to = 'FGFR3_mutation') %>%
    map(mutate,
        legend_lab = paste0(FGFR3_mutation,
                            '\ntotal: n = ', n_total,
                            '\nevents: n = ', n_events))

  fgfr_os$legend_labs <- fgfr_os$n_groups %>%
    map(~set_names(.x$legend_lab, .x$FGFR3_mutation))

# Survfit objects --------

  insert_msg('Survfit objects')

  fgfr_os$survfit_objects <- fgfr_os$data %>%
    map(surv_fit,
        formula = Surv(os_days, death) ~ FGFR3_mutation)

# Median survival --------

  insert_msg('Median survival')

  fgfr_os$stats <- fgfr_os$survfit_objects %>%
    map(surv_median) %>%
    map(mutate,
        FGFR3_mutation = stri_replace(strata,
                                     regex = '.*=',
                                     replacement = '')) %>%
    map(as_tibble)

# Peto-Peto test for general differences in survival -----

  insert_msg('Testing')

  fgfr_os$test <- fgfr_os$survfit_objects %>%
    map(surv_pvalue, method = 'S1') %>%
    map(re_adjust, 'pval', method = 'BH') %>%
    map(as_tibble)

# Plotting -------

  insert_msg('Kaplan-Meier plots')

  ## bare plots

  fgfr_os$km_plots <-
    list(fit = fgfr_os$survfit_objects,
         pval = map(fgfr_os$test, ~.x$significance)) %>%
    pmap(ggsurvplot,
         pval.size = 2.75)

  ## additional styling

  fgfr_os$km_plots <-
    list(x = fgfr_os$km_plots,
         y = globals$cohort_labs[names(fgfr_os$km_plots)] %>%
           ifelse(is.na(.),
                  paste(globals$cohort_labs["tcga"], 'without pT4 cancers'),
                  .),
         w = fgfr_os$n_captions,
         z = fgfr_os$legend_labs) %>%
    pmap(function(x, y, w, z) x$plot +
           labs(title = y,
                subtitle = w,
                x = 'overall survival, days') +
           scale_color_manual(values = unname(globals$mut_status_color),
                              labels = unname(z),
                              name = html_italic('FGFR3')) +
           globals$common_theme +
           theme(legend.title = element_markdown()))

# Result table --------

  insert_msg('Result table')

  fgfr_os$result_tbl <- fgfr_os$stats %>%
    map(mutate,
        median = paste0(signif(median, 3),
                           ' [95% CI: ',
                           signif(lower, 3),
                           ' to ', signif(upper, 3), ']')) %>%
    map(select, FGFR3_mutation, median) %>%
    map(column_to_rownames, 'FGFR3_mutation') %>%
    map(t) %>%
    map(as_tibble)

  fgfr_os$result_tbl <-
    map2(fgfr_os$result_tbl,
         map(fgfr_os$test, ~.x[c('variable', 'significance')]),
         cbind) %>%
    map(as_tibble) %>%
    map(relocate, variable)

# END ------

  fgfr_os$data <- NULL
  fgfr_os$survfit_objects <- NULL

  fgfr_os <- compact(fgfr_os)

  insert_tail()
