# Differences in disease-specific survival and relapse-free survival
# in TCGA patients

  insert_head()

# container ------

  fgfr_rfs <- list()

# analysis data --------

  insert_msg('Analysis data')

  fgfr_rfs$data <- list(tss = tcga, rfs = tcga) %>%
    map(~.x$clinic)

  fgfr_rfs$data$tss <-
    fgfr_rfs$data$tss[c('sample_id', 'tss_days', 'tumor_death', 'pt_stage')]

  fgfr_rfs$data$rfs <-
    fgfr_rfs$data$rfs[c('sample_id', 'rfs_days', 'relapse', 'pt_stage')]

  fgfr_rfs$data <- fgfr_rfs$data %>%
    map(inner_join,
        fgfr_globals$genetic_groups$tcga,
        by = 'sample_id')

  ## variants of the data set without T4 cancers

  fgfr_rfs$data[c('tss_wo_t4', 'rfs_wo_t4')] <-
    fgfr_rfs$data[c("tss", "rfs")] %>%
    map(filter, pt_stage != 'T4')

  fgfr_rfs$data <- fgfr_rfs$data %>%
    map(select, -pt_stage) %>%
    map(~filter(.x, complete.cases(.x)))

# N numbers -------

  insert_msg('N numbers')

  ## complete observations and events

  fgfr_rfs$n_numbers <- fgfr_rfs$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x[[3]])))

  fgfr_rfs$n_captions <- fgfr_rfs$n_numbers %>%
    map_chr(~paste0('total: n = ', .x[[1]], ', events: n = ', .x[[2]]))

  ## numbers of cases in the genetic groups

  fgfr_rfs$n_groups <- fgfr_rfs$data %>%
    map(blast, FGFR3_mutation) %>%
    map(map,
        ~tibble(n_total = nrow(.x),
                n_events = sum(.x[[3]]))) %>%
    map(compress,
        names_to = 'FGFR3_mutation') %>%
    map(mutate,
        legend_lab = paste0(FGFR3_mutation,
                            '\ntotal: n = ', n_total,
                            '\nevents: n = ', n_events))

  fgfr_rfs$legend_labs <- fgfr_rfs$n_groups %>%
    map(~set_names(.x$legend_lab, .x$FGFR3_mutation))

# Survfit objects --------

  insert_msg('Survfit objects')

  fgfr_rfs$survfit_objects <-
    list(formula = list(tss = Surv(tss_days, tumor_death) ~ FGFR3_mutation,
                        rfs = Surv(rfs_days, relapse) ~ FGFR3_mutation,
                        tss_wo_t4 = Surv(tss_days, tumor_death) ~ FGFR3_mutation,
                        rfs_wo_t4 = Surv(rfs_days, relapse) ~ FGFR3_mutation),
         data = fgfr_rfs$data) %>%
    pmap(surv_fit)

# Median survival --------

  insert_msg('Median survival')

  fgfr_rfs$stats <- fgfr_rfs$survfit_objects %>%
    map(surv_median) %>%
    map(mutate,
        FGFR3_mutation = stri_replace(strata,
                                     regex = '.*=',
                                     replacement = '')) %>%
    compress(names_to = 'response') %>%
    as_tibble

# Peto-Peto test for general differences in survival -----

  insert_msg('Testing')

  fgfr_rfs$test <- fgfr_rfs$survfit_objects %>%
    map(surv_pvalue, method = 'S1') %>%
    map(re_adjust, 'pval', method = 'BH') %>%
    compress(names_to = 'response') %>%
    as_tibble

# Plotting -------

  insert_msg('Kaplan-Meier plots')

  ## bare plots

  fgfr_rfs$km_plots <-
    list(fit = fgfr_rfs$survfit_objects,
         pval = fgfr_rfs$test$significance) %>%
    pmap(ggsurvplot,
         pval.size = 2.75)

  ## additional styling

  fgfr_rfs$km_plots <-
    list(x = fgfr_rfs$km_plots,
         y = c(paste('Disease-specific survival,',
                     globals$cohort_labs["tcga"]),
               paste('Relapse-free survival,',
                     globals$cohort_labs["tcga"]),
               paste('Disease-specific survival,',
                     globals$cohort_labs["tcga"],
                     'without pT4 cancers'),
               paste('Relapse-free survival,',
                     globals$cohort_labs["tcga"],
                     'without pT4 cancers')),
         w = fgfr_rfs$n_captions,
         z = fgfr_rfs$legend_labs,
         v = rep(c('disease specific survival, days',
                   'relapse-free survival, days'),
                 2)) %>%
    pmap(function(x, y, w, z, v) x$plot +
           labs(title = y,
                subtitle = w,
                x = v) +
           scale_color_manual(values = unname(globals$mut_status_color),
                              labels = unname(z),
                              name = '') +
           globals$common_theme)

# Result table -------

  insert_msg('Result table')

  fgfr_rfs$result_tbl <- fgfr_rfs$stats %>%
    mutate(median = paste0(signif(median, 3),
                           ' [95% CI: ',
                           signif(lower, 3),
                           ' to ', signif(upper, 3), ']')) %>%
    blast(response) %>%
    map(select, FGFR3_mutation, median) %>%
    map(column_to_rownames, 'FGFR3_mutation') %>%
    map(t) %>%
    map(as_tibble) %>%
    compress(names_to = 'response')

  fgfr_rfs$result_tbl <- fgfr_rfs$result_tbl %>%
    mutate(response = car::recode(response,
                                  "'rfs' = 'relapse-free survival';
                                  'os' = 'overall survival'"),
           p_value = fgfr_rfs$test$significance) %>%
    relocate(response)

# END -------

  fgfr_rfs$data <- NULL
  fgfr_rfs$survfit_objects <- NULL

  fgfr_rfs <- compact(fgfr_rfs)

  insert_tail()
