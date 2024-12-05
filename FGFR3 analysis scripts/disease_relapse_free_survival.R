# Differences in disease-specific survival and relapse-free survival
# in TCGA patients

  insert_head()

# container ------

  fgfr_rfs <- list()

# analysis data --------

  insert_msg('Analysis data')

  fgfr_rfs$data <- list(tss_tcga = tcga,
                        rfs_tcga = tcga,
                        tss_bcan = bcan) %>%
    map(~.x$clinic)

  fgfr_rfs$data[c("tss_tcga", "tss_bcan")] <-
    fgfr_rfs$data[c("tss_tcga", "tss_bcan")] %>%
    map(~.x[c('sample_id',
              'tss_days',
              'tumor_death',
              'pt_stage')]) %>%
    map(set_names,
        c('sample_id', 'time', 'event', 'pt_stage'))

  fgfr_rfs$data$rfs_tcga <-
    fgfr_rfs$data$rfs[c('sample_id',
                        'rfs_days',
                        'relapse',
                        'pt_stage')] %>%
    set_names(c('sample_id', 'time', 'event', 'pt_stage'))

  fgfr_rfs$data <-
    map2(fgfr_rfs$data,
         fgfr_globals$genetic_groups[c("tcga", "tcga", "bcan")],
         inner_join, by = 'sample_id')

  ## variants of the data set without T4 cancers

  fgfr_rfs$data[c('tss_tcga_wo_t4', 'rfs_tcga_wo_t4', 'tss_bcan_wo_t4')] <-
    fgfr_rfs$data %>%
    map(filter, pt_stage != 'T4')

  fgfr_rfs$data <- fgfr_rfs$data %>%
    map(select, -pt_stage) %>%
    map(~filter(.x, complete.cases(.x)))

# N numbers -------

  insert_msg('N numbers')

  ## complete observations and events

  fgfr_rfs$n_numbers <- fgfr_rfs$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x[['event']])))

  fgfr_rfs$n_captions <- fgfr_rfs$n_numbers %>%
    map_chr(~paste0('total: n = ', .x[[1]],
                    ', events: n = ', .x[[2]]))

  ## numbers of cases in the genetic groups

  fgfr_rfs$n_groups <- fgfr_rfs$data %>%
    map(blast, FGFR3_mutation) %>%
    map(map,
        ~tibble(n_total = nrow(.x),
                n_events = sum(.x[['event']]))) %>%
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

  fgfr_rfs$survfit_objects <- fgfr_rfs$data %>%
    map(surv_fit,
        formula = Surv(time, event) ~ FGFR3_mutation)

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
                     globals$cohort_labs["bcan"]),
               paste('Disease-specific survival,',
                     globals$cohort_labs["tcga"],
                     'without pT4 cancers'),
               paste('Relapse-free survival,',
                     globals$cohort_labs["tcga"],
                     'without pT4 cancers'),
               paste('Disease-specific survival,',
                     globals$cohort_labs["bcan"],
                     'without pT4 cancers')),
         w = fgfr_rfs$n_captions,
         z = fgfr_rfs$legend_labs,
         v = rep(c('disease specific survival, days',
                   'relapse-free survival, days',
                   'disease specific survival, days'),
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

  ## numbers of patients and events

  fgfr_rfs$result_tbl$n_numbers <- fgfr_rfs$n_groups %>%
    map(column_to_rownames, 'FGFR3_mutation') %>%
    map(select, n_total, n_events) %>%
    map(t) %>%
    map(as_tibble) %>%
    map(mutate, variable = c('Patients, N', 'Events, N')) %>%
    map(relocate, variable)

  ## median survival times and test results

  fgfr_rfs$result_tbl$stats <- fgfr_rfs$stats %>%
    blast(response) %>%
    map(mutate,
        median = ifelse(is.na(median),
                        '',
                        paste0(signif(median, 3),
                               ' [95% CI: ',
                               signif(lower, 3),
                               ' to ', signif(upper, 3), ']'))) %>%
    map(select, FGFR3_mutation, median) %>%
    map(column_to_rownames, 'FGFR3_mutation') %>%
    map(t) %>%
    map(as_tibble)

  fgfr_rfs$result_tbl$stats <-
    map2(fgfr_rfs$result_tbl$stats,
         fgfr_rfs$test %>%
           blast(response) %>%
           map(~.x['significance']),
         cbind) %>%
    map(mutate, variable = 'median survival, days') %>%
    map(as_tibble)

  ## the entire result table

  fgfr_rfs$result_tbl <-
    map2(fgfr_rfs$result_tbl$n_numbers,
         fgfr_rfs$result_tbl$stats,
         full_rbind) %>%
    compress(names_to = 'response_cohort') %>%
    mutate(response = stri_split_fixed(response_cohort,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
           response = car::recode(response,
                                  "'tss' = 'disease-specific survival';
                                  'rfs' = 'relapse-specific survival'"),
           cohort = stri_split_fixed(response_cohort,
                                     pattern = '_',
                                     simplify = TRUE)[, 2],
           cohort = globals$cohort_labs[cohort],
           subset = ifelse(stri_detect(response_cohort,
                                       regex = 't4$'),
                           'without pT4 cancers', 'all patients')) %>%
    relocate(response, cohort, subset) %>%
    select(-response_cohort) %>%
    set_names(c('Response', 'Cohort', 'Subset', 'Variable',
                levels(fgfr_rfs$data[[1]]$FGFR3_mutation),
                'Significance'))

# END -------

  fgfr_rfs$data <- NULL
  fgfr_rfs$survfit_objects <- NULL

  fgfr_rfs <- compact(fgfr_rfs)

  insert_tail()
