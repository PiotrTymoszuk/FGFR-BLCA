# Differences in disease-specific survival and relapse-free survival between
# the genetic clusters in the TCGA cohort

  insert_head()

# container ------

  lca_rfs <- list()

# analysis data --------

  insert_msg('Analysis data')

  lca_rfs$data <- list(tss = tcga, rfs = tcga) %>%
    map(~.x$clinic)

  lca_rfs$data$tss <-
    lca_rfs$data$tss[c('sample_id', 'tss_days', 'tumor_death', 'pt_stage')]

  lca_rfs$data$rfs <-
    lca_rfs$data$rfs[c('sample_id', 'rfs_days', 'relapse', 'pt_stage')]

  lca_rfs$data <- lca_rfs$data %>%
    map(inner_join,
        lca_pred$assignment$tcga,
        by = 'sample_id')

  ## data set variants without T4 patients

  lca_rfs$data[c('tss_wo_t4',
                 'rfs_wo_t4')] <-
    lca_rfs$data[c("tss",
                   "rfs")] %>%
    map(filter, pt_stage != 'T4')

  lca_rfs$data <- lca_rfs$data %>%
    map(select, -pt_stage) %>%
    map(~filter(.x, complete.cases(.x)))

# N numbers -------

  insert_msg('N numbers')

  ## complete observations and events

  lca_rfs$n_numbers <- lca_rfs$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x[[3]])))

  lca_rfs$n_captions <- lca_rfs$n_numbers %>%
    map_chr(~paste0('total: n = ', .x[[1]], ', events: n = ', .x[[2]]))

  ## numbers of cases in the genetic groups

  lca_rfs$n_groups <- lca_rfs$data %>%
    map(blast, clust_id) %>%
    map(map,
        ~tibble(n_total = nrow(.x),
                n_events = sum(.x[[3]]))) %>%
    map(compress,
        names_to = 'clust_id') %>%
    map(mutate,
        legend_lab = paste0(clust_id,
                            '\ntotal: n = ', n_total,
                            '\nevents: n = ', n_events))

  lca_rfs$legend_labs <- lca_rfs$n_groups %>%
    map(~set_names(.x$legend_lab, .x$clust_id))

# Survfit objects --------

  insert_msg('Survfit objects')

  lca_rfs$survfit_objects <-
    list(formula = list(tss = Surv(tss_days, tumor_death) ~ clust_id,
                        rfs = Surv(rfs_days, relapse) ~ clust_id,
                        tss_wo_t4 = Surv(tss_days, tumor_death) ~ clust_id,
                        rfs_wo_t4 = Surv(rfs_days, relapse) ~ clust_id),
         data = lca_rfs$data) %>%
    pmap(surv_fit)

# Median survival --------

  insert_msg('Median survival')

  lca_rfs$stats <- lca_rfs$survfit_objects %>%
    map(surv_median) %>%
    map(mutate,
        clust_id = stri_replace(strata,
                                     regex = '.*=',
                                     replacement = '')) %>%
    compress(names_to = 'response') %>%
    as_tibble

# Peto-Peto test for general differences in survival -----

  insert_msg('Testing')

  lca_rfs$test <- lca_rfs$survfit_objects %>%
    map(surv_pvalue, method = 'S1') %>%
    map(re_adjust, 'pval', method = 'BH') %>%
    compress(names_to = 'response') %>%
    as_tibble

# Plotting -------

  insert_msg('Kaplan-Meier plots')

  ## bare plots

  lca_rfs$km_plots <-
    list(fit = lca_rfs$survfit_objects,
         pval = lca_rfs$test$significance) %>%
    pmap(ggsurvplot,
         pval.size = 2.75)

  ## additional styling

  lca_rfs$km_plots <-
    list(x = lca_rfs$km_plots,
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
         w = lca_rfs$n_captions,
         z = lca_rfs$legend_labs,
         v = rep(c('disease specific survival, days',
                   'relapse-free survival, days'),
                 2)) %>%
    pmap(function(x, y, w, z, v) x$plot +
           labs(title = y,
                subtitle = w,
                x = v) +
           scale_color_manual(values = unname(globals$genet_colors),
                              labels = unname(z),
                              name = '') +
           globals$common_theme)

# Result table -------

  insert_msg('Result table')

  lca_rfs$result_tbl <- lca_rfs$stats %>%
    mutate(median = paste0(signif(median, 3),
                           ' [95% CI: ',
                           signif(lower, 3),
                           ' to ', signif(upper, 3), ']')) %>%
    blast(response) %>%
    map(select, clust_id, median) %>%
    map(column_to_rownames, 'clust_id') %>%
    map(t) %>%
    map(as_tibble) %>%
    compress(names_to = 'response')

  lca_rfs$result_tbl <- lca_rfs$result_tbl %>%
    mutate(response = car::recode(response,
                                  "'rfs' = 'relapse-free survival';
                                  'os' = 'overall survival'"),
           p_value = lca_rfs$test$significance) %>%
    relocate(response)

# END -------

  lca_rfs$data <- NULL
  lca_rfs$survfit_objects <- NULL

  lca_rfs <- compact(lca_rfs)

  insert_tail()
