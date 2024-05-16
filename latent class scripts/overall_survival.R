# Differences in overall survival in the genetic clusters.

  insert_head()

# container --------

  lca_os <- list()

# analysis data --------

  insert_msg('Analysis data')

  ## analysis data: GENIE and TCGA

  lca_os$data <- list(genie = genie,
                      msk = msk,
                      tcga = tcga) %>%
    map(~.x$clinic) %>%
    map(select, sample_id, os_days, death, any_of('pt_stage'))

  lca_os$data <- map2(lca_pred$assignment[names(lca_os$data)],
                      lca_os$data,
                      inner_join, by = 'sample_id')

  ## MSK and TCGA: variants without T4 tumors

  lca_os$data[c('msk_wo_t4', 'tcga_wo_t4')] <-
    lca_os$data[c('msk', 'tcga')] %>%
    map(filter, pt_stage != 'T4')

  lca_os$data <- lca_os$data %>%
    map(~.x[names(.x) != 'pt_stage']) %>%
    map(~filter(.x, complete.cases(.x)))

# N numbers -------

  insert_msg('N numbers')

  ## complete observations and events

  lca_os$n_numbers <- lca_os$data %>%
    map(~tibble(n_total = nrow(.x),
                n_events = sum(.x$death)))

  lca_os$n_captions <- lca_os$n_numbers %>%
    map_chr(~paste0('total: n = ', .x[[1]], ', events: n = ', .x[[2]]))

  ## numbers of cases in the genetic groups

  lca_os$n_groups <- lca_os$data %>%
    map(blast, clust_id) %>%
    map(map,
        ~tibble(n_total = nrow(.x),
                n_events = sum(.x$death))) %>%
    map(compress,
        names_to = 'clust_id') %>%
    map(mutate,
        legend_lab = paste0(clust_id,
                            '\ntotal: n = ', n_total,
                            '\nevents: n = ', n_events))

  lca_os$legend_labs <- lca_os$n_groups %>%
    map(~set_names(.x$legend_lab, .x$clust_id))

# Survfit objects --------

  insert_msg('Survfit objects')

  lca_os$survfit_objects <- lca_os$data %>%
    map(surv_fit,
        formula = Surv(os_days, death) ~ clust_id)

# Median survival --------

  insert_msg('Median survival')

  lca_os$stats <- lca_os$survfit_objects %>%
    map(surv_median) %>%
    map(mutate,
        clust_id = stri_replace(strata,
                                regex = '.*=',
                                replacement = '')) %>%
    map(as_tibble)

# Peto-Peto test for general differences in survival -----

  insert_msg('Testing')

  lca_os$test <- lca_os$survfit_objects %>%
    map(surv_pvalue, method = 'S1') %>%
    map(re_adjust, 'pval', method = 'BH') %>%
    map(as_tibble)

# Plotting -------

  insert_msg('Kaplan-Meier plots')

  ## bare plots

  lca_os$km_plots <-
    list(fit = lca_os$survfit_objects,
         pval = map(lca_os$test, ~.x$significance)) %>%
    pmap(ggsurvplot,
         pval.size = 2.75)

  ## additional styling

  lca_os$km_plots <-
    list(x = lca_os$km_plots,
         y = c(globals$cohort_labs[c("genie", "msk", "tcga")],
               paste(globals$cohort_labs[c("msk", "tcga")],
                     ' without pT4 cancers')),
         w = lca_os$n_captions,
         z = lca_os$legend_labs) %>%
    pmap(function(x, y, w, z) x$plot +
           labs(title = y,
                subtitle = w,
                x = 'overall survival, days') +
           scale_color_manual(values = unname(globals$genet_colors),
                              labels = unname(z),
                              name = '') +
           globals$common_theme)

# Result table --------

  insert_msg('Result table')

  lca_os$result_tbl <- lca_os$stats %>%
    map(mutate,
        median = paste0(signif(median, 3),
                           ' [95% CI: ',
                           signif(lower, 3),
                           ' to ', signif(upper, 3), ']')) %>%
    map(select, clust_id, median) %>%
    map(column_to_rownames, 'clust_id') %>%
    map(t) %>%
    map(as_tibble)

  lca_os$result_tbl <-
    map2(lca_os$result_tbl,
         map(lca_os$test, ~.x[c('variable', 'significance')]),
         cbind) %>%
    map(as_tibble) %>%
    map(relocate, variable)

# END ------

  lca_os$data <- NULL
  lca_os$survfit_objects <- NULL

  lca_os <- compact(lca_os)

  insert_tail()
