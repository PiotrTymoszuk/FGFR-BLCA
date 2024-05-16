# Differences in total mutation burdens between the genetic clusters

  insert_head()

# container -------

  lca_tmb <- list()

# analysis data --------

  insert_msg('Analysis data')

  ## TMB determined by the study authors

  lca_tmb$data$tmb <- list(genie = genie,
                           msk = msk,
                           tcga = tcga) %>%
    map(~.x$clinic) %>%
    map(select, sample_id, tmb_per_mb, any_of('msi_score'))

  lca_tmb$data$tmb <-
    map2(lca_pred$assignment[names(lca_tmb$data$tmb)],
         lca_tmb$data$tmb,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

  ## counts of all investigated mutations, amplifications and deletions

  lca_tmb$data$counts <- list(genie = genie,
                              msk = msk,
                              tcga = tcga) %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map(map, rowSums) %>%
    map(~map2(.x, c('mutations', 'deletions', 'amplifications'),
             ~compress(.x,
                       names_to = 'sample_id',
                       values_to = .y))) %>%
    map(reduce, left_join, by = 'sample_id')

  lca_tmb$data <-
    map2(lca_tmb$data$tmb,
         lca_tmb$data$counts,
         left_join, by = 'sample_id')

  ## variable lexicon

  lca_tmb$lexicon <-
    tibble(variable = c('tmb_per_mb',
                        'mutations',
                        'deletions',
                        'amplifications',
                        'msi_score'),
           label = c('TMB',
                     'Mutation count',
                     'Deletion count',
                     'Amplification count',
                     'MSI score'),
           ax_label = c('mutations/MB',
                        rep('count', 3),
                        'score'))

  ## variable vectors specific for particular data sets

  lca_tmb$variables <- lca_tmb$data %>%
    map(names) %>%
    map(~filter(lca_tmb$lexicon, variable %in% .x)) %>%
    map(~.x$variable)

# Descriptive stats ------

  insert_msg('Descriptive stats')

  lca_tmb$stats <-
    map2(lca_tmb$data,
         lca_tmb$variables,
         ~explore(.x,
                  variables = .y,
                  split_factor = 'clust_id',
                  what = 'table',
                  pub_styled = TRUE)) %>%
    map(format_desc)

# Testing for differences between the genetic strata -------

  insert_msg('Testing for differences between the genetic strata')

  lca_tmb$test <-
    map2(lca_tmb$data,
         lca_tmb$variables,
         ~compare_variables(.x,
                            variables = .y,
                            split_factor = 'clust_id',
                            what = 'eff_size',
                            types = 'kruskal_etasq',
                            ci = FALSE,
                            exact = FALSE,
                            pub_styled = TRUE,
                            adj_method = 'BH')) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

# Plotting ------

  insert_msg('Box plots')

  for(i in names(lca_tmb$data)) {

    lca_tmb$plots[[i]] <-
      list(variable = lca_tmb$test[[i]]$variable,
           plot_title = lca_tmb$test[[i]]$variable %>%
             exchange(lca_tmb$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = lca_tmb$test[[i]]$plot_cap,
           y_lab = lca_tmb$test[[i]]$variable %>%
             exchange(lca_tmb$lexicon,
                      value = 'ax_label')) %>%
      pmap(plot_variable,
           lca_tmb$data[[i]],
           split_factor = 'clust_id',
           type = 'box',
           cust_theme = globals$common_theme,
           point_hjitter = 0,
           x_n_labs = TRUE) %>%
      map(~.x +
            scale_fill_manual(values = globals$genet_colors)) %>%
      set_names(lca_tmb$test[[i]]$variable)

  }

  lca_tmb$plots$genie$tmb_per_mb <-
    lca_tmb$plots$genie$tmb_per_mb +
    labs(y = 'fraction of genome')

# Result tables -------

  insert_msg('Result table')

  lca_tmb$result_tbl <-
    map2(lca_tmb$stats,
         lca_tmb$test,
         merge_stat_test) %>%
    map(format_res_tbl,
        dict = lca_tmb$lexicon)

# END -----

  insert_tail()
