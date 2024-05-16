# differences in total mutation burden between the cancers stratified
# as specified in 'globals.R'

  insert_head()

# container -------

  fgfr_tmb <- list()

# analysis data --------

  insert_msg('Analysis data')

  ## TMB provided by the study authors

  fgfr_tmb$data$tmb <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic) %>%
    map(select, sample_id, tmb_per_mb, any_of('msi_score'))

  fgfr_tmb$data$tmb <-
    map2(fgfr_globals$genetic_groups[names(fgfr_tmb$data$tmb)],
         fgfr_tmb$data$tmb,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

  ## counts of mutations, amplifications and deletions
  ## Amplifications and deletions in the IMvigor cohort are copies
  ## of the mutation counts

  fgfr_tmb$data$counts <-  globals$cohort_expr %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    map(compact) %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map(map, ~.x[names(.x) != 'patient_id']) %>%
    map(map, rowSums) %>%
    map(~map2(.x, c('mutations', 'deletions', 'amplifications'),
              ~compress(.x,
                        names_to = 'sample_id',
                        values_to = .y))) %>%
    map(reduce, left_join, by = 'sample_id')

  fgfr_tmb$data <-
    map2(fgfr_tmb$data$tmb,
         fgfr_tmb$data$counts,
         left_join, by = 'sample_id')

  ## variable lexicon

  fgfr_tmb$lexicon <-
    tibble(variable = c('tmb_per_mb',
                        'mutations',
                        'deletions',
                        'amplifications',
                        'msi_score'),
           label = c('TMB',
                     'Mutations',
                     'Deletions',
                     'Amplifications',
                     'MSI score'),
           ax_label = c('mutations/MB',
                        rep('count', 3),
                        'score'))

  ## data set-specific variable vectors

  fgfr_tmb$variables <- fgfr_tmb$data %>%
    map(names) %>%
    map(~filter(fgfr_tmb$lexicon,
                variable %in% .x)) %>%
    map(~.x$variable)

# Descriptive stats ------

  insert_msg('Descriptive stats')

  fgfr_tmb$stats <-
    map2(fgfr_tmb$data,
         fgfr_tmb$variables,
         ~explore(.x,
                  variables = .y,
                  split_factor = 'FGFR3_mutation',
                  what = 'table',
                  pub_styled = TRUE)) %>%
    map(format_desc)

# Testing for differences between the genetic strata -------

  insert_msg('Testing for differences between the genetic strata')

  fgfr_tmb$test <-
    map2(fgfr_tmb$data,
         fgfr_tmb$variables,
         ~compare_variables(.x,
                            variables = .y,
                            split_factor = 'FGFR3_mutation',
                            what = 'eff_size',
                            types = 'wilcoxon_r',
                            ci = FALSE,
                            exact = FALSE,
                            pub_styled = TRUE,
                            adj_method = 'BH')) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

# Plotting ------

  insert_msg('Box plots')

  for(i in names(fgfr_tmb$data)) {

    fgfr_tmb$plots[[i]] <-
      list(variable = fgfr_tmb$test[[i]]$variable,
           plot_title = fgfr_tmb$test[[i]]$variable %>%
             exchange(fgfr_tmb$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = fgfr_tmb$test[[i]]$plot_cap,
           y_lab = fgfr_tmb$test[[i]]$variable %>%
             exchange(fgfr_tmb$lexicon,
                      value = 'ax_label')) %>%
      pmap(plot_variable,
           fgfr_tmb$data[[i]],
           split_factor = 'FGFR3_mutation',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(axis.title.x = element_markdown(),
                   legend.position = 'none'),
           point_hjitter = 0,
           x_n_labs = TRUE,
           x_lab = html_italic('FGFR3')) %>%
      map(~.x +
            scale_fill_manual(values = globals$mut_status_color)) %>%
      set_names(fgfr_tmb$test[[i]]$variable)

  }

  fgfr_tmb$plots$genie$tmb_per_mb <-
    fgfr_tmb$plots$genie$tmb_per_mb +
    labs(y = 'fraction of genome')

# Result tables -------

  insert_msg('Result table')

  fgfr_tmb$result_tbl <-
    map2(fgfr_tmb$stats,
         fgfr_tmb$test,
         merge_stat_test) %>%
    map(format_res_tbl,
        dict = fgfr_tmb$lexicon)

# END -----

  insert_tail()
