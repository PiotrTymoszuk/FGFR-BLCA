# demographic and clinical characteristic of the consensus subsets of MIBC

  insert_head()

# container --------

  lca_clinic <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis variables and data --------

  insert_msg('Analysis variables and data')

  ## variables: removal of the ones analyzed by other scripts

  lca_clinic$lexicon <- globals$clinic_lexicon %>%
    filter(!variable %in% c('tmb_per_mb', 'bi_reponse', 'death')) %>%
    mutate(test_type = ifelse(format == 'numeric',
                              'kruskal_etasq', 'cramer_v'),
           plot_type = ifelse(format == 'numeric',
                              'box', 'stack'),
           ax_label = ifelse(is.na(unit),
                             '% of subset', unit))

  ## analysis data

  lca_clinic$data <-
    list(tcga = tcga, imvigor = imvigor, bcan = bcan) %>%
    map(~.x$clinic) %>%
    map(select, sample_id, any_of(lca_clinic$lexicon$variable))

  lca_clinic$data <-
    map2(subtypes$assignment[names(lca_clinic$data)],
         lca_clinic$data,
         inner_join, by = 'sample_id') %>%
    map(filter, !is.na(consensusClass)) %>%
    map(map_dfc, function(x) if(is.factor(x)) droplevels(x) else x)

  ## data set-specific variable vectors. No point at testing
  ## the tissue variable for the TCGA BLCA data set (bladder only!)

  lca_clinic$variables <- lca_clinic$data %>%
    map(names) %>%
    map(~filter(lca_clinic$lexicon, variable %in% .x)) %>%
    map(~.x$variable)

  lca_clinic$variables$tcga <-
    lca_clinic$variables$tcga[lca_clinic$variables$tcga != 'tissue']

# N numbers --------

  insert_msg('N numbers of samples')

  lca_clinic$n_numbers <- lca_clinic$data %>%
    map(count, consensusClass) %>%
    map(column_to_rownames, 'consensusClass') %>%
    map(t) %>%
    map(as_tibble) %>%
    map(mutate,
        variable = 'Samples, N') %>%
    map(relocate, variable)

# Descriptive stats -------

  insert_msg('Analysis stats')

  lca_clinic$stats <-
    future_map2(lca_clinic$data,
                lca_clinic$variables,
                ~explore(.x,
                         variables = .y,
                         split_factor = 'consensusClass',
                         what = 'table',
                         pub_styled = TRUE),
                .options = furrr_options(seed = TRUE)) %>%
    map(format_desc)

# Testing for differences between the genetic clusters ------

  insert_msg('Testing for differences between the genetic clusters')

  lca_clinic$test <-
    future_map2(lca_clinic$data,
                lca_clinic$variables,
                ~compare_variables(.x,
                                   variables = .y,
                                   split_factor = 'consensusClass',
                                   what = 'eff_size',
                                   types = exchange(.y,
                                                    lca_clinic$lexicon,
                                                    value = 'test_type'),
                                   exact = FALSE,
                                   ci = FALSE,
                                   pub_styled = TRUE,
                                   adj_method = 'BH'),
                .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

# Plots --------

  insert_msg('Plots')

  for(i in names(lca_clinic$data)) {

    lca_clinic$plots[[i]] <-
      list(variable = lca_clinic$test[[i]]$variable,
           plot_title = lca_clinic$test[[i]]$variable %>%
             exchange(lca_clinic$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = lca_clinic$test[[i]]$plot_cap,
           y_lab = lca_clinic$test[[i]]$variable %>%
             exchange(lca_clinic$lexicon,
                      value = 'ax_label'),
           type = lca_clinic$test[[i]]$variable %>%
             exchange(lca_clinic$lexicon,
                      value = 'plot_type')) %>%
      future_pmap(plot_variable,
                  lca_clinic$data[[i]],
                  split_factor = 'consensusClass',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  x_n_labs = TRUE,
                  .options = furrr_options(seed = TRUE)) %>%
      set_names(lca_clinic$test[[i]]$variable)

    lca_clinic$plots[[i]]$age <- lca_clinic$plots[[i]]$age +
      scale_fill_manual(values = globals$genet_colors)

    lca_clinic$plots[[i]][names(lca_clinic$plots[[i]]) != 'age'] <-
      lca_clinic$plots[[i]][names(lca_clinic$plots[[i]]) != 'age'] %>%
      map(~.x + scale_fill_brewer(palette = 'Reds'))

  }

# Result table --------

  insert_msg('Result table')

  lca_clinic$result_tbl <-
    map2(lca_clinic$stats,
         lca_clinic$test %>%
           map(filter, !stri_detect(eff_size, fixed = 'Inf')),
         merge_stat_test) %>%
    map(format_res_tbl,
        dict = lca_clinic$lexicon) %>%
    map2(., lca_clinic$n_numbers,
         ~full_rbind(.y, .x)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(lca_clinic$data[[1]]$consensusClass),
                'Significance',
                'Effect size'))

# END -------

  plan('sequential')

  rm(i)

  lca_clinic$data <- NULL

  lca_clinic <- compact(lca_clinic)

  insert_tail()
