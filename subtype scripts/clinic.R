# demographic and clinical characteristic of the consensus subsets of MIBC

  insert_head()

# container --------

  sub_clinic <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis variables and data --------

  insert_msg('Analysis variables and data')

  ## variables: removal of the ones analyzed by other scripts

  sub_clinic$lexicon <- globals$clinic_lexicon %>%
    filter(!variable %in% c('tmb_per_mb', 'bi_reponse', 'death',
                            'tissue', 'invasiveness')) %>%
    mutate(test_type = ifelse(format == 'numeric',
                              'kruskal_etasq', 'cramer_v'),
           plot_type = ifelse(format == 'numeric',
                              'box', 'stack'),
           table_label = ifelse(!is.na(unit),
                                paste(label, unit, sep = ', '),
                                label),
           ax_label = ifelse(is.na(unit),
                             '% of subset', unit))

  ## analysis data

  sub_clinic$data <-
    list(tcga = tcga, imvigor = imvigor, bcan = bcan) %>%
    map(~.x$clinic) %>%
    map(select, sample_id, any_of(sub_clinic$lexicon$variable))

  sub_clinic$data <-
    map2(subtypes$assignment[names(sub_clinic$data)],
         sub_clinic$data,
         inner_join, by = 'sample_id') %>%
    map(filter, !is.na(consensusClass)) %>%
    map(map_dfc, function(x) if(is.factor(x)) droplevels(x) else x)

  ## consensus classes to test: only the LumP, LumU, Stroma-rich, and Ba/Sq
  ## classes

  sub_clinic$data <- sub_clinic$data %>%
    map(filter,
        !consensusClass %in% c('not assigned', 'NE-like', 'LumNS')) %>%
    map(mutate,
        consensusClass = droplevels(consensusClass))

  ## data set-specific variable vectors. No point at testing
  ## the tissue variable (bladder only)

  sub_clinic$variables <- sub_clinic$data %>%
    map(names) %>%
    map(~filter(sub_clinic$lexicon, variable %in% .x)) %>%
    map(~.x$variable)

  sub_clinic$variables$tcga <-
    sub_clinic$variables$tcga[sub_clinic$variables$tcga != 'tissue']

  ## characters to vectors

  sub_clinic$data <- sub_clinic$data %>%
    map(select, -sample_id) %>%
    map(map_dfc, function(x) if(is.character(x)) factor(x) else x)

# N numbers --------

  insert_msg('N numbers of samples')

  sub_clinic$n_numbers <- sub_clinic$data %>%
    map(count, consensusClass) %>%
    map(column_to_rownames, 'consensusClass') %>%
    map(t) %>%
    map(as_tibble) %>%
    map(mutate,
        variable = 'Samples, N') %>%
    map(relocate, variable)

# Descriptive stats -------

  insert_msg('Analysis stats')

  sub_clinic$stats <-
    future_map2(sub_clinic$data,
                sub_clinic$variables,
                ~explore(.x,
                         variables = .y,
                         split_factor = 'consensusClass',
                         what = 'table',
                         pub_styled = TRUE),
                .options = furrr_options(seed = TRUE)) %>%
    map(format_desc)

# Testing for differences between the consensus classes ------

  insert_msg('Testing for differences between the consensus classes')

  sub_clinic$test <-
    future_map2(sub_clinic$data,
                sub_clinic$variables,
                ~compare_variables(.x,
                                   variables = .y,
                                   split_factor = 'consensusClass',
                                   what = 'eff_size',
                                   types = exchange(.y,
                                                    sub_clinic$lexicon,
                                                    value = 'test_type'),
                                   exact = FALSE,
                                   ci = FALSE,
                                   pub_styled = TRUE,
                                   adj_method = 'BH'),
                .options = furrr_options(seed = TRUE)) %>%
    map(p_formatter, text = TRUE) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '))

# Plots --------

  insert_msg('Plots')

  for(i in names(sub_clinic$data)) {

    sub_clinic$plots[[i]] <-
      list(variable = sub_clinic$test[[i]]$variable,
           plot_title = sub_clinic$test[[i]]$variable %>%
             exchange(sub_clinic$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = sub_clinic$test[[i]]$plot_cap,
           y_lab = sub_clinic$test[[i]]$variable %>%
             exchange(sub_clinic$lexicon,
                      value = 'ax_label'),
           type = sub_clinic$test[[i]]$variable %>%
             exchange(sub_clinic$lexicon,
                      value = 'plot_type')) %>%
      future_pmap(plot_variable,
                  sub_clinic$data[[i]],
                  split_factor = 'consensusClass',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  x_n_labs = TRUE,
                  .options = furrr_options(seed = TRUE)) %>%
      set_names(sub_clinic$test[[i]]$variable)

    sub_clinic$plots[[i]]$age <- sub_clinic$plots[[i]]$age +
      scale_fill_manual(values = globals$sub_colors)

    sub_clinic$plots[[i]][names(sub_clinic$plots[[i]]) != 'age'] <-
      sub_clinic$plots[[i]][names(sub_clinic$plots[[i]]) != 'age'] %>%
      map(~.x + scale_fill_brewer(palette = 'Reds'))

  }

# Result table --------

  insert_msg('Result table')

  sub_clinic$result_tbl <-
    map2(sub_clinic$stats,
         sub_clinic$test %>%
           map(filter, !stri_detect(eff_size, fixed = 'Inf')),
         merge_stat_test) %>%
    map(format_res_tbl,
        dict = sub_clinic$lexicon,
        value = "table_label") %>%
    map2(., sub_clinic$n_numbers,
         ~full_rbind(.y, .x)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(sub_clinic$data[[1]]$consensusClass),
                'FDR p value',
                'Effect size'))

# END -------

  plan('sequential')

  rm(i)

  sub_clinic$data <- NULL

  sub_clinic <- compact(sub_clinic)

  insert_tail()
