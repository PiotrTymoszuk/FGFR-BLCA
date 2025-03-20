# Clinical and pathological characteristic of the FGFR3 mutation strata

  insert_head()

# container ------

  fgfr_clinic <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis variables and data --------

  insert_msg('Analysis variables and data')

  ## variables: removal of the ones analyzed by other scripts

  fgfr_clinic$lexicon <- globals$clinic_lexicon %>%
    filter(!variable %in% c('tmb_per_mb', 'bi_response', 'death')) %>%
    mutate(test_type = ifelse(format == 'numeric',
                              'wilcoxon_r', 'cramer_v'),
           plot_type = ifelse(format == 'numeric',
                              'box', 'stack'),
           ax_label = ifelse(is.na(unit),
                             '% of cluster', unit))

  ## analysis data

  fgfr_clinic$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$clinic) %>%
    map(select, sample_id, any_of(fgfr_clinic$lexicon$variable))

  fgfr_clinic$data <-
    map2(fgfr_globals$genetic_groups[names(fgfr_clinic$data)],
         fgfr_clinic$data,
         inner_join, by = 'sample_id') %>%
    map(map_dfc, function(x) if(is.factor(x)) droplevels(x) else x)

  fgfr_clinic$data$bcan <- fgfr_clinic$data$bcan %>%
    mutate(tissue = ifelse(tissue == 'bladder',
                           'bladder', 'non-bladder'),
           tissue = factor(tissue, c('bladder', 'non-bladder')))

  ## data set-specific variable vectors

  fgfr_clinic$variables <- fgfr_clinic$data %>%
    map(names) %>%
    map(~filter(fgfr_clinic$lexicon, variable %in% .x)) %>%
    map(~.x$variable)

# N numbers -------

  insert_msg('Total numbers of samples')

  fgfr_clinic$n_numbers <- fgfr_clinic$data %>%
    map(count, FGFR3_mutation) %>%
    map(column_to_rownames, 'FGFR3_mutation') %>%
    map(t) %>%
    map(as.data.frame) %>%
    map(mutate, variable = 'Patients, N') %>%
    map(relocate, variable) %>%
    map(as_tibble)

# Descriptive stats -------

  insert_msg('Analysis stats')

  fgfr_clinic$stats <-
    future_map2(fgfr_clinic$data,
                fgfr_clinic$variables,
                ~explore(.x,
                         variables = .y,
                         split_factor = 'FGFR3_mutation',
                         what = 'table',
                         pub_styled = TRUE),
                .options = furrr_options(seed = TRUE)) %>%
    map(format_desc)

# Testing for differences between the mutation strata ------

  insert_msg('Testing for differences between the mutation strata')

  fgfr_clinic$test <-
    future_map2(fgfr_clinic$data,
                fgfr_clinic$variables,
                ~compare_variables(.x,
                                   variables = .y,
                                   split_factor = 'FGFR3_mutation',
                                   what = 'eff_size',
                                   types = exchange(.y,
                                                    fgfr_clinic$lexicon,
                                                    value = 'test_type'),
                                   exact = FALSE,
                                   ci = FALSE,
                                   pub_styled = TRUE,
                                   adj_method = 'BH'),
                .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '),
        significance = ifelse(stri_detect(eff_size, fixed = 'Inf'),
                              NA, significance),
        plot_cap = ifelse(stri_detect(eff_size, fixed = 'Inf'),
                          NA, plot_cap),
        eff_size = ifelse(stri_detect(eff_size, fixed = 'Inf'),
                          NA, eff_size))

# Plots --------

  insert_msg('Plots')

  for(i in names(fgfr_clinic$data)) {

    fgfr_clinic$plots[[i]] <-
      list(variable = fgfr_clinic$test[[i]]$variable,
           plot_title = fgfr_clinic$test[[i]]$variable %>%
             exchange(fgfr_clinic$lexicon) %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = fgfr_clinic$test[[i]]$plot_cap,
           y_lab = fgfr_clinic$test[[i]]$variable %>%
             exchange(fgfr_clinic$lexicon,
                      value = 'ax_label'),
           type = fgfr_clinic$test[[i]]$variable %>%
             exchange(fgfr_clinic$lexicon,
                      value = 'plot_type')) %>%
      future_pmap(plot_variable,
                  fgfr_clinic$data[[i]],
                  split_factor = 'FGFR3_mutation',
                  scale = 'percent',
                  cust_theme = globals$common_theme +
                    theme(axis.title.x = element_markdown(),
                          legend.title = element_markdown()),
                  x_lab = html_italic('FGFR3'),
                  x_n_labs = TRUE,
                  .options = furrr_options(seed = TRUE)) %>%
      set_names(fgfr_clinic$test[[i]]$variable)

    fgfr_clinic$plots[[i]]$age <- fgfr_clinic$plots[[i]]$age +
      scale_fill_manual(values = globals$mut_status_color,
                        name = html_italic('FGFR3'))

    fgfr_clinic$plots[[i]][names(fgfr_clinic$plots[[i]]) != 'age'] <-
      fgfr_clinic$plots[[i]][names(fgfr_clinic$plots[[i]]) != 'age'] %>%
      map(~.x +
            scale_fill_brewer(palette = 'Purples',
                              name = html_italic('FGFR3')))

  }

# Result table --------

  insert_msg('Result table')

  fgfr_clinic$result_tbl <-
    map2(fgfr_clinic$stats,
         fgfr_clinic$test %>%
           map(filter, !stri_detect(eff_size, fixed = 'Inf')),
         merge_stat_test) %>%
    map(format_res_tbl,
        dict = fgfr_clinic$lexicon) %>%
    map2(., fgfr_clinic$n_numbers,
         ~full_rbind(.y, .x)) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    relocate(cohort) %>%
    set_names(c('Cohort',
                'Variable',
                levels(fgfr_clinic$data[[1]]$FGFR3_mutation),
                'Significance',
                'Effect size'))

# END -------

  plan('sequential')

  rm(i)

  fgfr_clinic$data <- NULL

  fgfr_clinic <- compact(fgfr_clinic)

  insert_tail()
