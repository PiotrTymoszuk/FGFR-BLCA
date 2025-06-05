# Demographic, clinical and pathological characteristic of the cohorts.

  insert_tail()

# container ------

  expl_cohorts <- list()

# analysis variables and data -------

  insert_msg('Analysis variables and data')

  ## for BCAN: collapsing the non-bladder and metastasis organs

  expl_cohorts$data <- globals$cohort_expr %>%
    eval %>%
    c(list(hpa = hpa)) %>%
    map(~.x$clinic)

  expl_cohorts$data$bcan <- expl_cohorts$data$bcan %>%
    mutate(tissue = ifelse(tissue == 'bladder',
                           'bladder', 'non-bladder'),
           tissue = factor(tissue, c('bladder', 'non-bladder')))

  ## variable lexicon: simplified variable labels for gender, pM and pN stage

  expl_cohorts$lexicon <- globals$clinic_lexicon %>%
    filter(!variable %in% c("tissue", "invasiveness")) %>%
    mutate(label = ifelse(variable == "sex",
                          "Male gender", label),
           label = ifelse(variable == "pm_stage",
                          "Distant metastasis", label),
           label = ifelse(variable == "pn_stage",
                          "Lymph node metastasis", label)) %>%
    mutate(unit = ifelse(variable == 'tmb_per_mb',
                         'mutations/MB or fraction of genome',
                         unit),
           table_label = ifelse(!is.na(unit),
                                paste(label, unit, sep = ', '),
                                label),
           axis_lab = ifelse(is.na(unit), '% of cohort', unit),
           test_type = ifelse(format == 'numeric',
                              'kruskal_etasq', 'cramer_v'),
           plot_type = ifelse(format == 'numeric',
                              'box', 'stack'))

  ## variable vectors

  expl_cohorts$variables <- expl_cohorts$data %>%
    map(names) %>%
    map(~.x[.x %in% expl_cohorts$lexicon$variable])

  expl_cohorts$level_numbers <- expl_cohorts$variables %>%
    unlist %>%
    unname %>%
    table %>%
    as.data.frame %>%
    set_names(c('variable', 'n_levels'))

  expl_cohorts$lexicon <-
    left_join(expl_cohorts$lexicon,
              expl_cohorts$level_numbers,
              by = 'variable') %>%
    mutate(test_type = ifelse(n_levels == 1, 'none', test_type),
           test_type = ifelse(test_type == 'kruskal_etasq' & n_levels == 2,
                              'wilcoxon_r', test_type),
           test_type = ifelse(variable %in% c('tmb_per_mb', 'invasiveness'),
                              'none', test_type))

  ## conversion of characters to factors

  expl_cohorts$data <- expl_cohorts$data %>%
    map(map_dfc, function(x) if(is.character(x)) factor(x) else x)

  ## a common analysis data frame:
  ## re-coding the mortality index as factor,
  ## simplifying the gender, pN and pM variables

  expl_cohorts$data <- expl_cohorts$data %>%
    map2_dfr(., names(.),
             ~mutate(.x, cohort = .y)) %>%
    mutate(cohort = unname(globals$cohort_labs[cohort]),
           cohort = factor(cohort, unname(globals$cohort_labs)),
           death = as.character(death),
           death = car::recode(death, "'0' = 'no'; '1' = 'yes'"),
           death = factor(death, c('no', 'yes')),
           sex = ifelse(sex == "male", "yes", "no"),
           sex = factor(sex, c("no", "yes")),
           pn_stage = ifelse(pn_stage == "N0", "no", "yes"),
           pn_stage = factor(pn_stage, c("no", "yes")),
           pm_stage = ifelse(pm_stage == "M0", "no", "yes"),
           pm_stage = factor(pm_stage, c("no", "yes")))

  expl_cohorts$data <- expl_cohorts$data %>%
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x) %>%
    map_dfc(unname)

# Total numbers of cancer samples ------

  insert_msg('Total numbers of cancer samples')

  expl_cohorts$n_numbers <- expl_cohorts$data %>%
    count(cohort) %>%
    column_to_rownames('cohort') %>%
    t %>%
    as_tibble %>%
    mutate(variable = 'Samples, N') %>%
    relocate(variable)

# Descriptive stats -------

  insert_msg('Descriptive stats')

  expl_cohorts$stats <- expl_cohorts$data %>%
    explore(variables = expl_cohorts$lexicon$variable,
            split_factor = 'cohort',
            what = 'table',
            pub_styled = TRUE) %>%
    map(mutate,
        statistic = ifelse(stri_detect(statistic, regex = 'NA|NaN'),
                           NA, statistic)) %>%
    format_desc

# Testing for differences between the cohorts --------

  insert_msg('Testing')

  expl_cohorts$test_lexicon <- expl_cohorts$lexicon %>%
    filter(test_type != 'none')

  expl_cohorts$test <-
    list(x = expl_cohorts$test_lexicon$variable,
         y = expl_cohorts$test_lexicon$test_type) %>%
    pmap_dfr(function(x, y) expl_cohorts$data[, c('cohort', x)] %>%
               filter(complete.cases(.)) %>%
               compare_variables(variables = x,
                                 split_factor = 'cohort',
                                 what = 'eff_size',
                                 types = y,
                                 ci = FALSE,
                                 exact = FALSE,
                                 pub_styled = TRUE)) %>%
    re_adjust %>%
    p_formatter(text = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Result table -------

  insert_msg('Result table')

  expl_cohorts$result_tbl <-
    merge_stat_test(expl_cohorts$stats,
                    expl_cohorts$test %>%
                      filter(!stri_detect(eff_size, fixed = 'Inf'))) %>%
    format_res_tbl(dict = expl_cohorts$lexicon,
                   value = 'table_label')

  ## appending with the N numbers of samples

  expl_cohorts$result_tbl <-
    full_rbind(expl_cohorts$n_numbers,
               expl_cohorts$result_tbl) %>%
    set_names(c('Variable', levels(expl_cohorts$data$cohort),
                'FDR p value', 'Effect size'))

# Plots ---------

  insert_msg('Plots')

  expl_cohorts$plots <-
    list(variable = expl_cohorts$test$variable,
         plot_title = expl_cohorts$test$variable %>%
           exchange(expl_cohorts$lexicon),
         plot_subtitle = expl_cohorts$test$plot_cap,
         type = expl_cohorts$test$variable %>%
           exchange(expl_cohorts$lexicon,
                    value = 'plot_type'),
         y_lab = expl_cohorts$test$variable %>%
           exchange(expl_cohorts$lexicon,
                    value = 'axis_lab')) %>%
    pmap(plot_variable,
         expl_cohorts$data,
         split_factor = 'cohort',
         cust_theme = globals$common_theme,
         x_n_labs = TRUE,
         point_hjitter = 0,
         scale = 'percent',
         x_lab = 'Cohort') %>%
  set_names(expl_cohorts$test$variable)

  ## styling

  expl_cohorts$plots$age <- expl_cohorts$plots$age +
    scale_fill_manual(values = unname(globals$cohort_colors[c("genie",
                                                              "msk",
                                                              "tcga",
                                                              "bcan")]))

  expl_cohorts$plots[names(expl_cohorts$plots) != 'age'] <-
    expl_cohorts$plots[names(expl_cohorts$plots) != 'age'] %>%
    map(~.x +
          scale_fill_brewer('Blues'))

# END ------

  insert_tail()
