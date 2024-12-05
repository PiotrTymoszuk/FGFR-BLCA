# Differences in distribution of the consensus molecular subtypes between
# the cohorts. Done for the TCGA BLCA, Imvigor, and BCAN cohorts

  insert_head()

# container -----

  expl_sub <- list()

# analysis data -------

  insert_msg('Analysis data')

  expl_sub$data <- subtypes$assignment %>%
    map(~.x[c('sample_id', 'consensusClass')]) %>%
    map(filter, !is.na(consensusClass), consensusClass != 'not assigned') %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = unname(globals$cohort_labs[cohort]),
           cohort = factor(cohort, unname(globals$cohort_labs)),
           cohort = droplevels(cohort),
           consensusClass = droplevels(consensusClass)) %>%
    filter(complete.cases(.))

# Descriptive stats ------

  insert_msg('Descriptive stats')

  expl_sub$stats <- expl_sub$data %>%
    explore(variables = 'consensusClass',
            split_factor = 'cohort',
            what = 'table',
            pub_styled = TRUE) %>%
    format_desc

# Testing for differences between the cohorts ------

  insert_msg('Test')

  expl_sub$test <- expl_sub$data %>%
    compare_variables(variables = 'consensusClass',
                      split_factor = 'cohort',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Stack plot ------

  insert_msg('Stack plot')

  expl_sub$plot <-
    plot_variable(expl_sub$data,
                  variable = 'consensusClass',
                  split_factor = 'cohort',
                  type = 'stack',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  plot_title = 'Consensus molecular subtypes',
                  plot_subtitle = expl_sub$test$plot_cap,
                  x_lab = 'Cohort',
                  y_lab = '% of cohort',
                  x_n_labs = TRUE) +
    scale_fill_manual(values = globals$sub_colors)

# Result table -------

  insert_msg('Result table')

  expl_sub$result_tbl <-
    merge_stat_test(expl_sub$stats,
                    expl_sub$test) %>%
    format_res_tbl %>%
    mutate(variable = 'Consensus molecular subtype')

# END -----

  expl_sub$data <- NULL

  expl_sub <- compact(expl_sub)

  insert_tail()
