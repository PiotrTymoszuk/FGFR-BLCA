# Best overall response to ICB in the IMvigor cohort stratified by
# prsence/absence of FGFR3 mutations

  insert_head()

# container -------

  fgfr_icb <- list()

# analysis data -------

  insert_msg('Analysis data')

  fgfr_icb$data <-
    inner_join(fgfr_globals$genetic_groups$imvigor,
               imvigor$clinic[c('sample_id', 'bi_response')],
               by = 'sample_id') %>%
    filter(complete.cases(.))

# descriptive stats ------

  insert_msg('Descriptive stats')

  fgfr_icb$stats <- fgfr_icb$data %>%
    explore(variable = 'bi_response',
            split_factor = 'FGFR3_mutation',
            what = 'table',
            pub_styled = TRUE) %>%
    format_desc()

# Test ------

  insert_msg('Test')

  fgfr_icb$test <- fgfr_icb$data %>%
    compare_variables(variables = 'bi_response',
                      split_factor = 'FGFR3_mutation',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Stack plot --------

  insert_msg('Stack plot')

  fgfr_icb$plot <- fgfr_icb$data %>%
    plot_variable(variable = 'bi_response',
                  split_factor = 'FGFR3_mutation',
                  type = 'stack',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  plot_title = paste('Best overall response to ICB,',
                                     globals$cohort_labs["imvigor"]),
                  plot_subtitle = fgfr_icb$test$plot_cap,
                  x_lab = html_italic('FGFR3'),
                  y_lab = '% of observations',
                  x_n_labs = TRUE) +
    scale_fill_brewer(palette = 'Purples') +
    theme(axis.title.x = element_markdown())

# Result table -------

  insert_msg('Result table')

  fgfr_icb$result_tbl <-
    merge_stat_test(fgfr_icb$stats, fgfr_icb$test) %>%
    format_res_tbl(dict = NULL) %>%
    mutate(variable = 'Best overall response to ICB') %>%
    set_names(c('Variable', levels(fgfr_icb$data$FGFR3_mutation),
                'Significance', 'Effect size'))

# END -------

  fgfr_icb$data <- NULL

  fgfr_icb <- compact(fgfr_icb)

  insert_tail()
