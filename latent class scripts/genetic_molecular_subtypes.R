# distribution of the genetic clusters between the consensus molecular classes.
# the analysis is possible only for the TCGA cohort.

  insert_head()

# container --------

  lca_cons <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  lca_cons$data <-
    inner_join(lca_pred$assignment$tcga,
               subtypes$assignment$tcga[c('sample_id', 'consensusClass')],
               by = 'sample_id') %>%
    filter(complete.cases(.))

# Descriptive stats -------

  insert_msg('Descriptive stats')

  lca_cons$stats <- lca_cons$data %>%
    explore(variables = 'clust_id',
            split_factor = 'consensusClass',
            what = 'table',
            pub_styled = TRUE) %>%
    format_desc

# Testing for differences in cluster distribution between the consensus classes ------

  insert_msg('Testing')

  lca_cons$test <- lca_cons$data %>%
    compare_variables(variables = 'clust_id',
                      split_factor = 'consensusClass',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

# Stack plot ------

  insert_msg('Stack plot')

  lca_cons$plot <- lca_cons$data %>%
    plot_variable(variable = 'consensusClass',
                  split_factor = 'clust_id',
                  type = 'stack',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  plot_title = paste('Genetic clusters,',
                                     globals$cohort_labs["tcga"]),
                  plot_subtitle = lca_cons$test$plot_cap,
                  x_lab = 'MIBC consensus class',
                  y_lab = '% of subset',
                  x_n_labs = TRUE) +
    scale_fill_manual(values = globals$sub_colors)

# Result table --------

  insert_msg('Result table')

  lca_cons$result_tbl <-
    merge_stat_test(lca_cons$stats, lca_cons$test) %>%
    format_res_tbl(dict = NULL) %>%
    set_names(c('Variable',
                levels(lca_cons$data$consensusClass),
                'Significance', 'Effect size'))

# END ------

  lca_cons$data <- NULL

  lca_cons <- compact(lca_cons)

  insert_tail()
