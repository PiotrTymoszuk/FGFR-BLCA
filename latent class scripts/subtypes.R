# Frequency of genetic features in the molecular subtypes.
# The subtype information is done only for the TCGA cohort, where the expression
# and genetics information was provided.

  insert_head()

# container ------

  lca_sub <- list()

# analysis globals: data and variables -------

  insert_msg('Analysis variables')

  lca_sub$variables <- lca_globals$variables

  lca_sub$lexicon <- lca_globals$lexicon

  lca_sub$data <-
    inner_join(subtypes$assignment$tcga[c('sample_id', 'consensusClass')],
               lca_globals$data$tcga,
               by = 'sample_id')

# N numbers -------

  insert_msg('Number od observations')

  ## numbers of samples in the molecular classes

  lca_sub$n_numbers <- lca_sub$data %>%
    count(consensusClass)

  ## a ready to use plot caption

  lca_sub$n_caption <-
    map2_chr(lca_sub$n_numbers[[1]],
             lca_sub$n_numbers[[2]],
             paste, sep = ': n = ') %>%
    paste(collapse = ', ')

# Descriptive stats ------

  insert_msg('Descriptive stats')

  lca_sub$stats <- lca_sub$data %>%
    select(-sample_id) %>%
    count_binary(split_fct = 'consensusClass')

# Testing for differences between the consensus classes -------

  insert_msg('Testing for differences between the consenus classes')

  lca_sub$test <-
    compare_variables(lca_sub$data,
                      variables = lca_sub$variables,
                      split_factor = 'consensusClass',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH',
                      .parallel = TRUE,
                      .paropts = furrr_options(seed = TRUE,
                                               globals = 'lca_sub')) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '),
           gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           plot_title = paste(html_italic(gene_symbol), alteration),
           hm_lab = paste(plot_title, plot_cap, sep = '<br>'))

  ## significant effects

  lca_sub$significant <- lca_sub$test %>%
    filter(p_adjusted < 0.05) %>%
    .$variable

# Stack plots for single genetic features ------

  insert_msg('Stack plots for single genetic features')

  lca_sub$plots <-
    list(variable = lca_sub$test$variable,
         plot_title = lca_sub$test$plot_title %>%
           paste(globals$cohort_labs["tcga"], sep = ', '),
         plot_subtitle = lca_sub$test$plot_cap) %>%
    future_pmap(plot_variable,
                lca_sub$data,
                split_factor = 'consensusClass',
                type = 'stack',
                scale = 'percent',
                cust_theme = globals$common_theme +
                  theme(plot.title = element_markdown()),
                x_n_labs = TRUE,
                y_lab = '% of subset',
                x_lab = 'MIBC consesus subset',
                .options = furrr_options(seed = TRUE)) %>%
    set_names(lca_sub$test$variable)

  ## styling

  lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_mut')] <-
    lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_mut')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$mut_status_color),
                            labels = names(globals$mut_status_color),
                            name = ''))

  lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_del')] <-
    lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_del')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$del_status_color),
                            labels = names(globals$del_status_color),
                            name = ''))

  lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_amp')] <-
    lca_sub$plots[stri_detect(names(lca_sub$plots), fixed = '_amp')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$amp_status_color),
                            labels = names(globals$amp_status_color),
                            name = ''))

# Result table --------

  insert_msg('Result table')

  lca_sub$result_tbl <-
    merge_stat_test(lca_sub$stats,
                    lca_sub$test) %>%
    mutate(gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           percent = signif(percent, 2)) %>%
    arrange(variable, consensusClass) %>%
    select(gene_symbol,
           alteration,
           consensusClass,
           percent,
           n,
           n_total,
           significance,
           eff_size) %>%
    set_names(c('Gene symbol', 'Alteration',
                'MIBC consensus class',
                'Percentage of samples with alterations',
                'Samples with alterations, N',
                'Total samples, N',
                'Significance',
                'Effects size'))

# Heat map representation of the significant alterations ------

  insert_msg('Heat map')

  ## classification of the variables for plotting by specificity
  ## for the four biggest consensus class

  lca_sub$hm_classification <- lca_sub$data %>%
    filter(consensusClass %in% c('LumP', 'LumU', 'Stroma-rich', 'Ba/Sq')) %>%
    mutate(consensusClass = droplevels(consensusClass)) %>%
    classify(variables = lca_sub$significant,
             split_fct = 'consensusClass') %>%
    .$classification %>%
    select(variable, consensusClass, auc) %>%
    set_names(c('variable', 'marker', 'auc'))

  lca_sub$hm_data <- lca_sub$data %>%
    select(sample_id, consensusClass, all_of(lca_sub$significant)) %>%
    pivot_longer(cols = all_of(lca_sub$significant),
                 names_to = 'variable',
                 values_to = 'alteration') %>%
    left_join(lca_sub$hm_classification, by = 'variable')

  ## heat maps

  lca_sub$hm_plot <- lca_sub$hm_data %>%
    ggplot(aes(x = reorder(sample_id, -alteration),
               y = variable,
               fill = factor(alteration))) +
    geom_tile() +
    facet_grid(marker ~ consensusClass,
               scales = 'free',
               space = 'free') +
    scale_y_discrete(labels = set_names(lca_sub$test$hm_lab,
                                        lca_sub$test$variable)) +
    scale_fill_manual(values = c('0' = 'gray70',
                                 '1' = 'orangered3'),
                      labels = c('0' = 'absent',
                                 '1' = 'present'),
                      name = '') +
    globals$common_theme +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_markdown(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(hjust = 0, angle = 90),
          strip.text.y = element_blank(),
          strip.background.y = element_blank()) +
    labs(title = globals$cohort_labs["tcga"],
         subtitle = lca_sub$n_caption,
         x = 'sample')

# END -------

  lca_sub$data <- NULL
  lca_sub$hm_data <- NULL
  lca_sub$hm_classification <- NULL

  lca_sub <- compact(lca_sub)

  insert_tail()
