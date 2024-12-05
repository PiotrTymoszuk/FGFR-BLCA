# Frequency of genetic features in the molecular subtypes.
# The subtype information is done only for the TCGA cohort, where the expression
# and genetics information was provided.

  insert_head()

# container ------

  lca_sub <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis globals: data and variables -------

  insert_msg('Analysis variables')

  lca_sub$variables <- lca_globals$variables

  lca_sub$lexicon <- lca_globals$lexicon

  lca_sub$data <-
    inner_join(subtypes$assignment$tcga[c('sample_id', 'consensusClass')],
               lca_globals$data$tcga,
               by = 'sample_id') %>%
    filter(!is.na(consensusClass))

  ## background matrix with all alterations

  lca_sub$bcg_mtx <-
    tcga[c("mutation", "amplification", "deletion")] %>%
    map(column_to_rownames, 'sample_id') %>%
    map2(., c('_mutation', '_amplification', '_deletion'),
         ~set_names(.x,
                    paste0(names(.x), .y)))

  lca_sub$bcg_mtx <-  lca_sub$bcg_mtx %>%
    reduce(cbind)

  lca_sub$bcg_mtx <- lca_sub$bcg_mtx[rownames(lca_sub$data), ]

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

  ## enrichment testing with perich

  lca_sub$test <-
    perTest(lca_sub$data[lca_sub$variables],
            f = lca_sub$data$consensusClass,
            background = lca_sub$bcg_mtx,
            alternative = 'greater',
            n_iter = 10000,
            as_data_frame = TRUE,
            compress = TRUE,
            adj_method = 'BH',
            .parallel = TRUE)

  ## formatting the testing results

  lca_sub$test <- lca_sub$test %>%
    re_adjust('global_p_value') %>%
    ## global significance and plot caption
    mutate(eff_size = paste('V =', signif(cramer_v, 2)),
           global_significance = significance,
           plot_cap = paste(eff_size, global_significance, sep = ', ')) %>%
    ## significant enrichment
    mutate(regulation = ifelse(ES >= 1.2 & p_adjusted < 0.05,
                               'enriched', 'ns'),
           regulation = factor(regulation, c('enriched', 'ns'))) %>%
    ## gene symbol and alteration type
    mutate(gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           label = paste(html_italic(gene_symbol),
                         alteration,
                         sep = ', ')) %>%
    re_adjust %>%
    as_tibble

  ## significant effects: enrichment in at least one cluster

  lca_sub$significant <- lca_sub$test %>%
    filter(regulation == 'enriched') %>%
    .$variable %>%
    unique

  ## testing results to be presented in the plots

  lca_sub$stack_test <- lca_sub$test %>%
    filter(!duplicated(variable))

# Stack plots for single genetic features ------

  insert_msg('Stack plots for single genetic features')

  lca_sub$plots <-
    list(variable = lca_sub$stack_test$variable,
         plot_title = lca_sub$stack_test$gene_symbol %>%
           html_italic %>%
           paste(globals$cohort_labs["tcga"], sep = ', '),
         plot_subtitle = lca_sub$stack_test$plot_cap) %>%
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
    set_names(lca_sub$stack_test$variable)

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
    left_join(lca_sub$stats,
              lca_sub$stack_test[, c("variable",
                                     "alteration",
                                     "gene_symbol",
                                     "global_significance",
                                     "eff_size")],
              by = 'variable') %>%
    mutate(percent = signif(percent, 2)) %>%
    arrange(variable, consensusClass) %>%
    select(gene_symbol,
           alteration,
           consensusClass,
           percent,
           n,
           n_total,
           global_significance,
           eff_size) %>%
    set_names(c('Gene symbol', 'Alteration',
                'MIBC consensus class',
                'Percentage of samples with alterations',
                'Samples with alterations, N',
                'Total samples, N',
                'Significance',
                'Effects size'))

# Oncoplot of the significant alterations ------

  insert_msg('Oncoplot for snignificant alterations')

  lca_sub$onco_plot <-
    plot_bionco(lca_sub$data,
                variables = lca_sub$significant,
                split_fct = 'consensusClass',
                plot_title = globals$cohort_labs["tcga"],
                hide_x_axis_text = TRUE,
                cust_theme = globals$common_theme +
                  theme(axis.title.y = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        axis.line = element_blank(),
                        panel.grid.major.x = element_blank(),
                        strip.text.x = element_text(hjust = 0, angle = 90)),
                color_scale = c('gray85', 'orangered3'),
                one_plot = FALSE,
                x_lab = 'cancer sample')

  lca_sub$onco_plot <- lca_sub$onco_plot$main +
    theme(axis.text.y = element_markdown()) +
    scale_y_discrete(labels = set_names(lca_sub$test$label,
                                        lca_sub$test$variable))

# Caching --------

  insert_msg('Caching')

  lca_sub$data <- NULL
  lca_sub$stack_test <- NULL
  lca_sub$bcg_mtx <- NULL

  lca_sub <- compact(lca_sub)

  save(lca_sub, file = './cache/lca_sub.RData')

# END -------

  plan('sequential')

  insert_tail()
