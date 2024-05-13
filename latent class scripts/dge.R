# Differential expression of the FGF and FGFR genes in the the genetic clusters.
# The analysis is done only for the TCGA cohort

  insert_head()

# container -------

  lca_dge <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  lca_dge$data <-
    inner_join(lca_pred$assignment$tcga,
               tcga$expression,
               by = 'sample_id')

  ## analysis variables

  lca_dge$variables <- globals[c("receptors", "ligands")] %>%
    reduce(union)

  lca_dge$variables <-
    lca_dge$variables[lca_dge$variables %in% names(lca_dge$data)]

  lca_dge$data <- lca_dge$data %>%
    select(sample_id, clust_id, all_of(lca_dge$variables))

# Descriptive stats ------

  insert_msg('Descriptive stats')

  lca_dge$stats <- lca_dge$data %>%
    explore(variables = lca_dge$variables,
            split_factor = 'clust_id',
            what = 'table',
            pub_styled = TRUE) %>%
    format_desc

# One-way ANOVA ------

  insert_msg('One-way ANOVA')

  lca_dge$test <- lca_dge$data %>%
    test_anova(split_fct = 'clust_id',
               variables = lca_dge$variables,
               adj_method = 'BH') %>%
    .$anova %>%
    re_adjust(method = 'none') %>%
    mutate(eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
           plot_cap = paste(eff_size, significance, sep = ', '),
           variable = response,
           ax_lab = paste(html_italic(response),
                          plot_cap,
                          sep = '<br>'),
           ax_lab = ifelse(p_adjusted < 0.05,
                           html_bold(ax_lab), ax_lab))

  ## significant effects

  lca_dge$significant <- lca_dge$test %>%
    filter(p_adjusted < 0.05) %>%
    .$response

# Box plots for single variables --------

  insert_msg('Box plots for single variables')

  lca_dge$plots <-
    list(variable = lca_dge$test$variable,
         plot_title = lca_dge$test$variable %>%
           html_italic %>%
           paste(globals$cohort_labs["tcga"], sep = ', '),
         plot_subtitle = lca_dge$test$plot_cap) %>%
    pmap(plot_variable,
         lca_dge$data,
         split_factor = 'clust_id',
         type = 'box',
         cust_theme = globals$common_theme,
         x_n_labs = TRUE,
         x_lab = 'Genetic cluster',
         y_lab = expression('log'[2] * ' expression')) %>%
    map(~.x +
          theme(plot.title = element_markdown()) +
          scale_fill_manual(values = globals$genet_colors)) %>%
    set_names(lca_dge$test$variable)

# Result table ------

  insert_msg('Result tables')

  ## n numbers

  lca_dge$n_numbers <- lca_dge$data %>%
    count(clust_id) %>%
    column_to_rownames('clust_id') %>%
    t %>%
    as_tibble %>%
    mutate(variable = 'Samples, N')

  lca_dge$result_tbl <-
    merge_stat_test(lca_dge$stats, lca_dge$test) %>%
    format_res_tbl(dict = NULL,
                   remove_complete = TRUE) %>%
    full_rbind(lca_dge$n_numbers, .) %>%
    relocate(variable) %>%
    set_names(c('Variable',
                levels(lca_dge$data$clust_id),
                'Significance', 'Effect size'))

# END -------

  lca_dge$data <- NULL

  lca_dge <- compact(lca_dge)

  insert_tail()
