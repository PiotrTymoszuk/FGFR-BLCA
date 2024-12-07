# Differential expression of the FGF and FGFR genes in the the genetic clusters.
# The analysis is done only for the TCGA cohort.
#
# Statistical significance is determined by one-way ANOVA with eta-square effect
# size metrics. FDR correction.

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

  lca_dge$variables <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    reduce(union)

  lca_dge$variables <-
    lca_dge$variables[lca_dge$variables %in% names(lca_dge$data)]

  lca_dge$data <- lca_dge$data %>%
    select(sample_id, clust_id, all_of(lca_dge$variables))

# Descriptive stats ------

  insert_msg('Descriptive stats')

  lca_dge$stats <- lca_dge$data %>%
    fast_num_stats(split_fct = 'clust_id',
                   variables = lca_dge$variables)
# One-way ANOVA ------

  insert_msg('One-way ANOVA')

  lca_dge$test <-
    f_one_anova(lca_dge$data[lca_dge$variables],
                f = lca_dge$data$clust_id,
                safely = TRUE,
                as_data_frame = TRUE,
                adj_method = 'BH') %>%
    re_adjust %>%
    mutate(eff_size = paste('\u03B7\u00B2 =', signif(etasq, 2)),
           plot_cap = paste(eff_size, significance, sep = ', '),
           ax_lab = paste(html_italic(variable),
                          plot_cap,
                          sep = '<br>'),
           ax_lab = ifelse(p_adjusted < 0.05,
                           html_bold(ax_lab), ax_lab)) %>%
    as_tibble

  ## significant effects

  lca_dge$significant <- lca_dge$test %>%
    filter(p_adjusted < 0.05, etasq >= 0.06) %>%
    .$variable

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
         y_lab = globals$cohort_axis_labs[[i]]) %>%
    map(~.x +
          theme(plot.title = element_markdown(),
                axis.title.y = element_markdown()) +
          scale_fill_manual(values = globals$genet_colors)) %>%
    set_names(lca_dge$test$variable)

# Result table ------

  insert_msg('Result table')

  lca_dge$result_tbl <-
    merge_stat_test(lca_dge$stats, lca_dge$test) %>%
    format_res_tbl(dict = NULL,
                   remove_complete = TRUE) %>%
    relocate(variable) %>%
    set_names(c('Variable',
                levels(lca_dge$data$clust_id),
                'Significance', 'Effect size'))

# END -------

  lca_dge$data <- NULL

  lca_dge <- compact(lca_dge)

  insert_tail()
