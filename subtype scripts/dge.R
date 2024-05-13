# Differential expression of the FGFR and FGF genes between the molecular
# subtypes. The analysis is done for the TCGA and IMvigor.
# Differentially regulated genes are identified with the following rules:
#
# 1) Significant difference in expression between the classes by FDR-corrected
# one-way ANOVA (pFDR < 0.05).


  insert_head()

# container -------

  sub_dge <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data -------

  insert_msg('Analysis data')

  sub_dge$variables <- globals[c("receptors", "ligands")] %>%
    reduce(union)

  ## expression data

  sub_dge$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(select, sample_id, any_of(sub_dge$variables))

  sub_dge$variables <-
    sub_dge$variables[sub_dge$variables %in% names(sub_dge$data[[1]])]

  ## assignment to the molecular classes, setting the modeling
  ## intercept to LumP

  sub_dge$data <-
    map2(subtypes$assignment[names(sub_dge$data)],
         sub_dge$data,
         inner_join, by = 'sample_id') %>%
    map(mutate,
        consensusClass = fct_relevel(consensusClass, 'LumP'))

# Descriptive stats -------

  insert_msg('Descriptive stats')

  sub_dge$stats <- sub_dge$data %>%
    future_map(explore,
               variables = sub_dge$variables,
               split_factor = 'consensusClass',
               what = 'table',
               pub_styled = TRUE,
               .options = furrr_options(seed = TRUE)) %>%
    map(format_desc)

  plan('sequential')

# One-way ANOVA --------

  insert_msg('One-way ANOVA and linear modeling')

  sub_dge$test <- sub_dge$data %>%
    map(test_anova,
        split_fct = 'consensusClass',
        variables = sub_dge$variables,
        adj_method = 'BH')

  sub_dge$test <- transpose(sub_dge$test)

  ## formatting the ANOVA results to use them in plot captions
  ## and labels

  sub_dge$test$anova <- sub_dge$test$anova %>%
    map(re_adjust, method = 'none') %>%
    map(mutate,
        gene_symbol = html_italic(response),
        eff_size = paste('\u03B7\u00B2 =', signif(effect_size, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '),
        ax_lab = paste(gene_symbol, plot_cap, sep = '<br>'),
        ax_lab = ifelse(p_adjusted < 0.05,
                        html_bold(ax_lab), ax_lab))

  ## significant differences in ANOVA:

  sub_dge$anova_significant <- sub_dge$test$anova %>%
    map(filter,
        p_adjusted < 0.05,
        effect_size >= 0.14) %>%
    map(~.x$response)

  ## differences as compared with LumP class

  sub_dge$test$lm <- sub_dge$test$lm %>%
    map(filter, level != '(Intercept)') %>%
    map2(., sub_dge$anova_significant,
         ~mutate(.x,
                 regulation = ifelse(!response %in% .y | p_adjusted >= 0.05,
                                     'ns',
                                     ifelse(estimate >= log2(1.5),
                                            'upregulated',
                                            ifelse(estimate <= -log2(1.5),
                                                   'downregulated', 'ns'))),
                 regulation = factor(regulation,
                                     c('upregulated', 'downregulated', 'ns'))))

# Differentially regulated genes -------

  insert_msg('Differentially regulated genes')

  sub_dge$significant <- sub_dge$test$lm %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, regulation, level) %>%
    transpose %>%
    map(map, ~.x$response)

  ## significant effects shared by both cohorts

  sub_dge$common_significant <- sub_dge$significant %>%
    map(reduce, intersect)

# Box plots for single variables -------

  insert_msg('Box plots')

  for(i in names(sub_dge$data)) {

    sub_dge$plots[[i]] <-
      list(variable = sub_dge$test$anova[[i]]$response,
           plot_title = sub_dge$test$anova[[i]]$response %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = sub_dge$test$anova[[i]]$plot_cap) %>%
      pmap(plot_variable,
           sub_dge$data[[i]],
           split_factor = 'consensusClass',
           type = 'box',
           cust_theme = globals$common_theme,
           x_n_labs = TRUE,
           x_lab = 'MIBC consensus class',
           y_lab = expression('log'[2] * ' expression')) %>%
      map(~.x +
            theme(plot.title = element_markdown()) +
            scale_fill_manual(values = globals$sub_colors)) %>%
      set_names(sub_dge$test$anova[[i]]$response)

  }

# Heat map for the common differentially regulated genes ------

  insert_msg('Heat map for the selected genes')

  ## variables: common differentially regulated genes

  sub_dge$hm_variables <- reduce(sub_dge$common_significant, union)

  ## plotting order: classification for the four biggest classes
  ## in the TCGA cohort

  sub_dge$hm_classification <- sub_dge$data$tcga %>%
    filter(consensusClass %in% c('LumP', 'LumU', 'Stroma-rich', 'Ba/Sq')) %>%
    mutate(consensusClass = droplevels(consensusClass)) %>%
    classify(variables = sub_dge$hm_variables,
             split_fct = 'consensusClass')

  ## heat maps

  sub_dge$hm_plots <-
    list(data = sub_dge$data,
         plot_title = globals$cohort_labs[names(sub_dge$data)]) %>%
    pmap(heat_map,
         variables = sub_dge$hm_variables,
         split_fct = 'consensusClass',
         variable_classification = sub_dge$hm_classification$classification,
         cust_theme = globals$common_theme,
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales::squish,
         x_lab = 'sample')

  ## additional styling

  sub_dge$hm_plots <-
    list(x = sub_dge$hm_plots,
         y = sub_dge$test$anova) %>%
    pmap(function(x, y) x +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_markdown(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.x = element_text(hjust = 0, angle = 90),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank()) +
           scale_y_discrete(labels = set_names(y$ax_lab, y$response)))

# Result table -------

  insert_msg('Result table')

  sub_dge$n_number_tbl <- sub_dge$data %>%
    map(count, consensusClass) %>%
    map(column_to_rownames, 'consensusClass') %>%
    map(t) %>%
    map(as_tibble) %>%
    map(mutate, variable = 'Samples, N')

  sub_dge$result_tbl <-
    map2(sub_dge$stats,
         map(sub_dge$test$anova, mutate, variable = response),
         merge_stat_test) %>%
    map(format_res_tbl, dict = NULL, remove_complete = TRUE) %>%
    map2(sub_dge$n_number_tbl, .,
         full_rbind) %>%
    compress(names_to = 'cohort')

  sub_dge$result_tbl <- sub_dge$result_tbl %>%
    mutate(cohort = globals$cohort_labs[cohort]) %>%
    select(cohort, variable,
           all_of(levels(sub_dge$data[[1]]$consensusClass)),
           significance, eff_size) %>%
    set_names(c('Cohort', 'Variable',
                unname(globals$sub_labels),
                'Significance', 'Effect size'))

# END ------

  sub_dge$data <- NULL
  sub_dge$hm_classification <- NULL
  sub_dge$n_number_tbl <- NULL

  sub_dge <- compact(sub_dge)

  rm(i)

  plan('sequential')

  insert_tail()
