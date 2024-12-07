# Differential expression of the FGFR and FGF genes between the molecular
# subtypes. The analysis is done for the TCGA and IMvigor.
# Differentially regulated genes are identified with the following rules:
#
# 1) Significant difference in expression between the classes by FDR-corrected
# one-way ANOVA (pFDR < 0.05).
#
# 2) Differential expression in a given class as compared with LumP is
# investigated by two-tailed T test.
#
# significantly differentially regulated genes are identified by pFDR(ANOVA) < 0.05


  insert_head()

# container -------

  sub_dge <- list()

# analysis data -------

  insert_msg('Analysis data')

  sub_dge$variables <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    reduce(union)

  ## expression data

  sub_dge$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(select, sample_id, any_of(sub_dge$variables))

  sub_dge$variables <- sub_dge$data %>%
    map(names) %>%
    map(intersect, sub_dge$variables)

  ## assignment to the molecular classes, setting the modeling
  ## intercept to LumP

  sub_dge$data <-
    map2(subtypes$assignment[names(sub_dge$data)],
         sub_dge$data,
         inner_join, by = 'sample_id') %>%
    map(filter, consensusClass != 'not assigned') %>%
    map(mutate, consensusClass = droplevels(consensusClass)) %>%
    map(mutate,
        consensusClass = fct_relevel(consensusClass, 'LumP'))

  ## BCAN: only one observation in the NE-like subset, removal of the subset
  ## in this cohort

  sub_dge$data$bcan <- sub_dge$data$bcan %>%
    filter(consensusClass != 'NE-like') %>%
    mutate(consensusClass = droplevels(consensusClass))

  sub_dge$data <- sub_dge$data %>%
    map(filter, !is.na(consensusClass))

  ## levels for T tests

  sub_dge$t_levels <- sub_dge$data %>%
    map(filter, consensusClass != 'LumP') %>%
    map(mutate, consensusClass = droplevels(consensusClass)) %>%
    map(~.x$consensusClass) %>%
    map(levels)

  sub_dge$t_levels <- sub_dge$t_levels %>%
    map(~set_names(.x, .x))

# Descriptive stats -------

  insert_msg('Descriptive stats')

  sub_dge$stats <-
    list(data = sub_dge$data,
         variables = sub_dge$variables) %>%
    pmap(fast_num_stats,
         split_fct = 'consensusClass')

# One-way ANOVA --------

  insert_msg('One-way ANOVA')

  sub_dge$anova <-
    map2(sub_dge$data,
         sub_dge$variables,
         ~f_one_anova(.x[.y],
                      f = .x[['consensusClass']],
                      safely = TRUE,
                      as_data_frame = TRUE,
                      adj_method = 'BH')) %>%
    map(re_adjust) %>%
    map(mutate,
        eff_size = paste('\u03B7\u00B2 =', signif(etasq, 2)),
        plot_cap = paste(eff_size, significance, sep = ', ')) %>%
    map(as_tibble)

  ## formatting and significant ANOVA effects

  sub_dge$anova <- sub_dge$anova %>%
    map(mutate,
        regulation = ifelse(p_adjusted < 0.05 & etasq >= 0.14,
                            'regulated', 'ns'),
        regulation = factor(regulation, c('regulated', 'ns'))) %>%
    map(mutate,
        ax_lab = paste(html_italic(variable), plot_cap, sep = '<br>'),
        ax_lab = ifelse(regulation == 'regulated',
                        html_bold(ax_lab), ax_lab),
        ax_lab_short = ifelse(regulation == 'regulated',
                              html_bold(html_italic(variable)),
                              html_italic(variable)))

  sub_dge$anova_significant <- sub_dge$anova %>%
    map(filter, regulation == 'regulated') %>%
    map(~.x$variable)

# Post-hoc T test -------

  insert_msg('Post-hoc T test')

  for(i in names(sub_dge$data)) {

    sub_dge$test[[i]] <- sub_dge$t_levels[[i]] %>%
      map(~filter(sub_dge$data[[i]], consensusClass %in% c('LumP', .x))) %>%
      map(mutate, consensusClass = droplevels(consensusClass))

    sub_dge$test[[i]] <- sub_dge$test[[i]] %>%
      map(~f_t_test(.x[sub_dge$variables[[i]]],
                    .x[['consensusClass']],
                    type = 'welch',
                    safely = TRUE,
                    as_data_frame = TRUE,
                    adj_method = 'BH')) %>%
      compress(names_to = 'consensusClass') %>%
      re_adjust %>%
      as_tibble

  }

  ## formatting the testing results

  sub_dge$test <-
    map2(sub_dge$test,
         sub_dge$anova_significant,
         ~mutate(.x,
                 anova_significant = ifelse(variable %in% .y,
                                            'yes', 'no'),
                 regulation = ifelse(anova_significant == 'no' | p_adjusted >= 0.05,
                                     'ns',
                                     ifelse(estimate >= log2(1.25),
                                            'upregulated',
                                            ifelse(estimate <= -log2(1.25),
                                                   'downregulated', 'ns'))),
                 regulation = factor(regulation,
                                     c('upregulated', 'downregulated', 'ns'))))

# Differentially regulated genes -------

  insert_msg('Differentially regulated genes')

  sub_dge$significant <- sub_dge$test %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, consensusClass) %>%
    map(map, blast, regulation) %>%
    transpose %>%
    map(compact) %>%
    map(transpose) %>%
    map(map, map, ~.x$variable)

  ## significant effects shared by at least two cohorts

  sub_dge$common_significant <- sub_dge$significant %>%
    map(map, shared_features, m = 2) %>%
    map(map, as.character)

# Box plots for single variables -------

  insert_msg('Box plots')

  for(i in names(sub_dge$data)) {

    sub_dge$plots[[i]] <-
      list(variable = sub_dge$anova[[i]]$variable,
           plot_title = sub_dge$anova[[i]]$variable %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = sub_dge$anova[[i]]$plot_cap) %>%
      pmap(plot_variable,
           sub_dge$data[[i]],
           split_factor = 'consensusClass',
           type = 'box',
           cust_theme = globals$common_theme,
           x_n_labs = TRUE,
           x_lab = 'MIBC consensus class',
           y_lab = globals$cohort_axis_labs[[i]]) %>%
      map(~.x +
            theme(plot.title = element_markdown(),
                  axis.title.y = element_markdown()) +
            scale_fill_manual(values = globals$sub_colors)) %>%
      set_names(sub_dge$anova[[i]]$variable)

  }

# Heat map for the common differentially regulated genes ------

  insert_msg('Heat map for the selected genes')

  ## variables: common differentially regulated genes

  sub_dge$hm_variables <- sub_dge$common_significant %>%
    unlist(use.names = FALSE) %>%
    unique

  ## plotting order: classification for the four biggest classes
  ## in the TCGA cohort

  sub_dge$hm_classification <- sub_dge$data$tcga %>%
    filter(consensusClass %in% c('LumP', 'LumU', 'Stroma-rich', 'Ba/Sq')) %>%
    mutate(consensusClass = droplevels(consensusClass)) %>%
    classify(variables = sub_dge$hm_variables,
             split_fct = 'consensusClass')

  ## heat maps

  sub_dge$hm_plots <-
    list(data = sub_dge$data %>%
           map(mutate,
               consensusClass = fct_recode(consensusClass,
                                           `Stroma\nrich` = 'Stroma-rich')),
         plot_title = globals$cohort_labs[names(sub_dge$data)]) %>%
    pmap(heat_map,
         variables = sub_dge$hm_variables,
         split_fct = 'consensusClass',
         variable_classification = sub_dge$hm_classification$classification,
         cust_theme = globals$common_theme,
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales::squish,
         x_lab = 'cancer sample')

  ## additional styling

  sub_dge$hm_plots <-
    list(x = sub_dge$hm_plots,
         y = sub_dge$anova) %>%
    pmap(function(x, y) x +
           labs(subtitle = stri_replace(x$labels$subtitle,
                                        fixed = '\n',
                                        replacement = '-')) +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_markdown(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.x = element_text(hjust = 0, angle = 90),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank()) +
           scale_y_discrete(labels = set_names(y$ax_lab_short, y$variable)))

# Result table -------

  insert_msg('Result table')

  sub_dge$result_tbl <-
    map2(sub_dge$stats,
         sub_dge$anova,
         merge_stat_test) %>%
    map(format_res_tbl, dict = NULL) %>%
    map2(., names(.),
         ~mutate(.x, cohort = globals$cohort_labs[.y])) %>%
    reduce(full_rbind)

  sub_dge$result_tbl <- sub_dge$result_tbl %>%
    select(cohort, variable,
           all_of(levels(sub_dge$data[[1]]$consensusClass)),
           significance, eff_size) %>%
    set_names(c('Cohort', 'Variable',
                globals$sub_labels[levels(sub_dge$data[[1]]$consensusClass)],
                'Significance', 'Effect size'))

# END ------

  sub_dge$data <- NULL
  sub_dge$hm_classification <- NULL
  sub_dge$n_number_tbl <- NULL

  sub_dge <- compact(sub_dge)

  rm(i)

  plan('sequential')

  insert_tail()
