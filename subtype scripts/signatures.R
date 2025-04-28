# Comparison of ssGSEA scores of signatures of FGFR-related genes between
# consensus molecular subtypes of MIBC.
#
# Statistical significance of overall differences between the subsets is
# determined by one-way ANOVA with eta-square effect size statistic.
# Differences between single subsets (LumU, Stroma-rich, Ba/Sq) and the LumP
# class are assessed by Welch's T test with Cohen's d effect size statistic.
# Differentially regulated signatures are defined by pFDR(ANOVA) < 0.05,
# eta-square >= 0.14, pFDR(T test) < 0.05, and Cohen's d >= 0.5.

  insert_tail()

# container ------

  sub_sig <- list()

# analysis globals ----

  insert_msg('Analysis globals')

  ## variables and their lexicon

  sub_sig$lexicon <- sig$lexicon %>%
    mutate(gene_string = map_chr(genes, paste, collapse = ', '))

  sub_sig$variables <- sub_sig$lexicon$variable

  ## ssGSEA scores appended with the consensus class assignment

  sub_sig$data <- sig$scores

  sub_sig$data <-
    map2(subtypes$assignment[names(sub_sig$data)],
         sub_sig$data,
         inner_join, by = 'sample_id') %>%
    map(filter,
        !consensusClass %in% c('NE-like', 'LumNS', 'not assigned')) %>%
    map(mutate, consensusClass = droplevels(consensusClass)) %>%
    map(mutate, consensusClass = fct_relevel(consensusClass, 'LumP')) %>%
    map(filter, !is.na(consensusClass))

  ## levels for T tests

  sub_sig$t_levels <- sub_sig$data %>%
    map(filter, consensusClass != 'LumP') %>%
    map(mutate, consensusClass = droplevels(consensusClass)) %>%
    map(~.x$consensusClass) %>%
    map(levels)

  sub_sig$t_levels <- sub_sig$t_levels %>%
    map(~set_names(.x, .x))

# descriptive stats --------

  insert_msg('Descriptive stats')

  sub_sig$stats <- sub_sig$data %>%
    map(fast_num_stats,
        variables = sub_sig$variables,
        split_fct = 'consensusClass')

# One-way ANOVA --------

  insert_msg('One-way ANOVA')

  sub_sig$anova <- sub_sig$data %>%
    map(~f_one_anova(.x[sub_sig$variables],
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

  sub_sig$anova <- sub_sig$anova %>%
    map(mutate,
        regulation = ifelse(p_adjusted < 0.05 & etasq >= 0.14,
                            'regulated', 'ns'),
        regulation = factor(regulation, c('regulated', 'ns'))) %>%
    map(mutate,
        ax_lab = exchange(variable, sub_sig$lexicon),
        ax_lab = paste(ax_lab, plot_cap, sep = '<br>'),
        ax_lab = ifelse(regulation == 'regulated',
                        html_bold(ax_lab), ax_lab),
        ax_lab_short = exchange(variable, sub_sig$lexicon),
        ax_lab_short = ifelse(regulation == 'regulated',
                              html_bold(ax_lab_short), ax_lab_short),
        star_lab = ifelse(regulation == 'regulated',
                          html_bold(plot_cap),
                          paste0('<span style = "color: #585858">',
                                 plot_cap, '</span>')))

  sub_sig$anova_significant <- sub_sig$anova %>%
    map(filter, regulation == 'regulated') %>%
    map(~.x$variable)

# Post-hoc T test -------

  insert_msg('Post-hoc T test')

  for(i in names(sub_sig$data)) {

    sub_sig$test[[i]] <- sub_sig$t_levels[[i]] %>%
      map(~filter(sub_sig$data[[i]], consensusClass %in% c('LumP', .x))) %>%
      map(mutate, consensusClass = droplevels(consensusClass))

    sub_sig$test[[i]] <- sub_sig$test[[i]] %>%
      map(~f_t_test(.x[sub_sig$variables],
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

  sub_sig$test <-
    map2(sub_sig$test,
         sub_sig$anova_significant,
         ~mutate(.x,
                 anova_significant = ifelse(variable %in% .y,
                                            'yes', 'no'),
                 regulation = ifelse(anova_significant == 'no' | p_adjusted >= 0.05,
                                     'ns',
                                     ifelse(cohen_d >= 0.5,
                                            'upregulated',
                                            ifelse(cohen_d <= -0.5,
                                                   'downregulated', 'ns'))),
                 regulation = factor(regulation,
                                     c('upregulated', 'downregulated', 'ns'))))

# Differentially regulated signatures -------

  insert_msg('Differentially regulated signatures')

  sub_sig$significant <- sub_sig$test %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, consensusClass) %>%
    map(map, blast, regulation) %>%
    transpose %>%
    map(compact) %>%
    map(transpose) %>%
    map(map, map, ~.x$variable)

  ## significant effects shared by at least two cohorts

  sub_sig$common_significant <- sub_sig$significant %>%
    map(map, shared_features, m = 2) %>%
    map(map, as.character)

  ## top common significant effects:
  ## the largest effect sizes in ANOVA

  sub_sig$top_common_significant <- sub_sig$anova %>%
    map(filter, variable %in% unique(unlist(sub_sig$common_significant))) %>%
    map(slice_max, etasq, n = 15) %>%
    map(~.x$variable) %>%
    reduce(intersect)

# Heat map for the common significant signatures ------

  insert_msg('Heat map for the common significant signatures')

  sub_sig$hm_variables <- sub_sig$common_significant %>%
    unlist(use.names = FALSE) %>%
    unique

  ## heat maps

  sub_sig$hm_plots <-
    list(data = sub_sig$data %>%
           map(mutate,
               consensusClass = fct_recode(consensusClass,
                                           `Stroma\nrich` = 'Stroma-rich')),
         plot_title = globals$cohort_labs[names(sub_sig$data)]) %>%
    pmap(heat_map,
         variables = sub_sig$hm_variables,
         split_fct = 'consensusClass',
         variable_classification = NULL,
         normalize = FALSE,
         cust_theme = globals$common_theme,
         midpoint = 0,
         limits = c(-0.8, 0.8),
         oob = scales::squish,
         x_lab = 'cancer sample',
         fill_lab = 'ssGSEA')

  ## additional styling

  sub_sig$hm_plots <-
    list(x = sub_sig$hm_plots,
         y = sub_sig$anova) %>%
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
           scale_y_discrete(labels = set_names(y$star_lab,
                                               y$variable)))

# Result table -------

  insert_msg('Result table')

  sub_sig$result_tbl <-
    map2(sub_sig$stats,
         sub_sig$anova,
         merge_stat_test) %>%
    map(format_res_tbl, dict = NULL) %>%
    map2(., names(.),
         ~mutate(.x, cohort = globals$cohort_labs[.y])) %>%
    reduce(full_rbind)

  ## appending the samples with gene members of the signatures
  ## and user friendly signature names

  sub_sig$result_tbl <- sub_sig$result_tbl %>%
    mutate(genes = exchange(variable, sub_sig$lexicon, value = 'gene_string'),
           variable = exchange(variable, sub_sig$lexicon),
           variable = ifelse(is.na(variable), 'Samples, N', variable))

  ## final styling

  sub_sig$result_tbl <- sub_sig$result_tbl %>%
    select(cohort, variable, genes,
           all_of(levels(sub_sig$data[[1]]$consensusClass)),
           significance, eff_size) %>%
    set_names(c('Cohort', 'Variable', 'Member genes',
                globals$sub_labels[levels(sub_sig$data[[1]]$consensusClass)] %>%
                  stri_capitalize_first,
                'Significance', 'Effect size'))

# END ------

  sub_sig$data <- NULL
  sub_sig$hm_classification <- NULL
  sub_sig$n_number_tbl <- NULL

  sub_sig <- compact(sub_sig)

  rm(i)

  plan('sequential')

  insert_tail()
