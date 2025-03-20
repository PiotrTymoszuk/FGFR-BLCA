# Comparison of metrics of predicted resistance (IC50 and AUC) to anti-cancer
# drugs between consensus classes of bladder cancer in the TCGA, IMvigor, and
# BCAN cohorts.
#
# Statistical significance is determined by Kruskal-Wallis test with eta-square
# effect size statistic (multiple predictions in the BCAN cohort are clearly
# non normally distributed).
# Post-hoc test: Mann-Whitney test against the LumP subset with biserial r
# effect size metric.
#
# Significant differences: pFDR(Kruskal-Wallis) < 0.05 and eta-square >= 0.14
# and significant differences between the consensus classes with
# pFDR(Mann-Whitney) < 0.05 and r >= 0.3.

  insert_head()

# container ------

  sub_drugs <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## drug variables and their lexicon

  sub_drugs$lexicon <- drugs$lexicon %>%
    mutate(axis_title = paste('predicted resistance,', unit))

  sub_drugs$variables <- sub_drugs$lexicon$variable

  ## drugs that target FGFR

  sub_drugs$fgfr_drugs <- drugs$fgfr_drugs

  ## the predicted resistance metrics, appended with
  ## the consensus class assignment

  sub_drugs$data <- drugs$predictions

  sub_drugs$data <-
    map2(subtypes$assignment[names(sub_drugs$data)],
         sub_drugs$data,
         inner_join, by = 'sample_id')

  ## removal of the NE-like and LumNS subsets: not enough observations
  ## to evaluate

  sub_drugs$data <- sub_drugs$data %>%
    map(filter,
        !consensusClass %in% c('NE-like', 'LumNS', 'not assigned')) %>%
    map(mutate,
        consensusClass = droplevels(consensusClass))

  sub_drugs$data <- sub_drugs$data %>%
    map(filter, !is.na(consensusClass))

  ## levels for Mann-Whitney post-hoc tests

  sub_drugs$posthoc_levels <- sub_drugs$data %>%
    map(filter, consensusClass != 'LumP') %>%
    map(mutate, consensusClass = droplevels(consensusClass)) %>%
    map(~.x$consensusClass) %>%
    map(levels)

  sub_drugs$posthoc_levels <- sub_drugs$posthoc_levels %>%
    map(~set_names(.x, .x))

  ## title suffixes

  sub_drugs$title_suffixes <-
    c(all = 'all drugs',
      common_significant = 'significant in at least two cohorts',
      fgfr_drugs = 'FGFR-targeting drugs')

# Descriptive stats -------

  insert_msg('Descriptive stats')

  sub_drugs$stats <- sub_drugs$data %>%
    map(fast_num_stats,
        variables = sub_drugs$variables,
        split_fct = 'consensusClass')

# Kruskal-Wallis tests ------

  insert_msg('Kruskal-Wallis tests')

  sub_drugs$kruskal <- sub_drugs$data %>%
    map(~f_kruskal_test(.x[sub_drugs$variables],
                        .x$consensusClass,
                        as_data_frame = TRUE,
                        adj_method = 'BH',
                        safely = TRUE)) %>%
    map(re_adjust) %>%
    map(mutate,
        eff_size = paste('\u03B7\u00B2 =', signif(etasq, 2)),
        plot_cap = paste(eff_size, significance, sep = ', ')) %>%
    map(as_tibble)

  ## formatting and significant ANOVA effects

  sub_drugs$kruskal <- sub_drugs$kruskal %>%
    map(mutate,
        regulation = ifelse(p_adjusted < 0.05 & etasq >= 0.14,
                            'regulated', 'ns'),
        regulation = factor(regulation, c('regulated', 'ns'))) %>%
    map(mutate,
        title_label = exchange(variable,
                               sub_drugs$lexicon,
                               value = 'plot_title'),
        hm_label = exchange(variable,
                            sub_drugs$lexicon,
                            value = 'hm_label'),
        hm_label_long = paste(hm_label, plot_cap, sep = '<br>'),
        hm_label = ifelse(regulation == 'regulated',
                          html_bold(hm_label),
                          hm_label),
        hm_label_long = ifelse(regulation == 'regulated',
                               html_bold(hm_label_long),
                               hm_label_long))

  sub_drugs$kruskal_significant <- sub_drugs$kruskal %>%
    map(filter, regulation == 'regulated') %>%
    map(~.x$variable)

# Post-hoc Mann-Whitney tests --------

  insert_msg('Post-hoc Mann-whitney tests')

  for(i in names(sub_drugs$data)) {

    sub_drugs$test[[i]] <- sub_drugs$posthoc_levels[[i]] %>%
      map(~filter(sub_drugs$data[[i]], consensusClass %in% c('LumP', .x))) %>%
      map(mutate, consensusClass = droplevels(consensusClass))

    sub_drugs$test[[i]] <- sub_drugs$test[[i]] %>%
      map(~f_wilcox_test(.x[sub_drugs$variables],
                         .x$consensusClass,
                         exact = FALSE,
                         safely = TRUE,
                         as_data_frame = TRUE,
                         adj_method = 'BH')) %>%
      compress(names_to = 'consensusClass') %>%
      re_adjust %>%
      as_tibble

  }

  ## formatting the testing results: regulation

  sub_drugs$test <-
    map2(sub_drugs$test,
         sub_drugs$kruskal_significant,
         ~mutate(.x,
                 kruskal_significant = ifelse(variable %in% .y,
                                            'yes', 'no'),
                 regulation = ifelse(kruskal_significant == 'no' | p_adjusted >= 0.05,
                                     'ns',
                                     ifelse(biserial_r >= 0.3,
                                            'resistant',
                                            ifelse(biserial_r <= -0.3,
                                                   'sensitive', 'ns'))),
                 regulation = factor(regulation,
                                     c('resistant', 'sensitive', 'ns'))))

  ## formatting the testing results: levels of the consensusClasses

  sub_drugs$test <-
    map2(sub_drugs$test,
         sub_drugs$posthoc_levels,
         ~mutate(.x,
                 consensusClass = factor(consensusClass, levels = .y)))

  ## formatting the testing results: plotting labels

  sub_drugs$test <- sub_drugs$test %>%
    map(mutate,
        eff_size = paste('r =', signif(biserial_r, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '),
        title_label = exchange(variable,
                               sub_drugs$lexicon,
                               value = 'plot_title'),
        hm_label = exchange(variable,
                            sub_drugs$lexicon,
                            value = 'hm_label'),
        hm_label_long = paste(hm_label, plot_cap, sep = '<br>'),
        hm_label = ifelse(regulation %in% c('resistant', 'sensitive'),
                          html_bold(hm_label),
                          hm_label),
        hm_label_long = ifelse(regulation %in% c('resistant', 'sensitive'),
                               html_bold(hm_label_long),
                               hm_label_long))

# Significant differences in predicted drug resistance -------

  insert_msg('Significant differences in predicted drug resistance')

  sub_drugs$significant <- sub_drugs$test %>%
    map(filter, regulation %in% c('resistant', 'sensitive')) %>%
    map(blast, consensusClass) %>%
    map(map, blast, regulation) %>%
    transpose %>%
    map(compact) %>%
    map(transpose) %>%
    map(map, map, ~.x$variable)

  ## significant effects shared by at least two cohorts

  sub_drugs$common_significant <- sub_drugs$significant %>%
    map(map, shared_features, m = 2) %>%
    map(map, as.character)

# Volcano plots for the common significant effects ------

  insert_msg('Volcano plots')

  for(i in names(sub_drugs$test)) {

    sub_drugs$volcano_plots[[i]] <-
      list(data = map2(sub_drugs$test[[i]] %>%
                         blast(consensusClass),
                       sub_drugs$common_significant %>%
                         map(reduce, union),
                       ~filter(.x, variable %in% .y)),
           plot_title = paste0(globals$cohort_labs[[i]], ', ',
                               levels(sub_drugs$test[[i]]$consensusClass),
                               ' vs LumP, ',
                               sub_drugs$title_suffixes["common_significant"])) %>%
      pmap(plot_volcano,
           regulation_variable = 'biserial_r',
           p_variable = 'p_adjusted',
           regulation_level = 0.3,
           x_lab = 'effect size, biserial r, vs LumP',
           y_lab = expression('-log'[10] * ' pFDR'),
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   legend.title = element_markdown()),
           label_variable = 'title_label',
           top_significant = 20,
           txt_size = 2.5,
           txt_face = 'plain',
           label_type = 'text') %>%
      map(tag2subtitle) %>%
      map(~.x +
            scale_fill_manual(values = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray60'),
                              labels = c(upregulated = 'resistant',
                                         downregulated = 'sensitive',
                                         ns = 'ns'),
                              name = 'resistance vs LumP'))

  }

# Heat maps for the common significant effects and FGFR inhibitors ----

  insert_msg('Heat maps for the common signifiant effects and FGFR inhibitors')

  ## the common significant effects

  sub_drugs$hm_plots$common_significant <-
    list(data = sub_drugs$data,
         plot_title = globals$cohort_labs[names(sub_drugs$data)] %>%
           paste(sub_drugs$title_suffixes["common_significant"],
                 sep = ', ')) %>%
    pmap(heat_map,
         variables = unique(unlist(sub_drugs$common_significant)),
         split_fct = 'consensusClass',
         normalize = TRUE,
         x_lab = 'cancer sample',
         hide_x_axis_text = TRUE,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank(),
                 axis.text.y = element_markdown()),
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales::squish)

  sub_drugs$hm_plots$fgfr_drugs <-
    list(data = sub_drugs$data,
         plot_title = globals$cohort_labs[names(sub_drugs$data)] %>%
           paste(sub_drugs$title_suffixes["fgfr_drugs"],
                 sep = ', ')) %>%
    pmap(heat_map,
         variables = sub_drugs$fgfr_drugs,
         split_fct = 'consensusClass',
         normalize = TRUE,
         x_lab = 'cancer sample',
         hide_x_axis_text = TRUE,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank(),
                 axis.text.y = element_markdown()),
         midpoint = 0,
         limits = c(-3, 3),
         oob = scales::squish)

  ## additional styling
  ## color of the Y axis label codes for the training data set
  ## significant effects are highlighted in bold
  ## for the FGFR drugs we're showing the p values and effect sizes

  sub_drugs$hm_plots$common_significant <-
    map2(sub_drugs$hm_plots$common_significant,
         sub_drugs$kruskal,
         ~.x +
           scale_y_discrete(labels = set_names(.y$hm_label, .y$variable)))

  sub_drugs$hm_plots$fgfr_drugs <-
    map2(sub_drugs$hm_plots$fgfr_drugs,
         sub_drugs$kruskal,
         ~.x +
           scale_y_discrete(labels = set_names(.y$hm_label_long, .y$variable)))

# Box plots for the FGFR-targeting drugs -------

  insert_msg('Box plots for the FGFR-trageting drugs')

  ## plotting meta-data

  sub_drugs$box_data <- sub_drugs$kruskal %>%
    map(filter, variable %in% sub_drugs$fgfr_drugs)

  ## box plots

  for(i in names(sub_drugs$data)) {

    sub_drugs$fgfr_plots[[i]] <-
      list(variable = sub_drugs$box_data[[i]]$variable,
           plot_title = sub_drugs$box_data[[i]]$title_label,
           plot_subtitle = sub_drugs$box_data[[i]]$plot_cap,
           y_lab = sub_drugs$box_data[[i]]$variable %>%
             exchange(sub_drugs$lexicon,
                      value = 'axis_title')) %>%
      pmap(plot_variable,
           sub_drugs$data[[i]],
           split_factor = 'consensusClass',
           type = 'box',
           cust_theme = globals$common_theme,
           x_lab = 'MIBC consensus class',
           x_n_labs = TRUE) %>%
      map(~.x +
            scale_fill_manual(values = globals$sub_colors)) %>%
      set_names(sub_drugs$box_data[[i]]$variable)

  }

# Result summary table ------

  insert_msg('Result summary table')

  sub_drugs$result_tbl <-
    map2(sub_drugs$stats,
         sub_drugs$kruskal %>%
           map(~.x[c('variable', 'title_label', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    transmute(Cohort = globals$cohort_labs[cohort],
              Variable = title_label,
              Variable = ifelse(is.na(Variable), 'Samples, N', Variable),
              `Resistance metric` = exchange(variable,
                                             sub_drugs$lexicon,
                                             value = 'unit'),
              LumP = LumP,
              LumU = LumU,
              `Stroma-rich` = `Stroma-rich`,
              `Ba/Sq` = `Ba/Sq`,
              Significance = significance,
              `Effect size` = eff_size)

# END ---------

  sub_drugs$data <- NULL
  sub_drugs$variables <- NULL
  sub_drugs$fgfr_drugs <- NULL
  sub_drugs$posthoc_levels <- NULL
  sub_drugs$box_data <- NULL

  sub_drugs <- compact(sub_drugs)

  insert_tail()
