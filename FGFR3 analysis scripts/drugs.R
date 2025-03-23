# Effect of FGFR3 mutations on predicted resistance to anti-cancer drugs
# in the TCGA, IMvigor, and BCAN cohorts.
#
# log IC50 and AUC values are compared between samples with and without
# FGFR3 mutation by Mann-Whitney test with biserial r effect size statistics.
# The p values are corrected for multiple testing with the FDR method.
# Significant effects of at least moderate effect size (pFDR < 0.05 and r >= 0.3)
# are of interests.
# Drugs found to differ in the predicted resistance between the FGFR3 mutation
# strata in at least two cohorts are analyzed in more detail.

  insert_head()

# container -------

  fgfr_drugs <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## drug variables and their lexicon

  fgfr_drugs$lexicon <- drugs$lexicon %>%
    mutate(axis_title = paste('predicted resistance,', unit))

  fgfr_drugs$variables <- fgfr_drugs$lexicon$variable

  ## drugs that target FGFR

  fgfr_drugs$fgfr_drugs <- drugs$fgfr_drugs

  ## the predicted resistance metrics, appended with the FGFR3 mutation status

  fgfr_drugs$data <- drugs$predictions

  fgfr_drugs$data <-
    map2(fgfr_globals$genetic_groups[names(fgfr_drugs$data)],
         fgfr_drugs$data,
         inner_join, by = 'sample_id')

  ## title suffixes

  fgfr_drugs$title_suffixes <-
    c(all = 'all drugs',
      common_significant = 'significant in at least two cohorts',
      fgfr_drugs = 'FGFR-targeting drugs')

# Descriptive stats: median, interquartile ranges, and ranges -------

  insert_msg('Descriptive stats')

  fgfr_drugs$stats <- fgfr_drugs$data %>%
    map(fast_num_stats,
        split_fct = 'FGFR3_mutation',
        variables = fgfr_drugs$variables) %>%
    map(mutate)

# Testing -------

  insert_msg('Testing')

  fgfr_drugs$test <- fgfr_drugs$data %>%
    map(~f_wilcox_test(.x[fgfr_drugs$variables],
                       .x$FGFR3_mutation,
                       exact = FALSE,
                       safely = TRUE,
                       as_data_frame = TRUE,
                       adj_method = 'BH')) %>%
    map(as_tibble)

  ## formatting the testing results

  fgfr_drugs$test <- fgfr_drugs$test %>%
    map(re_adjust) %>%
    map(mutate,
        eff_size = paste('r =', signif(biserial_r, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '),
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(biserial_r >= 0.3, 'resistant',
                                   ifelse(biserial_r <= -0.3,
                                          'sensitive', 'ns'))),
        regulation = factor(regulation,
                            c('resistant', 'sensitive', 'ns')),
        title_label = exchange(variable,
                               fgfr_drugs$lexicon,
                               value = 'plot_title'),
        hm_label = exchange(variable,
                            fgfr_drugs$lexicon,
                            value = 'hm_label'),
        hm_label_long = paste(title_label, plot_cap, sep = '<br>'),
        hm_label = ifelse(regulation %in% c('resistant', 'sensitive'),
                          html_bold(hm_label),
                          hm_label),
        hm_label_long = ifelse(regulation %in% c('resistant', 'sensitive'),
                               html_bold(hm_label_long),
                               paste0('<span style = "color: #585858">',
                                      hm_label_long, '</span>')))

# significant differences in predicted drug resistance -------

  insert_msg('significant differences in predicted resistance')

  fgfr_drugs$significant <- fgfr_drugs$test %>%
    map(filter, regulation %in% c('resistant', 'sensitive')) %>%
    map(blast, regulation) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## common significant effects: differences in at least two cohorts

  fgfr_drugs$common_significant <- fgfr_drugs$significant %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Volcano plots --------

  insert_msg('volcano plots')

  ## the plotting data

  fgfr_drugs$volcano_plots <-
    list(all = fgfr_drugs$test,
         common_significant = fgfr_drugs$test %>%
           map(filter,
               variable %in% reduce(fgfr_drugs$common_significant, union)),
         fgfr_drugs = fgfr_drugs$test %>%
           map(filter,
               variable %in% fgfr_drugs$fgfr_drugs))

  ## plots

  for(i in names(fgfr_drugs$volcano_plots)) {

    fgfr_drugs$volcano_plots[[i]] <-
      list(data = fgfr_drugs$volcano_plot[[i]],
           plot_title = globals$cohort_labs[names(fgfr_drugs$volcano_plots[[i]])] %>%
             paste(fgfr_drugs$title_suffixes[[i]], sep = ', ')) %>%
      pmap(plot_volcano,
           regulation_variable = 'biserial_r',
           p_variable = 'p_adjusted',
           regulation_level = 0.3,
           x_lab = 'effect size, biserial r, mutated vs WT',
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
                              name = paste('resistance,', html_italic('FGFR3'),
                                           'mutated vs WT')))

  }

# Heat maps of resistance metrics, common and FGFR-targeting drugs -------

  insert_msg('Heat maps for the common significant and FGFR-targetting drugs')

  fgfr_drugs$hm_plots$common_significant <-
    list(data = fgfr_drugs$data,
         plot_title = globals$cohort_labs[names(fgfr_drugs$data)] %>%
           paste(fgfr_drugs$title_suffixes["common_significant"],
                 sep = ', ')) %>%
    pmap(heat_map,
         variables = reduce(fgfr_drugs$common_significant, union),
         split_fct = 'FGFR3_mutation',
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

  fgfr_drugs$hm_plots$fgfr_drugs <-
    list(data = fgfr_drugs$data,
         plot_title = globals$cohort_labs[names(fgfr_drugs$data)] %>%
           paste(fgfr_drugs$title_suffixes["fgfr_drugs"],
                 sep = ', ')) %>%
    pmap(heat_map,
         variables = fgfr_drugs$fgfr_drugs,
         split_fct = 'FGFR3_mutation',
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

  fgfr_drugs$hm_plots$common_significant <-
    map2(fgfr_drugs$hm_plots$common_significant,
         fgfr_drugs$test,
         ~.x +
           scale_y_discrete(labels = set_names(.y$hm_label, .y$variable)))

  fgfr_drugs$hm_plots$fgfr_drugs <-
    map2(fgfr_drugs$hm_plots$fgfr_drugs,
         fgfr_drugs$test,
         ~.x +
           scale_y_discrete(labels = set_names(.y$hm_label_long, .y$variable)))

# Box plots for the FGFR-targeting drugs --------

  insert_msg('Box plots for the anti-FGFR drugs')

  ## plotting meta-data

  fgfr_drugs$box_data <- fgfr_drugs$test %>%
    map(filter, variable %in% fgfr_drugs$fgfr_drugs)

  ## box plots

  for(i in names(fgfr_drugs$data)) {

    fgfr_drugs$fgfr_plots[[i]] <-
      list(variable = fgfr_drugs$box_data[[i]]$variable,
           plot_title = fgfr_drugs$box_data[[i]]$title_label,
           plot_subtitle = fgfr_drugs$box_data[[i]]$plot_cap,
           y_lab = fgfr_drugs$box_data[[i]]$variable %>%
             exchange(fgfr_drugs$lexicon,
                      value = 'axis_title')) %>%
      pmap(plot_variable,
           fgfr_drugs$data[[i]],
           split_factor = 'FGFR3_mutation',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(axis.title.x = element_markdown()),
           x_lab = html_italic('FGFR3'),
           x_n_labs = TRUE) %>%
      map(~.x +
            scale_fill_manual(values = globals$mut_status_color)) %>%
      set_names(fgfr_drugs$box_data[[i]]$variable)

  }

# Result summary table --------

  insert_msg('Result summary')

  fgfr_drugs$result_tbl <-
    map2(fgfr_drugs$stats,
         fgfr_drugs$test %>%
           map(~.x[c('variable', 'title_label', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    compress(names_to = 'cohort') %>%
    transmute(Cohort = globals$cohort_labs[cohort],
              Variable = title_label,
              Variable = ifelse(is.na(Variable), 'Samples, N', Variable),
              `Resistance metric` = exchange(variable,
                                             fgfr_drugs$lexicon,
                                             value = 'unit'),
              `FGFR WT` = WT,
              `FGFR mutated` = mutated,
              Significance = significance,
              `Effect size` = eff_size)

# END ---------

  fgfr_drugs$data <- NULL
  fgfr_drugs$variables <- NULL
  fgfr_drugs$fgfr_drugs <- NULL
  fgfr_drugs$title_suffixes <- NULL
  fgfr_drugs$box_data <- NULL

  fgfr_drugs <- compact (fgfr_drugs)

  insert_tail()
