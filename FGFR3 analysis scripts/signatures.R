# Differences of ssGSEA scores of FGFR-related gene signatures between
# cancers with and without FGFR3 mutations.
#
# Statistical significance is determined by Welch's two-tailed T test with
# Cohen's d effect size statistic.
# Differential gene expression in single cohorts is considered for pFDR < 0.05
# and d >= 0.5 (at least moderate regulation).
# Differences shared by at least two cohorts are investigeted in detail.

  insert_head()

# container ------

  fgfr_sig <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## signatures of interest and their lexicon

  fgfr_sig$lexicon <- sig$lexicon

  fgfr_sig$variables <- fgfr_sig$lexicon$variable

  ## signature SSGSEA scores

  fgfr_sig$data <- sig$scores

  fgfr_sig$data <-
    map2(fgfr_globals$genetic_groups[names(fgfr_sig$data)],
        fgfr_sig$data,
        inner_join, by = 'sample_id')

# descriptive stats -------

  insert_msg('Descriptive stats')

  fgfr_sig$stats <- fgfr_sig$data %>%
    map(fast_num_stats,
        variables = fgfr_sig$variables,
        split_fct = 'FGFR3_mutation')

# Testing for differences in ssGSEA scores between mutant and WT samples -----

  insert_msg('Testing for differences')

  ## Welch's T test

  fgfr_sig$test <- fgfr_sig$data %>%
    map(~f_t_test(.x[fgfr_sig$variables],
                  .x$FGFR3_mutation,
                  type = 'welch',
                  as_data_frame = TRUE,
                  safely = TRUE,
                  adj_method = 'BH')) %>%
    map(as_tibble)

  ## formatting, significant regulation

  fgfr_sig$test <- fgfr_sig$test %>%
    map(re_adjust) %>%
    map(p_formatter, text = TRUE) %>%
    map(mutate,
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(cohen_d >= 0.5, 'upregulated',
                                   ifelse(cohen_d <= -0.5,
                                          'downregulated', 'ns'))),
        regulation = factor(regulation,
                            c('upregulated', 'downregulated', 'ns')))

  ## formatting, plot labels

  fgfr_sig$test <- fgfr_sig$test %>%
    map(mutate,
        variable_label = exchange(variable, fgfr_sig$lexicon),
        eff_size = paste('d =', signif(cohen_d, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '),
        hm_label = paste(variable_label, plot_cap, sep = '<br>'),
        hm_label = ifelse(regulation %in% c('upregulated', 'downregulated'),
                          html_bold(hm_label), hm_label))

# Significant effects -----

  insert_msg('Significant effects')

  ## in single cohorts

  fgfr_sig$significant <- fgfr_sig$test %>%
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>%
    map(blast, regulation) %>%
    transpose %>%
    map(map, ~.x$variable)

  ## the common significant effects shared by at least two cohorts

  fgfr_sig$common_significant <- fgfr_sig$significant %>%
    map(shared_features, m = 2) %>%
    map(as.character)

# Volcano plots --------

  insert_msg('Volcano plots')

  ## for the common differentially regulated signatures

  fgfr_sig$volcano_plots <-
    list(data = fgfr_sig$test %>%
           map(filter, variable %in% reduce(fgfr_sig$common_significant, union)),
         plot_title = globals$cohort_labs[names(fgfr_sig$test)] %>%
           paste(', FGFR-related gene signatures')) %>%
    pmap(plot_volcano,
         regulation_variable = 'estimate',
         p_variable = 'p_adjusted',
         regulation_level = 0,
         x_lab = '\u0394 ssGSEA, mutated vs WT',
         y_lab = expression('-log'[10] * ' pFDR'),
         cust_theme = globals$common_theme +
           theme(plot.title = element_markdown(),
                 legend.title = element_markdown()),
         label_variable = 'variable_label',
         top_significant = 20,
         txt_size = 2.5,
         txt_face = 'plain',
         label_type = 'text') %>%
    map(~.x +
          labs(subtitle = paste('total signatures: n =',
                                length(fgfr_sig$variables)))) %>%
    map(tag2subtitle)

# Heat map representation ------

  insert_msg('Heat maps for the common significant effects')

  fgfr_sig$hm_plots <-
    list(data = fgfr_sig$data,
         plot_title = globals$cohort_labs[names(fgfr_sig$data)] %>%
           paste(', FGFR-related gene signatures')) %>%
    pmap(heat_map,
         variables = reduce(fgfr_sig$common_significant, union),
         split_fct = 'FGFR3_mutation',
         normalize = FALSE,
         x_lab = 'cancer sample',
         hide_x_axis_text = TRUE,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 strip.background.y = element_blank(),
                 strip.text.y = element_blank(),
                 axis.text.y = element_markdown()),
         limits = c(-0.75, 0.75),
         midpoint = 0,
         oob = scales::squish)

  fgfr_sig$hm_plots <-
    map2(fgfr_sig$hm_plots,
         fgfr_sig$test,
         ~.x +
           scale_y_discrete(labels = set_names(.y$hm_label,
                                               .y$variable)))

# Result summary -------

  insert_msg('Result summary')

  fgfr_sig$result_tbl <-
    map2(fgfr_sig$stats,
         fgfr_sig$test %>%
           map(p_formatter, text = FALSE) %>%
           map(~.x[c('variable', 'variable_label', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(format_res_tbl, dict = NULL) %>%
    compress(names_to = 'cohort') %>%
    transmute(Cohort = globals$cohort_labs[cohort],
              Variable = variable_label,
              Variable = ifelse(is.na(Variable), 'Samples, N', Variable),
              `FGFR WT` = WT,
              `FGFR mutated` = mutated,
              `FDR p value` = significance,
              `Effect size, Cohen's d` = eff_size)

# END -------

  fgfr_sig$lexicon <- NULL
  fgfr_sig$data <- NULL

  fgfr_sig <- compact(fgfr_sig)

  insert_tail()
