# Differential expression of FGF and FGFR genes in cancers with WT and mutated
# alleles of FGFR1 - 4.
#
# Differential gene expression is investigated by FDR-corrected two-tailed
# T test with Cohen's d effect size statistic. Differentially expressed genes
# are defined by pFDR < 0.05 and fold regulation of at least 1.5.

  insert_head()

# container -------

  expl_dge <- list()

# parallel backend ------

  insert_msg('Parallel backend')

  plan('multisession')

# Analysis data -------

  insert_msg('Analysis data')

  ## gene expression variables

  expl_dge$variables <- globals[c("receptors", "ligands")]

  expl_dge$split_factors <- expl_dge$variables$receptors %>%
    paste0('_mut')

  expl_dge$split_labels <- expl_dge$variables$receptors %>%
    html_italic %>%
    set_names(expl_dge$split_factors)

  ## expression data

  expl_dge$expression <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(select, sample_id, any_of(reduce(expl_dge$variables, union)))

  expl_dge$variables <- expl_dge$variables %>%
    map(~.x[.x %in% names(expl_dge$expression[[1]])])

  ## mutation data

  expl_dge$mutation <- globals$cohort_expr %>%
    eval %>%
    map(~.x$mutation) %>%
    map(select, sample_id, all_of(expl_dge$variables$receptors)) %>%
    map(set_names, c('sample_id', expl_dge$split_factors))

  for(i in names(expl_dge$mutation)) {

    expl_dge$mutation[[i]][expl_dge$split_factors] <-
      expl_dge$mutation[[i]][expl_dge$split_factors] %>%
      map_dfc(car::recode, "0 = 'WT'; 1 = 'mutated'") %>%
      map_dfc(factor, c('WT', 'mutated'))

  }

  ## merging the expression and mutation data

  expl_dge$data <-
    map2(expl_dge$mutation[names(expl_dge$expression)],
         expl_dge$expression,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

  expl_dge$variables <- expl_dge$variables %>%
    reduce(union)

# N numbers ------

  insert_msg('N numbers')

  for(i in names(expl_dge$data)) {

    expl_dge$n_numbers[[i]] <- expl_dge$split_factors %>%
      map(~count(expl_dge$data[[i]], .data[[.x]])) %>%
      set_names(expl_dge$split_factors)

  }

  expl_dge$n_numbers <- transpose(expl_dge$n_numbers)

# Descriptive stats -------

  insert_msg('Descriptive stats')

  for(i in names(expl_dge$data)) {

    expl_dge$stats[[i]] <- expl_dge$split_factors %>%
      future_map(~explore(expl_dge$data[[i]],
                          variables = expl_dge$variables,
                          split_factor = .x,
                          what = 'table',
                          pub_styled = TRUE),
                 .options = furrr_options(seed = TRUE)) %>%
      map(format_desc) %>%
      set_names(expl_dge$split_factors)

  }

  expl_dge$stats <- transpose(expl_dge$stats)

  plan('sequential')

# Testing -------

  insert_msg('Testing')

  for(i in names(expl_dge$data)) {

    expl_dge$test[[i]] <- expl_dge$split_factors %>%
      map(~safely(test_two_groups)(expl_dge$data[[i]],
                                   variables = expl_dge$variables,
                                   split_fct = .x,
                                   type = 't',
                                   adj_method = 'BH')) %>%
      map(~.x$result) %>%
      set_names(expl_dge$split_factors)

  }

  expl_dge$test <- expl_dge$test %>%
    transpose %>%
    map(compact)

# Formatting of the results ------

  insert_msg('Formatting of the testing results')

  expl_dge$test <- expl_dge$test %>%
    map(compact) %>%
    map(map, re_adjust, method = 'none') %>%
    map(map,
        mutate,
        variable = response,
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(estimate >= log2(1.5), 'upregulated',
                                   ifelse(estimate <= -log2(1.5),
                                          'downregulated', 'ns'))),
        regulation = factor(regulation,
                            c('upregulated', 'downregulated', 'ns')),
        eff_size = paste('d =', signif(effect_size, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '))

# Differentially regulated genes --------

  insert_msg('Differentially regulated genes')

  for(i in names(expl_dge$test)) {

    expl_dge$significant[[i]] <- expl_dge$test[[i]] %>%
      map(filter, regulation %in% c('upregulated','downregulated')) %>%
      map(blast, regulation) %>%
      map(map, ~.x$response)

  }

  expl_dge$significant <- expl_dge$significant %>%
    map(transpose)

  ## differentially regulated genes shared by both cohorts

  expl_dge$common_significant <- expl_dge$significant %>%
    map(map, reduce, intersect)

# Volcano plots ------

  insert_msg('Volcano plots')

  for(i in names(expl_dge$test)) {

    expl_dge$volcano_plots[[i]] <-
      list(data = expl_dge$test[[i]],
           plot_title = expl_dge$test[[i]] %>%
             names %>%
             globals$cohort_labs[.] %>%
             paste(expl_dge$split_labels[i], .,
                   sep = ' mutations, '),
           plot_subtitle = expl_dge$n_numbers[[i]] %>%
             map(~paste0('WT: n = ', .x$n[1],
                         ', mutated: n = ', .x$n[2]))) %>%
      pmap(plot_volcano,
           regulation_variable = 'estimate',
           p_variable = 'p_adjusted',
           regulation_level = log2(1.5),
           x_lab = expression('log'[2] * ' fold-regulation, mutated vs WT'),
           y_lab = expression('-log'[10] * ' pFDR'),
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   plot.tag = element_blank()),
           label_variable = 'response',
           top_significant = 20,
           txt_size = 2.5,
           txt_face = 'italic',
           label_type = 'text')

  }

# Box plots for FGFR3 mutations ---------

  insert_msg('Box plots')

  for(i in names(expl_dge$data)) {

    expl_dge$box_plots[[i]] <-
      list(variable = expl_dge$test$FGFR3_mut[[i]]$response,
           plot_title = expl_dge$test$FGFR3_mut[[i]]$response %>%
             html_italic %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = expl_dge$test$FGFR3_mut[[i]]$plot_cap) %>%
      pmap(plot_variable,
           expl_dge$data[[i]],
           split_factor = 'FGFR3_mut',
           type = 'box',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   axis.title.x = element_markdown(),
                   legend.position = 'none'),
           x_n_labs = TRUE,
           x_lab = html_italic('FGFR3'),
           y_lab = expression('log'[2] * ' expression')) %>%
      map(~.x +
            scale_fill_manual(values = globals$mut_status_color)) %>%
      set_names(expl_dge$test$FGFR3_mut[[i]]$response)

  }

# Result table --------

  insert_msg('Result tables')

  for(i in names(expl_dge$stats)) {

    expl_dge$result_tbl[[i]] <-
      map2(expl_dge$stats[[i]],
           map(expl_dge$test[[i]],
               ~.x[c('variable', 'significance', 'eff_size')]),
           left_join, by = 'variable') %>%
      map(format_res_tbl,
          remove_complete = TRUE,
          dict = NULL) %>%
      map2(expl_dge$n_numbers[[i]], .,
           ~full_rbind(tibble(variable = 'Samples, N',
                              WT = .x$n[1],
                              mutated = .x$n[2]),
                       .y)) %>%
      map(mutate,
          mutated = as.character(mutated)) %>%
      compress(names_to = 'cohort') %>%
      mutate(cohort = globals$cohort_labs[cohort]) %>%
      select(cohort, variable,
             WT, mutated,
             significance, eff_size) %>%
      set_names(c('Cohort', 'Variable',
                  'Wild-type', 'mutated',
                  'Significance', 'Effect size'))

  }

# END ------

  rm(i)

  expl_dge$expression <- NULL
  expl_dge$mutation <- NULL

  expl_dge <- compact(expl_dge)

  insert_tail()
