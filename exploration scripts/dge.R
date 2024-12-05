# Differential expression of FGF and FGFR genes in cancers with WT and mutated
# alleles of FGFR1 - 4.
#
# Differential gene expression is investigated by FDR-corrected two-tailed
# T test with Cohen's d effect size statistic. Differentially expressed genes
# are defined by pFDR < 0.05 and fold regulation of at least 1.5.

  insert_head()

# container -------

  expl_dge <- list()

# Analysis data -------

  insert_msg('Analysis data')

  ## gene expression variables

  expl_dge$variables <- globals[c("receptors", "ligands", "binding_proteins")]

  expl_dge$split_factors <- paste0('FGFR', 1:4, '_mut')

  expl_dge$split_factors <- set_names(expl_dge$split_factors,
                                      expl_dge$split_factors)

  expl_dge$split_labels <- paste0('FGFR', 1:4) %>%
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
    map(select, sample_id, all_of(paste0('FGFR', 1:4))) %>%
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

  ## data set-specific expression lists

  expl_dge$variables <- expl_dge$variables %>%
    reduce(union)

  expl_dge$variables <- expl_dge$data %>%
    map(names) %>%
    map(intersect, expl_dge$variables)

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
      map(fast_num_stats,
          data = expl_dge$data[[i]],
          variables = expl_dge$variables[[i]])

  }

  expl_dge$stats <- transpose(expl_dge$stats)

# Testing -------

  insert_msg('Testing')

  for(i in names(expl_dge$data)) {

    expl_dge$test[[i]] <- expl_dge$split_factors %>%
      map(~safely(f_t_test)(expl_dge$data[[i]][expl_dge$variables[[i]]],
                            f = expl_dge$data[[i]][[.x]],
                            type = 'welch',
                            as_data_frame = TRUE,
                            safely = TRUE,
                            adj_method = 'BH')) %>%
      map(~.x$result) %>%
      compact %>%
      map(as_tibble)

  }

  expl_dge$test <- expl_dge$test %>%
    transpose %>%
    map(compact)

# Formatting of the results ------

  insert_msg('Formatting of the testing results')

  expl_dge$test <- expl_dge$test %>%
    map(map, re_adjust, method = 'BH') %>%
    map(map,
        mutate,
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(estimate >= log2(1.25), 'upregulated',
                                   ifelse(estimate <= -log2(1.25),
                                          'downregulated', 'ns'))),
        regulation = factor(regulation,
                            c('upregulated', 'downregulated', 'ns')),
        eff_size = paste('d =', signif(cohen_d, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '))

# Differentially regulated genes --------

  insert_msg('Differentially regulated genes')

  for(i in names(expl_dge$test)) {

    expl_dge$significant[[i]] <- expl_dge$test[[i]] %>%
      map(filter, regulation %in% c('upregulated','downregulated')) %>%
      map(blast, regulation) %>%
      map(map, ~.x$variable)

  }

  expl_dge$significant <- expl_dge$significant %>%
    map(transpose)

  ## differentially regulated genes shared by at least two cohorts

  expl_dge$common_significant <- expl_dge$significant %>%
    map(map, shared_features, m = 2) %>%
    map(map, as.character)

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
           plot_subtitle = expl_dge$n_numbers[[i]][names(expl_dge$test[[i]])] %>%
             map(~paste0('WT: n = ', .x$n[1],
                         ', mutated: n = ', .x$n[2]))) %>%
      pmap(plot_volcano,
           regulation_variable = 'estimate',
           p_variable = 'p_adjusted',
           regulation_level = log2(1.25),
           x_lab = expression('log'[2] * ' fold-regulation, mutated vs WT'),
           y_lab = expression('-log'[10] * ' pFDR'),
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   plot.tag = element_blank()),
           label_variable = 'variable',
           top_significant = 20,
           txt_size = 2.5,
           txt_face = 'italic',
           label_type = 'text')

  }

# Box plots for FGFR3 mutations ---------

  insert_msg('Box plots')

  for(i in names(expl_dge$data)) {

    expl_dge$box_plots[[i]] <-
      list(variable = expl_dge$test$FGFR3_mut[[i]]$variable,
           plot_title = expl_dge$test$FGFR3_mut[[i]]$variable %>%
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
      set_names(expl_dge$test$FGFR3_mut[[i]]$variable)

  }

# Result table --------

  insert_msg('Result tables')

  for(i in names(expl_dge$stats)) {

    expl_dge$result_tbl[[i]] <-
      map2(expl_dge$stats[[i]][names(expl_dge$test[[i]])],
           map(expl_dge$test[[i]],
               ~.x[c('variable', 'significance', 'eff_size')]),
           left_join, by = 'variable') %>%
      map(format_res_tbl,
          remove_complete = TRUE,
          dict = NULL) %>%
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
