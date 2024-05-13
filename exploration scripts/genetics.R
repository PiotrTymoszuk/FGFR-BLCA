# Genetic alterations in the cohorts of interest

  insert_head()

# container -----

  expl_genet <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## data

  expl_genet$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    transpose %>%
    map(compact) %>%
    map(map, select, -sample_id)

  ## common variables

  expl_genet$variables <- expl_genet$data %>%
    map(map, names) %>%
    map(reduce, intersect)

  for(i in names(expl_genet$data)) {

    expl_genet$data[[i]] <- expl_genet$data[[i]] %>%
      map(select, all_of(expl_genet$variables[[i]])) %>%
      compress(names_to = 'cohort') %>%
      mutate(cohort = factor(cohort, globals$analysis_cohorts)) %>%
      relocate(cohort)

  }

  ## variable lexicon with ready-to-use axis labels in HTML format
  ## components of the FGF/FGFR system are highlighted with bold font

  expl_genet$fgf_genes <-
    list(receptors = globals$receptors,
         ligands = globals$ligands)

  expl_genet$lexicon <-
    tibble(variable = reduce(expl_genet$variables, union)) %>%
    mutate(label = html_italic(variable),
           label = ifelse(variable %in% reduce(expl_genet$fgf_genes, union),
                          html_bold(label), label))

# N numbers -------

  insert_msg('N numbers')

  expl_genet$n_numbers <- expl_genet$data %>%
    map(count, cohort) %>%
    map(mutate, cohort = globals$cohort_labs[cohort])

  ## captions and figure legends

  expl_genet$n_captions <- expl_genet$n_numbers %>%
    map(~map2_chr(.x[[1]], .x[[2]],
                  paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

  expl_genet$n_legends <- expl_genet$n_numbers %>%
    map(~map2_chr(.x[[1]], .x[[2]],
                  paste, sep = ': n = '))

# Frequency of alterations -----

  insert_msg('Frequency of alterations')

  expl_genet$stats <- expl_genet$data %>%
    map(bin2stats,
        split_factor = 'cohort') %>%
    map(mutate,
        cohort = factor(cohort, globals$analysis_cohorts),
        cohort = droplevels(cohort))

  ## top alterations present in at least 5% of samples
  ## in all cohorts

  expl_genet$top_alterations <- expl_genet$stats %>%
    map(blast, cohort) %>%
    map(map, filter, percent >= 5) %>%
    map(map, ~.x$variable) %>%
    map(reduce, intersect)

# Testing for differences between the cohorts ------

  insert_msg('Testing for differences between the cohorts')

  ## done only for the top most frequent alterations

  expl_genet$test <-
    map2(expl_genet$data,
         expl_genet$top_alterations,
         ~compare_variables(.x,
                            variables = .y,
                            split_factor = 'cohort',
                            what = 'eff_size',
                            types = 'cramer_v',
                            exact = FALSE,
                            ci = FALSE,
                            pub_styled = FALSE,
                            adj_method = 'BH')) %>%
    map(mutate,
        eff_size = paste('V =', signif(estimate, 2)),
        plot_cap = paste(eff_size, significance, sep = ', '))

  ## significant effects

  expl_genet$significant <- expl_genet$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

# Bar plots of frequencies of the most common alterations -------

  insert_msg('Bar plots')

  expl_genet$plots <-
    list(x = map2(expl_genet$stats,
                  expl_genet$top_alterations,
                  ~filter(.x, variable %in% .y)),
         y = c('Somatic mutations',
               'Gene deletion',
               'Gene amplification'),
         z = expl_genet$n_captions,
         w = expl_genet$n_legends) %>%
    pmap(function(x, y, z, w) x %>%
           ggplot(aes(x = percent,
                      y = reorder(variable, percent),
                      fill = cohort)) +
           geom_bar(stat = 'identity',
                    color = 'black',
                    position = position_dodge(0.9)) +
           scale_y_discrete(labels = function(x) exchange(x, expl_genet$lexicon)) +
           scale_fill_manual(values = globals$cohort_colors,
                             labels = w,
                             name = 'Cohort') +
           globals$common_theme +
           theme(axis.text.y = element_markdown(),
                 axis.title.y = element_blank()) +
           labs(title = y,
                subtitle = z,
                x = '% of cohort'))

# Result table --------

  insert_msg('Result table')

  ## for the top most frequent alterations

  expl_genet$result_tbl <-
    map2(expl_genet$stats,
         map(expl_genet$test, ~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map2(., expl_genet$top_alterations,
         ~filter(.x, variable %in% .y)) %>%
    map(mutate,
        cohort = globals$cohort_labs[cohort],
        percent = signif(percent, 2)) %>%
    map(arrange, variable, cohort) %>%
    map(select,
        variable, cohort, percent, n, n_total, significance, eff_size) %>%
    map(set_names,
        c('Gene symbol', 'Cohort',
          'Percentage of alterations',
          'Samples with alterations, N',
          'Total samples, N',
          'Significance',
          'Effect size'))

# Detailed bar plots for the genes of interest -------

  insert_msg('Detailed bar plots for the genes of interest')

  ## plotting data

  expl_genet$fgf_stats <- expl_genet$stats %>%
    map(filter,
        variable %in% reduce(expl_genet$fgf_genes, union)) %>%
    map(mutate,
        gene_group = ifelse(variable %in% expl_genet$fgf_genes$receptors,
                            'receptors', 'ligands'),
        gene_group = factor(gene_group, c('receptors', 'ligands'))) %>%
    map(blast, variable) %>%
    map(map_dfr, function(x) if(all(x$n == 0)) NULL else x)

  ## bar plots

  expl_genet$fgf_plots <-
    list(x = expl_genet$fgf_stats,
         y = paste0(c('Somatic mutations',
                      'Gene deletion',
                      'Gene amplification'),
                    ', <em>FGFR/FGF</em> genes'),
         z = expl_genet$n_captions,
         w = expl_genet$n_legends) %>%
    pmap(function(x, y, z, w) x %>%
           ggplot(aes(x = percent,
                      y = reorder(variable, percent),
                      fill = cohort)) +
           facet_grid(gene_group ~ .,
                      scales = 'free',
                      space = 'free') +
           geom_bar(stat = 'identity',
                    color = 'black',
                    position = position_dodge(0.9)) +
           scale_fill_manual(values = globals$cohort_colors,
                             labels = w,
                             name = 'Cohort') +
           globals$common_theme +
           theme(axis.text.y = element_text(face = 'italic'),
                 axis.title.y = element_blank(),
                 plot.title = element_markdown()) +
           labs(title = y,
                subtitle = z,
                x = '% of cohort'))

# Result table for the genes of interest -------

  insert_msg('Result table for the genes of interest')

  expl_genet$fgf_result_tbl <- expl_genet$fgf_stats %>%
    map(mutate,
        percent = signif(percent, 2),
        cohort = globals$cohort_labs[cohort]) %>%
    map(arrange, gene_group, variable, cohort) %>%
    map(select,
        gene_group, variable, cohort,
        percent, n, n_total) %>%
    map(set_names,
        c('Gene group', 'Gene symbol', 'Cohort',
          'Percentage of alterations',
          'Samples with alterations, N',
          'Total samples, N'))

# END -------

  expl_genet$data <- NULL
  expl_genet$fgf_stats <- NULL

  expl_genet <- compact(expl_genet)

  insert_tail()
