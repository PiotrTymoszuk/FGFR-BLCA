# Details of mutation type and mutation location for FGFR1 - 4.
# The analysis is done for the GENIE, MSK, TCGA, and BCAN cohorts.

  insert_head()

# container -----

  expl_fgfr <- list()

# analysis data -------

  insert_msg('Analysis data')

  ## detailed data and numbers of complete cases

  expl_fgfr$data <- expl_globals$data

  expl_fgfr$n_totals <- expl_globals$n_totals

  ## domain information

  expl_fgfr$domain_lst <- expl_globals$domain_lst

# tallying the mutation types and domains -------

  insert_msg('Frequencies of mutation and variant types, and domains')

  for(i in c('mutation_type', 'variant_type', 'protein_domain')) {

    expl_fgfr$stats[[i]] <- expl_fgfr$data %>%
      map(blast, gene_symbol) %>%
      map(map, count, .data[[i]]) %>%
      map(map, arrange, desc(.data[[i]])) %>%
      map(map, mutate, plot_pos = cumsum(n) - 0.5 * n) %>%
      map(compress, names_to = 'gene_symbol') %>%
      map2(., expl_fgfr$n_totals,
           ~mutate(.x,
                   n_total = .y,
                   percent = n/n_total * 100,
                   plot_pos = plot_pos/n_total * 100))

  }

# Stack plots with mutation and variant types -------

  insert_msg('Stack plots with mutation and variant types')

  expl_fgfr$stack_plots <-
    list(data = unlist(expl_fgfr$stats, recursive = FALSE),
         fill_var = c(rep('mutation_type', 4),
                      rep('variant_type', 4),
                      rep('protein_domain', 4)),
         plot_title = c(paste('Mutation type,',
                              globals$cohort_labs[c("genie", "msk", "tcga", "bcan")]),
                        paste('Variant type,',
                              globals$cohort_labs[c("genie", "msk", "tcga", "bcan")]),
                        paste('Protein domain,',
                              globals$cohort_labs[c("genie", "msk", "tcga", "bcan")])),
         plot_subtitle = paste('total samples: n =',
                               rep(expl_fgfr$n_totals, 3)),
         fill_title = '',
         palette = globals[c(rep("mutation_colors", 4),
                             rep("variant_colors", 4),
                             rep("domain_colors", 4))],
         reverse_fill = c(rep(TRUE, 8), rep(FALSE, 4))) %>%
    pmap(plot_stack,
         x_var = 'percent',
         y_var = 'gene_symbol',
         y_limits = paste0('FGFR', c(1, 4, 2, 3))) %>%
    map(~.x +
          theme(axis.text.y = element_text(face = 'italic')))

# Position plots ---------

  insert_msg('Position plots')

  expl_fgfr$position_data <- expl_fgfr$data %>%
    map(blast, gene_symbol)

  for(i in names(expl_fgfr$position_data)) {

    expl_fgfr$position_plots[[i]] <-
      list(data = expl_fgfr$position_data[[i]],
           plot_title = expl_fgfr$position_data[[i]] %>%
             names %>%
             html_italic %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           protein_length = globals$receptor_lengths) %>%
      pmap(position_plot,
           fill_var = 'protein_domain',
           palette = globals$domain_colors[c("Ig-like C2-type 1",
                                             "Ig-like C2-type 2",
                                             "Ig-like C2-type 3",
                                             "Protein kinase")] %>%
             c('not assigned' = 'gray40'),
           y_limits = rev(c('missense',
                            'nonsense',
                            'in−frame deletion',
                            'in−frame insertion',
                            'frame shift deletion',
                            'frame shift insertion',
                            'splice region',
                            'splice site')),
           point_shape = 16,
           point_size = 1,
           point_alpha = 1)%>%
      map(~.x +
            theme(plot.title = element_markdown()) +
            labs(fill = 'Protein domain'))

  }

# Domain schemes -------

  insert_msg('Domain schemes')

  expl_fgfr$domain_plots <-
    list(domain_tbl = expl_fgfr$domain_lst,
         protein_name = names(expl_fgfr$domain_lst),
         protein_length = globals$receptor_lengths) %>%
    pmap(plot_domains)

  ## adding the domain schemes to the position plots

  for(i in names(expl_fgfr$position_plots)) {

    expl_fgfr$position_plots[[i]] <-
      list(position_plot = expl_fgfr$position_plots[[i]],
           domain_tbl = expl_fgfr$domain_lst,
           protein_name = names(expl_fgfr$domain_lst),
           protein_length = globals$receptor_lengths) %>%
      pmap(append_domain_scheme,
           align = 'v',
           rm_position_legend = TRUE,
           rm_domain_legend = TRUE)

  }

# Result tables --------

  insert_msg('Result tables')

  ## frequency of mutations split by type, variant type and protein domain

  expl_fgfr$result_tbl <- expl_fgfr$stats %>%
    map(compress, names_to = 'cohort') %>%
    map(mutate, cohort = globals$cohort_labs[cohort]) %>%
    map2(., names(.),
         ~select(.x, cohort, gene_symbol, all_of(.y), percent, n, n_total)) %>%
    map2(., names(.),
         ~arrange(.x, cohort, gene_symbol, .data[[.y]])) %>%
    map2(., c('Mutation type',
              'Variant type',
              'Protein domain'),
         ~set_names(.x,
                    c('Cohort', 'Gene symbol', .y,
                      'Percentage of samples',
                      'Samples with mutation, N',
                      'Total samples, N')))

# END ------

  rm(i)

  expl_fgfr$position_data <- NULL

  expl_fgfr <- compact(expl_fgfr)

  insert_tail()
