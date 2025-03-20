# Identification of mutation hotspots in FGFR genes

  insert_head()

# container --------

  expl_hot <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## detailed mutation data for FGFR genes and total numbers of samples
  ## and domain information

  expl_hot$data <- expl_globals$data %>%
    map(~.x[, c("gene_symbol",
                "protein_start_position",
                "mutation_type",
                "variant_type",
                "protein_domain",
                "protein_change")]) %>%
    map(filter, !is.na(protein_start_position))

  expl_hot$n_totals <- expl_globals$n_totals %>%
    compress(names_to = 'cohort',
             values_to = 'n_total')

  ## domain tables, schemes, and protein lengths

  expl_hot$domain_lst <- expl_globals$domain_lst

  expl_hot$domain_plots <- expl_fgfr$domain_plots %>%
    map(~.x + theme(plot.title = element_markdown()))

  expl_hot$receptor_lengths <- globals$receptor_lengths

  ## mutation category labs

  expl_hot$category_labs <-
    c(mutation_type = 'mutation\ntype',
      variant_type = 'variant\ntype',
      protein_domain = 'protein\ndomain')

# Counts of mutations per protein residue --------

  insert_msg('counts of mutations per protein residue')

  ## split by mutation type, variant type, and protein domain

  for(i in c('mutation_type', 'variant_type', 'protein_domain')) {

    expl_hot$frequency[[i]] <- expl_hot$data %>%
      map(count,
          gene_symbol,
          .data[[i]],
          protein_start_position,
          protein_change) %>%
      compress(names_to = 'cohort') %>%
      blast(gene_symbol) %>%
      map(left_join,
          expl_hot$n_totals,
                by = 'cohort') %>%
      map(mutate,
          cohort = factor(cohort, names(expl_hot$data)),
          percent = n/n_total * 100) %>%
      map(~filter(.x, complete.cases(.)))

  }

  ## protein domain stats, collapsing by protein changes

  expl_hot$frequency$protein_domain <-
    expl_hot$frequency$protein_domain %>%
    map(group_by,
        cohort,
        protein_domain,
        protein_start_position) %>%
    map(summarise,
        n = sum(n),
        n_total = mean(n_total),
        percent = n/n_total * 100,
        protein_change = paste(protein_change, collapse = '\n')) %>%
    map(ungroup) %>%
    map(arrange, -percent)

# Sliding window density of mutations --------

  insert_msg('Denstity of mutations')

  expl_hot$density <- expl_hot$data %>%
    map(blast, gene_symbol) %>%
    map2(., expl_hot$n_totals$n_total[names(.)],
         function(x, y) list(data = x,
                             n_total = y) %>%
           pmap(slide_window_density, k = 20))

  ## appending with the domain information

  expl_hot$density <- expl_hot$density %>%
    map(~map2(.x, expl_hot$domain_lst[names(.x)], append_domain_))

  ## merging the cohorts

  expl_hot$density <- expl_hot$density %>%
    transpose %>%
    map(compress, names_to = 'cohort') %>%
    map(left_join, expl_hot$n_totals, by = 'cohort') %>%
    map(mutate,
        cohort = factor(cohort, names(expl_hot$data)))

# Plots of mutation percentages and densities as a function of protein position -------

  insert_msg('Plots of mutation percentages and densities')

  ## plots of percentages

  for(i in names(expl_hot$frequency)) {

    expl_hot$percentage_plots[[i]] <-
      list(data = expl_hot$frequency[[i]],
           plot_title = names(expl_hot$frequency[[i]]) %>%
             html_italic %>%
             paste0(', mutation hotspots')) %>%
      pmap(plot_hot_spot,
           color_variable = i,
           color_palette = globals[c("mutation_colors",
                                     "variant_colors",
                                     "domain_colors")] %>%
             reduce(c),
           color_lab = expl_hot$category_labs[[i]],
           hot_cutoff = 1,
           point_shape = 16,
           point_size = 3,
           box.padding = 0.5,
           point.padding = 0.1,
           force = 2,
           direction = 'x',
           nudge_y = 6,
           segment.linetype = 'dashed',
           seed = 12345) %>%
      map(~.x +
            theme(plot.title = element_markdown()))

  }

  ## plots of mutation density

  expl_hot$density_plots <-
    list(data = expl_hot$density,
         plot_title = names(expl_hot$density) %>%
           html_italic %>%
           paste0(', mutation density')) %>%
    pmap(plot_hot_density) %>%
    map(~.x +
          theme(plot.title = element_markdown()) +
          scale_y_continuous(trans = 'sqrt'))

# Appeding the plots with domain schemes -------

  insert_msg('Appending the plots with domain schemes')

  for(i in names(expl_hot$percentage_plots)) {

    expl_hot$domain_percentage_plots[[i]] <-
      list(position_plot = expl_hot$percentage_plots[[i]],
           domain_tbl = expl_hot$domain_lst,
           protein_name = names(expl_hot$receptor_lengths),
           protein_length = expl_hot$receptor_lengths) %>%
      pmap(append_domain_scheme,
           rel_heights = c(1, 4.5),
           rm_domain_legend = TRUE)

  }

  expl_hot$domain_density_plots <-
    list(position_plot = expl_hot$density_plots,
         domain_tbl = expl_hot$domain_lst,
         protein_name = names(expl_hot$receptor_lengths),
         protein_length = expl_hot$receptor_lengths) %>%
    pmap(append_domain_scheme,
         rel_heights = c(1, 4.5),
         rm_domain_legend = TRUE)

# Versions of the plots for the paper figure (larger fonts) -------

  insert_msg('Appending the plots with domain schemes, paper versions')

  for(i in names(expl_hot$percentage_plots)) {

    expl_hot$domain_percentage_plots_paper[[i]] <-
      list(position_plot = expl_hot$percentage_plots[[i]],
           domain_tbl = expl_hot$domain_lst,
           protein_name = names(expl_hot$receptor_lengths),
           protein_length = expl_hot$receptor_lengths) %>%
      pmap(append_domain_scheme,
           rel_heights = c(1, 4.5),
           rm_domain_legend = TRUE,
           cust_theme = globals$figure_theme +
             theme(plot.title = element_markdown(size = 10),
                   axis.title.x = element_blank()))

  }

  expl_hot$domain_density_plots_paper <-
    list(position_plot = expl_hot$density_plots,
         domain_tbl = expl_hot$domain_lst,
         protein_name = names(expl_hot$receptor_lengths),
         protein_length = expl_hot$receptor_lengths) %>%
    pmap(append_domain_scheme,
         rel_heights = c(1, 4.5),
         rm_domain_legend = TRUE,
         cust_theme = globals$figure_theme +
           theme(plot.title = element_markdown(size = 10),
                 axis.title.x = element_blank()))

# END -------

  expl_hot <-
    expl_hot[c("frequency", "density",
               "percentage_plots", "density_plots",
               "domain_percentage_plots", "domain_density_plots",
               "domain_percentage_plots_paper",
               "domain_density_plots_paper")]

  insert_tail()
