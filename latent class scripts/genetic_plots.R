# Plots for results of the genetic alteration enrichment analyses
# in Genetic Clusters. We're concentrating on genes found to be significantly
# enriched in the corresponding clusters of all three cohorts, i.e. GENIE, MSK,
# and TCGA BLCA.

  insert_head()

# container -------

  lca_genplots <- list()

# Analysis globals ------

  insert_msg('Analysis globals')

  ## cluster and cohort levels

  lca_genplots$clust_levels <- levels(lca_pred$assignment[[1]]$clust_id)

  lca_genplots$cohort_levels <- names(lca_pred$assignment)

  lca_genplots$alteration_levels <- c('mutation', 'deletion', 'amplification')

  ## colors for alteration types

  lca_genplots$alteration_colors <-
    c(mutation = globals$mut_status_color[["mutated"]],
      deletion = globals$del_status_color[["deleted"]],
      amplification = globals$amp_status_color[["amplified"]])

  ## numbers of significant gene mutations, deletions and amplification

  lca_genplots$n_genes <- list(mutation = lca_mut,
                               deletion = lca_del,
                               amplification = lca_amp) %>%
    map(~.x$enriched) %>%
    map(map, map_dbl, length) %>%
    map(map,
        compress,
        names_to = 'cohort',
        values_to = 'n') %>%
    map(compress,
        names_to = 'clust_id') %>%
    compress(names_to = 'alteration') %>%
    mutate(clust_id = factor(clust_id, lca_genplots$clust_levels),
           cohort = factor(cohort, lca_genplots$cohort_levels),
           alteraton = factor(alteration, lca_genplots$alteration_levels))

  ## significant genetic alterations shared by all three cohorts

  lca_genplots$common_enriched <- list(mutation = lca_mut,
                                       deletion = lca_del,
                                       amplification = lca_amp) %>%
    map(~.x$enriched) %>%
    map(map, reduce, intersect) %>%
    map(reduce, union)

  ## binary data for the common significant alterations

  lca_genplots$data <- list(genie = genie,
                            msk = msk,
                            tcga = tcga) %>%
    map(~.x[c('mutation', 'deletion', 'amplification')])

  for(i in names(lca_genplots$data)) {

    lca_genplots$data[[i]] <-
      map2(lca_genplots$data[[i]],
           lca_genplots$common_enriched,
           ~select(.x, sample_id, all_of(.y))) %>%
      map(column_to_rownames, 'sample_id') %>%
      map2(., names(.),
           ~set_names(.x, paste(names(.x), .y, sep = '_'))) %>%
      reduce(cbind) %>%
      rownames_to_column('sample_id') %>%
      inner_join(lca_pred$assignment[[i]],
                 by = 'sample_id')

  }

  ## variables of interest and their lexicon

  lca_genplots$variables <- lca_genplots$data[[1]] %>%
    select(-sample_id, -clust_id) %>%
    names(.)

  lca_genplots$lexicon <-
    tibble(variable = lca_genplots$variables) %>%
    mutate(gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           alteration = factor(alteration,
                               lca_genplots$alteration_levels))

  ## enrichment status for particular clusters

  lca_genplots$common_enriched <- list(mutation = lca_mut,
                                       deletion = lca_del,
                                       amplification = lca_amp) %>%
    map(~.x$enriched) %>%
    map(map, reduce, intersect) %>%
    map2(., names(.),
         function(x, y) x %>%
           map(~paste(.x, y, sep = '_'))) %>%
    transpose %>%
    map(reduce, union) %>%
    map(~.x[!stri_detect(.x, regex = '^_')]) %>%
    compact

  ## testing results for the global differences between the clusters

  lca_genplots$test <- list(mutation = lca_mut,
                            deletion = lca_del,
                            amplification = lca_amp) %>%
    map(~.x$test) %>%
    map2(., names(.),
         function(x, y) x %>%
           map(mutate,
               gene_symbol = variable,
               variable = paste(variable = paste(variable, y, sep = '_')),
               alteration = y,
               plot_title = paste(html_italic(gene_symbol), y, sep = ', '),
               plot_cap = paste(global_eff_size, global_significance, sep = ', '))) %>%
    transpose %>%
    map(reduce, rbind)

  lca_genplots$test <- lca_genplots$test %>%
    map(filter,
        variable %in% lca_genplots$variables,
        !duplicated(variable)) %>%
    map(select, variable, plot_title, plot_cap)

# Plot of numbers of significant alterations -------

  insert_msg('Plot of numbers of significant alterations')

  lca_genplots$n_plot <- lca_genplots$n_genes %>%
    ggplot(aes(x = n,
               y = factor(clust_id, rev(levels(clust_id))),
               fill = alteration)) +
    facet_grid(alteration ~ cohort,
               labeller = labeller(.cols = globals$cohort_labs)) +
    geom_bar(stat = 'identity',
             color = 'white') +
    geom_vline(xintercept = 0,
               linetype = 'solid') +
    #scale_x_continuous(trans = 'pseudo_log') +
    scale_fill_manual(values = lca_genplots$alteration_colors) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Numbers of genetic alterations',
         x = '# genes')

# Stack plots and oncoplots for the common significant alterations in particular clusters -------

  insert_msg('Stack plots and oncoplots for single clusters')

  for(i in names(lca_genplots$data)) {

    ## stack plots

    lca_genplots$cluster_stack_plots[[i]] <-
      list(variables = lca_genplots$common_enriched,
           plot_title = names(lca_genplots$common_enriched) %>%
             paste(globals$cohort_labs[[i]], sep = ', ')) %>%
      pmap(plot_bistack,
           data = lca_genplots$data[[i]],
           split_fct = 'clust_id',
           scale = 'percentage',
           x_lab = '% of cluster',
           cust_theme = globals$common_theme +
             theme(axis.text.y = element_text(face = 'italic'),
                   axis.title.y = element_blank()),
           color_scale = c("gray85", "orangered3"),
           variable_classification = lca_genplots$lexicon[, c("variable",
                                                              "alteration")]) %>%
      map(~.x +
            scale_y_discrete(labels = function(x) exchange(x,
                                                           lca_genplots$lexicon,
                                                           value = 'gene_symbol')))

    ## oncoplots

    lca_genplots$cluster_onco_plots[[i]] <-
      list(variables = lca_genplots$common_enriched,
           plot_title = names(lca_genplots$common_enriched) %>%
             paste(globals$cohort_labs[[i]], sep = ', ')) %>%
      pmap(plot_bionco,
           data = lca_genplots$data[[i]],
           split_fct = 'clust_id',
           x_lab = 'cancer sample',
           one_plot = FALSE,
           hide_x_axis_text = TRUE,
           cust_theme = globals$common_theme +
             theme(axis.title.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_text(angle = 90,
                                               hjust = 0)),
           color_scale = c("gray85", "orangered3"),
           variable_classification = lca_genplots$lexicon[, c("variable",
                                                              "alteration")]) %>%
      map(~.x$main +
            theme(axis.text.y = element_text(face = 'italic')) +
            scale_y_discrete(labels = function(x) exchange(x,
                                                           lca_genplots$lexicon,
                                                           value = 'gene_symbol')))

  }

# Stack plots and oncoplots for the entire cohorts -------

  insert_msg('Stack plots and oncoplots for the entire cohorts')

  lca_genplots$stack_plots <-
    list(data = lca_genplots$data,
         plot_title = globals$cohort_labs[names(lca_genplots$data)]) %>%
    pmap(plot_bistack,
         split_fct = 'clust_id',
         variables = lca_genplots$variables,
         scale = 'percentage',
         x_lab = '% of cluster',
         cust_theme = globals$common_theme +
           theme(axis.text.y = element_text(face = 'italic'),
                 axis.title.y = element_blank()),
         color_scale = c("gray85", "orangered3"),
         variable_classification = lca_genplots$lexicon[, c("variable",
                                                            "alteration")]) %>%
    map(~.x +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         lca_genplots$lexicon,
                                                         value = 'gene_symbol')))

  lca_genplots$onco_plots <-
    list(data = lca_genplots$data,
         plot_title = globals$cohort_labs[names(lca_genplots$data)]) %>%
    pmap(plot_bionco,
         variable = lca_genplots$variables,
         split_fct = 'clust_id',
         x_lab = 'cancer sample',
         one_plot = FALSE,
         hide_x_axis_text = TRUE,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.ticks.x = element_blank(),
                 strip.text.x = element_text(angle = 90,
                                             hjust = 0)),
         color_scale = c("gray85", "orangered3"),
         variable_classification = lca_genplots$lexicon[, c("variable",
                                                            "alteration")]) %>%
    map(~.x$main +
          theme(axis.text.y = element_text(face = 'italic')) +
          scale_y_discrete(labels = function(x) exchange(x,
                                                         lca_genplots$lexicon,
                                                         value = 'gene_symbol')))

# Stack plots for single genes --------

  insert_msg('Stack plots for single genes')

  plan('multisession')

  for(i in names(lca_genplots$data)) {

    lca_genplots$plots[[i]] <-
      list(variable = lca_genplots$test[[i]]$variable,
           plot_title = lca_genplots$test[[i]]$plot_title %>%
             paste(globals$cohort_labs[i], sep = ', '),
           plot_subtitle = lca_genplots$test[[i]]$plot_cap) %>%
      future_pmap(plot_variable,
                  lca_genplots$data[[i]],
                  split_factor = 'clust_id',
                  scale = 'percent',
                  type = 'stack',
                  cust_theme = globals$common_theme +
                    theme(plot.title = element_markdown()),
                  x_lab = 'genetic cluster',
                  y_lab = '% of cluster',
                  x_n_labs = TRUE,
                  .options = furrr_options(seed = TRUE)) %>%
      set_names(lca_genplots$test[[i]]$variable)

    ## additional styling: color scales

    lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_mut')] <-
      lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_mut')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$mut_status_color),
                              labels = names(globals$mut_status_color),
                              name = ''))

    lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_del')] <-
      lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_del')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$del_status_color),
                              labels = names(globals$del_status_color),
                              name = ''))

    lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_amp')] <-
      lca_genplots$plots[[i]][stri_detect(names(lca_genplots$plots[[i]]), fixed = '_amp')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$amp_status_color),
                              labels = names(globals$amp_status_color),
                              name = ''))

  }

  plan('sequential')

# END -------

  lca_genplots <-
    lca_genplots[c("lexicon",
                   "cluster_stack_plots", "cluster_onco_plots",
                   "stack_plots", "onco_plots",
                   "plots")]

  insert_tail()
