# Quality check for the genetic clusters

  insert_head()

# container ------

  lca_eval <- list()

# cluster objects and assignment frames --------

  insert_msg('Cluster objects and assignment frames')

  lca_eval$clust_objects <- lca_pred$genet_objects

  lca_eval$assignment <- lca_pred$assignment

  ## tibbles with factor-coded genetic features

  lca_eval$variables <- lca_pred$variables

  lca_eval$data <-
    map2(lca_eval$assignment,
         lca_globals$data,
         inner_join,
         by = 'sample_id')

# N numbers and cluster sizes ------

  insert_msg('Cluster sizes')

  ## total numbers and numbers of samples in the clusters

  lca_eval$n_total_labs <- lca_eval$data %>%
    map_dbl(nrow) %>%
    map2_chr(names(.), .,
             ~paste(globals$cohort_labs[.x], .y, sep = '\nn = '))

  lca_eval$n_numbers <- lca_eval$clust_objects %>%
    map(ngroups) %>%
    map(mutate,
        n_total = sum(n),
        percent = n/n_total * 100) %>%
    map(arrange, desc(clust_id)) %>%
    map(mutate,
        plot_pos = cumsum(percent) - 0.5 * percent)

  lca_eval$n_clust_labs <- lca_eval$n_numbers %>%
    map(arrange, clust_id) %>%
    map(~map2_chr(.x[[1]], .x[[2]],
                  paste, sep = ': n = ')) %>%
    map2(., lca_eval$n_numbers,
         ~set_names(.x, .y[[1]]))

  ## stack plots

  lca_eval$n_plot <- lca_eval$n_numbers %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, globals$analysis_cohorts),
           cohort = droplevels(cohort)) %>%
    ggplot(aes(x = percent,
               y = cohort,
               fill = clust_id)) +
    geom_bar(stat = 'identity',
             color = 'black') +
    geom_label(aes(label = signif(percent, 2),
                   x = plot_pos),
               size = 2.5,
               show.legend = FALSE) +
    scale_fill_manual(values = globals$genet_colors,
                      name = 'Genetic cluster') +
    scale_y_discrete(labels = lca_eval$n_total_labs) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Genetic cluster distribution',
         x = '% of samples')

# Numeric stats of cluster performance -------

  insert_msg('Numeric stats of cluster performance')

  ## silhouette widths do not make much sense
  ## in the evidently non-Euclidean setting

  lca_eval$clust_stats <- lca_eval$clust_objects %>%
    map(summary) %>%
    compress(names_to = 'cohort')

# Clustering variables: frequency in the clusters --------

  insert_msg('Clustering features, frequency')

  lca_eval$stats <- lca_eval$data %>%
    map(select, -sample_id) %>%
    map(bin2stats,
        split_factor = 'clust_id') %>%
    map(mutate,
        gene_symbol = stri_split_fixed(variable,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
        alteration = stri_split_fixed(variable,
                                      pattern = '_',
                                      simplify = TRUE)[, 2])

# Chi-square tests -------

  insert_msg('Chi-square tests')

  plan('multisession')

  lca_eval$test <- lca_eval$data %>%
    future_map(compare_variables,
               variables = lca_eval$variables,
               split_factor = 'clust_id',
               what = 'eff_size',
               types = 'cramer_v',
               exact = FALSE,
               ci = FALSE,
               adj_method = 'BH',
               pub_styled = TRUE,
               .parallel = TRUE,
               .options = furrr_options(seed = TRUE)) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '),
        gene_symbol = stri_split_fixed(variable,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
        alteration = stri_split_fixed(variable,
                                      pattern = '_',
                                      simplify = TRUE)[, 2],
        plot_title = paste(html_italic(gene_symbol), alteration),
        ax_lab = paste(plot_title, plot_cap, sep = '<br>'),
        ax_lab = ifelse(p_adjusted < 0.05,
                        html_bold(ax_lab), ax_lab),
        gene_html = ifelse(p_adjusted < 0.05,
                           html_bold(html_italic(gene_symbol)),
                           html_italic(gene_symbol)))

  plan('sequential')

  ## significant differences

  lca_eval$significant <- lca_eval$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

  lca_eval$common_significant <- lca_eval$significant %>%
    reduce(intersect)

# Stack plots for single genetic features ---------

  insert_msg('Stack plots')

  for(i in names(lca_eval$data)) {

    lca_eval$plots[[i]] <-
      list(variable = lca_eval$test[[i]]$variable,
           plot_title = lca_eval$test[[i]]$plot_title %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           plot_subtitle = lca_eval$test[[i]]$plot_cap) %>%
      pmap(plot_variable,
           lca_eval$data[[i]],
           split_factor = 'clust_id',
           type = 'stack',
           scale = 'percent',
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown()),
           x_n_labs = TRUE,
           x_lab = 'Genetic cluster',
           y_lab = '% of cluster') %>%
      set_names(lca_eval$test[[i]]$variable)

    ## additional styling, setting the fill scales appropriate
    ## for the alteration type

    lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_mut')] <-
      lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_mut')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$mut_status_color),
                              labels = names(globals$mut_status_color),
                              name = ''))

    lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_del')] <-
      lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_del')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$del_status_color),
                              labels = names(globals$del_status_color),
                              name = ''))

    lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_amp')] <-
      lca_eval$plots[[i]][stri_detect(names(lca_eval$plots[[i]]), fixed = '_amp')] %>%
      map(~.x +
            scale_fill_manual(values = unname(globals$amp_status_color),
                              labels = names(globals$amp_status_color),
                              name = ''))

  }

# Heat maps of the clustering features -------

  insert_msg('Classification and ')

  ## classification of the features by their character and  general frequency
  ## in the genie cohort

  lca_eval$hm_classification <- lca_eval$stats$genie %>%
    group_by(variable, alteration) %>%
    summarise(n = sum(n)) %>%
    ungroup %>%
    mutate(alteration = factor(alteration,
                               c('mutation', 'deletion', 'amplification'))) %>%
    arrange(alteration, n)

  ## plotting data

  lca_eval$hm_data <- lca_eval$data %>%
    map(select, clust_id, all_of(lca_eval$variables)) %>%
    map(~mutate(.x, observation = paste0('sample_', 1:nrow(.x)))) %>%
    map(pivot_longer,
        cols = any_of(lca_eval$variables),
        names_to = 'variable',
        values_to = 'feature') %>%
    map(left_join,
        lca_eval$hm_classification[c('variable', 'alteration')],
        by = 'variable') %>%
    map(mutate,
        variable = factor(variable, lca_eval$hm_classification$variable))

  ## heat maps

  lca_eval$hm_plots <-
    list(x = lca_eval$hm_data,
         y = globals$cohort_labs[names(lca_eval$hm_data)],
         z = lca_eval$n_clust_labs,
         w = lca_eval$test) %>%
    pmap(function(x, y, z, w) x %>%
           ggplot(aes(x = reorder(observation, -feature),
                      y = variable,
                      fill = factor(feature))) +
           geom_tile() +
           facet_grid(alteration ~ clust_id,
                      scales = 'free',
                      space = 'free') +
           scale_y_discrete(labels = set_names(w$gene_html, w$variable)) +
           scale_fill_manual(values = c('0' = 'gray70',
                                        '1' = 'orangered3'),
                             labels = c('0' = 'absent',
                                        '1' = 'present'),
                             name = '') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_markdown(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.line = element_blank(),
                 panel.grid.major.x = element_blank(),
                 strip.text.x = element_text(hjust = 0, angle = 90)) +
           labs(title = y,
                subtitle = paste(z, collapse = ', '),
                x = 'sample'))

# Result table -------

  insert_msg('Result table')

  lca_eval$result_tbl <-
    map2(lca_eval$stats,
         lca_eval$test,
         merge_stat_test) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = globals$cohort_labs[cohort],
           percent = signif(percent, 3),
           variable = stri_replace(variable,
                                   fixed = '_',
                                   replacement = ' ')) %>%
    select(cohort, variable, clust_id,
           percent, n, n_total,
           significance, eff_size) %>%
    set_names(c('Cohort', 'Variable', 'Genetic subset',
                'Samples with alterations, percentage',
                'Samples with alterations, N',
                'Tota samples, N',
                'Significance', 'Effect size'))

# END ------

  rm(i)

  lca_eval$data <- NULL
  lca_eval$hm_classification <- NULL
  lca_eval$hm_data <- NULL

  lca_eval <- compact(lca_eval)

  insert_tail()
