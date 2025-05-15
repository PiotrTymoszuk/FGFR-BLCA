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
      mutate(cohort = factor(cohort, globals$analysis_cohorts),
             cohort = droplevels(cohort)) %>%
      relocate(cohort)

  }

  ## variable lexicon with ready-to-use axis labels in HTML format
  ## components of the FGF/FGFR system are highlighted with bold font

  expl_genet$fgf_genes <-
    globals[c("receptors", "ligands", "binding_proteins")]

  expl_genet$lexicon <-
    tibble(variable = reduce(expl_genet$variables, union)) %>%
    mutate(label = html_italic(variable),
           label = ifelse(variable %in% reduce(expl_genet$fgf_genes, union),
                          html_bold(label), label))

  ## data for the FGF, FGFR, and FGFBP-coding genes
  ## detected in at least two cohorts:
  ## if a gene is not present in a cohort it is assumed all WT

  expl_genet$fgf_data <- globals$cohort_expr %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    transpose %>%
    map(compact) %>%
    map(map,
        select,
        any_of(reduce(expl_genet$fgf_genes, union)))

  expl_genet$fgf_genes <- expl_genet$fgf_data %>%
    map(map, colSums) %>%
    map(map, ~.x[.x > 0]) %>%
    map(map, names) %>%
    map(shared_features, m = 2) %>%
    map(as.character)

  expl_genet$fgf_data <-
    map2(expl_genet$fgf_data,
         expl_genet$fgf_genes,
         function(x, y) x %>%
           map(select, any_of(y))) %>%
    map(~map2_dfr(.x, names(.x),
                  function(dat, coh) mutate(dat, cohort = coh))) %>%
    map(relocate, cohort)

  expl_genet$fgf_data <- expl_genet$fgf_data %>%
    map(map_dfc, ~ifelse(is.na(.x), 0, .x))

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

  ## frequency of all alterations,
  ## and alterations in the FGF, FGFR, and FGFBP genes

  expl_genet[c("stats", "fgf_stats")] <-
    expl_genet[c("data", "fgf_data")] %>%
    map(map,
        count_binary,
        split_fct = 'cohort') %>%
    map(map,
        mutate,
        cohort = factor(cohort, globals$analysis_cohorts),
        cohort = droplevels(cohort))

  expl_genet$fgf_stats <- expl_genet$fgf_stats %>%
    map(mutate,
        gene_group = ifelse(variable %in% globals$receptors,
                            'receptor',
                            ifelse(variable %in% globals$ligands,
                                   'ligand', 'BP')),
        gene_group = factor(gene_group,
                            c('receptor', 'ligand', 'BP')))

  ## top alterations present in at least 5% of samples
  ## in all cohorts

  expl_genet$top_alterations <- expl_genet$stats %>%
    map(blast, cohort) %>%
    map(map, filter, percent >= 5) %>%
    map(map, ~.x$variable) %>%
    map(reduce, intersect)

  ## top FGF/FGFBP/FGFR alterations present in at least 1% of samples
  ## in all cohorts

  expl_genet$fgf_top_alterations <- expl_genet$fgf_stats %>%
    map(blast, cohort) %>%
    map(map, filter, percent >= 1) %>%
    map(map, ~.x$variable) %>%
    map(reduce, intersect)

# Testing for differences between the cohorts ------

  insert_msg('Testing for differences between the cohorts')

  ## done only for the top most frequent alterations

  expl_genet$test <-
    map2(expl_genet$data,
         expl_genet$top_alterations,
         ~f_chisq_test(.x[.y],
                       f = .x$cohort,
                       safely = TRUE,
                       as_data_frame = TRUE,
                       adj_method = 'BH')) %>%
    map(re_adjust) %>%
    map(p_formatter, text = TRUE) %>%
    map(mutate,
        eff_size = paste('V =', signif(cramer_v, 2)),
        plot_cap = paste(eff_size, significance, sep = ', ')) %>%
    map(as_tibble)

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
                    color = 'white',
                    position = position_dodge(0.9)) +
           geom_vline(xintercept = 0,
                      linetype = 'solid',
                      linewidth = 0.15) +
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
         expl_genet$test %>%
           map(p_formatter, text = FALSE) %>%
           map(~.x[c('variable', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map2(., expl_genet$top_alterations,
         ~filter(.x, variable %in% .y)) %>%
    map(format_genetics)

# Detailed bar plots for the genes of interest -------

  insert_msg('Detailed bar plots for the genes of interest')

  ## bar plots

  expl_genet$fgf_plots <-
    list(x = expl_genet$fgf_stats,
         y = paste0(c('Somatic mutations',
                      'Gene deletion',
                      'Gene amplification'),
                    ', <em>FGFR/FGF/FGFBP</em> genes'),
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
                    color = 'white',
                    position = position_dodge(0.9),
                    linewidth = 0.01) +
           geom_vline(xintercept = 0,
                      linetype = 'solid',
                      linewidth = 0.15) +
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
    map(format_genetics) %>%
    map(mutate,
        `Gene group` = car::recode(as.character(gene_group),
                                 "'BP' = 'binding protein'")) %>%
    map(relocate, `Gene group`) %>%
    map(select, -gene_group)

# END -------

  expl_genet$data <- NULL
  expl_genet$fgf_stats <- NULL
  expl_genet$fgf_data <- NULL

  expl_genet <- compact(expl_genet)

  insert_tail()
