# Frequency of FGFR and FGF alterations in the molecular subtypes.
# The subtype information is available only for the TCGA and IMvigor subsets.

  insert_head()

# container ------

  sub_genet <- list()

# analysis data -------

  insert_msg('Analysis data')

  ## variables

  sub_genet$variables <- globals[c("receptors", "ligands")] %>%
    reduce(union)

  ## data, appending with the subtype information

  sub_genet$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    map(compact) %>%
    map(map, select, sample_id, any_of(sub_genet$variables)) %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map(~map2(., names(.),
              ~set_names(.x, paste(names(.x), .y, sep = '_')))) %>%
    map(map, rownames_to_column, 'sample_id') %>%
    map(function(x) if(length(x) > 1) reduce(x, left_join, by = 'sample_id') else x[[1]])

  sub_genet$data <-
    map2(map(subtypes$assignment, ~.x[c('sample_id', 'consensusClass')]),
         sub_genet$data[names(subtypes$assignment)],
         inner_join, by = 'sample_id') %>%
    map(select, -sample_id)

  ## removal of features completely absent from the data sets

  for(i in names(sub_genet$data)) {

    sub_genet$data[[i]] <- sub_genet$data[[i]] %>%
      map_dfc(function(x) if(sum(as.numeric(x)) == 0) NULL else x)

  }

  ## data set-specific variable vectors

  sub_genet$variables <- sub_genet$data %>%
    map(select, -consensusClass) %>%
    map(names)

# N numbers -------

  insert_msg('Number od observations')

  ## numbers of samples in the molecular classes

  sub_genet$n_numbers <- sub_genet$data %>%
    map(count, consensusClass)

  ## ready to use plot captions

  sub_genet$n_captions <- sub_genet$n_numbers %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ')) %>%
    map(paste, collapse = ', ')

# Descriptive stats -------

  insert_msg('Descriptive stats')

  sub_genet$stats <- sub_genet$data %>%
    map(count_binary,
        split_fct = 'consensusClass')

# Testing for differences between the cohorts ----

  insert_msg('Testing for differences between the molecular subsets')

  sub_genet$test <-
    list(x = sub_genet$data,
         y = sub_genet$variables,
         z = c(TRUE, FALSE)) %>%
    pmap(function(x, y, z) x %>%
           compare_variables(variables = y,
                             split_factor = 'consensusClass',
                             what = 'eff_size',
                             types = 'cramer_v',
                             exact = FALSE,
                             ci = FALSE,
                             pub_styled = TRUE,
                             adj_method = 'BH',
                             .parallel = z,
                             .paropts = furrr_options(seed = TRUE,
                                                      globals = 'sub_genet'))) %>%
    map(mutate,
        plot_cap = paste(eff_size, significance, sep = ', '),
        gene_symbol = stri_split_fixed(variable,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
        alteration = stri_split_fixed(variable,
                                      pattern = '_',
                                      simplify = TRUE)[, 2],
        plot_title = paste(html_italic(gene_symbol), alteration),
        hm_lab = paste(plot_title, plot_cap, sep = '<br>'),
        hm_lab = ifelse(p_adjusted < 0.05,
                        html_bold(hm_lab), hm_lab))

  ## significant effects

  sub_genet$significant <- sub_genet$test %>%
    map(filter, p_adjusted < 0.05) %>%
    map(~.x$variable)

  ## common significant effects

  sub_genet$common_significant <- sub_genet$significant %>%
    reduce(intersect)

# Stack plots for single genetic features --------

  insert_msg('Stack plots')

  plan('multisession')

  for(i in names(sub_genet$data)) {

    sub_genet$plots[[i]] <-
      list(variable = sub_genet$test[[i]]$variable,
           plot_title = sub_genet$test[[i]]$plot_title %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           plot_subtitle = sub_genet$test[[i]]$plot_cap) %>%
      future_pmap(plot_variable,
                  sub_genet$data[[i]],
                  split_factor = 'consensusClass',
                  type = 'stack',
                  scale = 'percent',
                  cust_theme = globals$common_theme +
                    theme(plot.title = element_markdown()),
                  x_n_labs = TRUE,
                  y_lab = '% of subset',
                  x_lab = 'MIBC consesus subset',
                  .options = furrr_options(seed = TRUE)) %>%
      set_names(sub_genet$test[[i]]$variable)

  }

  plan('sequential')

  ## styling

  sub_genet$plots$imvigor <-  sub_genet$plots$imvigor %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$mut_status_color),
                            labels = names(globals$mut_status_color),
                            name = ''))

  sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_mut')] <-
    sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_mut')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$mut_status_color),
                            labels = names(globals$mut_status_color),
                            name = ''))

  sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_del')] <-
    sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_del')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$del_status_color),
                            labels = names(globals$del_status_color),
                            name = ''))

  sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_amp')] <-
    sub_genet$plots$tcga[stri_detect(names(sub_genet$plots$tcga), fixed = '_amp')] %>%
    map(~.x +
          scale_fill_manual(values = unname(globals$amp_status_color),
                            labels = names(globals$amp_status_color),
                            name = ''))

# Result table --------

  insert_msg('Result table')

  sub_genet$result_tbl <-
    map2(sub_genet$stats,
         sub_genet$test,
         merge_stat_test) %>%
    compress(names_to = 'cohort') %>%
    mutate(gene_symbol = stri_split_fixed(variable,
                                          pattern = '_',
                                          simplify = TRUE)[, 1],
           alteration = stri_split_fixed(variable,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           cohort = globals$cohort_labs[cohort],
           percent = signif(percent, 2)) %>%
    arrange(cohort, variable, consensusClass) %>%
    select(cohort,
           gene_symbol,
           alteration,
           consensusClass,
           percent,
           n,
           n_total,
           significance,
           eff_size) %>%
    set_names(c('Cohort', 'Gene symbol', 'Alteration',
                'MIBC consensus class',
                'Percentage of samples with alterations',
                'Samples with alterations, N',
                'Total samples, N',
                'Significance',
                'Effects size'))

# Heat map representation of alterations of interest -------

  insert_msg('Heat map')

  ## these are: FGFR3 mutations and amplifications of FGF3, FGF4, and FGF19
  ## the full representation is possible only for the TCGA cohort

  sub_genet$hm_variables <-
    c('FGFR1_mutation',
      'FGFR2_mutation',
      'FGFR3_mutation',
      'FGFR4_mutation',
      'FGFR1_amplification',
      'FGFR2_amplification',
      'FGFR3_amplification',
      'FGFR4_amplification',
      'FGF3_amplification',
      'FGF4_amplification',
      'FGF19_amplification')

  sub_genet$hm_data <- sub_genet$data %>%
    map(select, consensusClass, any_of(sub_genet$hm_variables)) %>%
    map(~mutate(.x, observation = paste0('sample_', 1:nrow(.x)))) %>%
    map(pivot_longer,
        cols = any_of(sub_genet$hm_variables),
        names_to = 'variable',
        values_to = 'alteration') %>%
    map(mutate,
        variable = factor(variable, rev(sub_genet$hm_variables)),
        gene_symbol = stri_split_fixed(variable,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
        gene_group = ifelse(gene_symbol %in% globals$receptors,
                            'receptor', 'ligand'),
        gene_group = factor(gene_group, c('receptor', 'ligand')))

  ## heat maps

  sub_genet$hm_plots <-
    list(x = sub_genet$hm_data,
         y = globals$cohort_labs[names(sub_genet$hm_data)],
         z = sub_genet$n_captions,
         w = sub_genet$test) %>%
    pmap(function(x, y, z, w) x %>%
           ggplot(aes(x = reorder(observation, -alteration),
                      y = variable,
                      fill = factor(alteration))) +
           geom_tile() +
           facet_grid(gene_group ~ consensusClass,
                      scales = 'free',
                      space = 'free') +
           scale_y_discrete(labels = set_names(w$hm_lab, w$variable)) +
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
                subtitle = z,
                x = 'sample'))

# END ------

  sub_genet$data <- NULL
  sub_genet$hm_data <- NULL

  sub_genet <- compact(sub_genet)

  rm(i)

  insert_tail()
