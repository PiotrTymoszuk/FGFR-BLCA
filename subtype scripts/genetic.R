# Frequency of FGFR and FGF alterations in the molecular subtypes.
# The subtype information is available only for the TCGA, IMvigor,
# and BCAN cohorts.
#
# Statistical significance is terted by weighted permutations implemented
# by perich.

  insert_head()

# container ------

  sub_genet <- list()

# parallel backend --------

  insert_msg('Parallel backend')

  plan('multisession')

# analysis data -------

  insert_msg('Analysis data')

  ## variables

  sub_genet$variables <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    reduce(union)

  ## data, appending with the subtype information

  sub_genet$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x[c('mutation', 'deletion', 'amplification')]) %>%
    map(compact) %>%
    map(map, column_to_rownames, 'sample_id') %>%
    map(map, ~.x[!names(.x) %in% c('sample_id', 'patient_id', 'tissue_type')]) %>%
    map(~map2(., names(.),
              ~set_names(.x, paste(names(.x), .y, sep = '_')))) %>%
    map(map, rownames_to_column, 'sample_id') %>%
    map(function(x) if(length(x) > 1) reduce(x, left_join, by = 'sample_id') else x[[1]])

  sub_genet$data <-
    map2(map(subtypes$assignment, ~.x[c('sample_id', 'consensusClass')]),
         sub_genet$data[names(subtypes$assignment)],
         inner_join, by = 'sample_id') %>%
    map(select, -sample_id) %>%
    map(filter, !is.na(consensusClass))

  ## removal of features completely absent from the data sets

  for(i in names(sub_genet$data)) {

    sub_genet$data[[i]] <- sub_genet$data[[i]] %>%
      map_dfc(function(x) if(sum(as.numeric(x), na.rm = TRUE) == 0) NULL else x)

  }

  ## removal of the LumNS, NE-like, and not assigned subsets,
  ## not enought cases for an analysis in the IMvigor and BCAN cohorts

  sub_genet$data <- sub_genet$data %>%
    map(filter,
        !consensusClass %in% c('not assigned', 'LumNS', 'NE-like')) %>%
    map(mutate,
        consensusClass = droplevels(consensusClass))

  ## data set-specific variable vectors

  sub_genet$variables <- sub_genet$data %>%
    map(select, -consensusClass) %>%
    map(names) %>%
    map(~.x[stri_detect(.x,
                        regex = sub_genet$variables %>%
                          paste(collapse = '|') %>%
                          paste0('^(', ., ')_'))])

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

  sub_genet$stats <-
    list(data = sub_genet$data,
         variables = sub_genet$variables) %>%
    pmap(count_binary,
         split_fct = 'consensusClass')

# Testing for differences between the cohorts ----

  insert_msg('Testing for differences between the molecular subsets')

  ## done with permutation testing

  sub_genet$test <-
    list(x = sub_genet$data,
         y = sub_genet$variables) %>%
    pmap(function(x, y) perTest(x = x[y],
                                f = x[['consensusClass']],
                                background = x[, -1],
                                alternative = 'greater',
                                n_iter = 40000,
                                as_data_frame = TRUE,
                                compress = TRUE,
                                adj_method = 'BH'))

  ## formatting the results:
  ## plot captions with the global effect size and the global p value
  ## identification of significantly enriched alterations

  sub_genet$test <- sub_genet$test %>%
    map(re_adjust, 'global_p_value') %>%
    map(p_formatter, text = TRUE) %>%
    map(as_tibble) %>%
    map(group_by, variable) %>%
    map(mutate,
        cramer_v = sqrt(chisq/sum(n_total)/df),
        global_significance = significance,
        global_eff_size = paste('V =', signif(cramer_v, 2)),
        plot_cap = paste(global_eff_size, global_significance, sep = ', ')) %>%
    map(ungroup)

  sub_genet$test <- sub_genet$test %>%
    map(re_adjust) %>%
    map(p_formatter, text = TRUE) %>%
    map(mutate,
        cluster_significance = significance,
        regulation = ifelse(p_adjusted < 0.05 & ES >= 1.2,
                            'enriched', 'ns'),
        regulation = factor(regulation, c('enriched', 'ns')))

  ## alteration type

  sub_genet$test <- sub_genet$test %>%
    map(mutate,
        gene_symbol = stri_split_fixed(variable,
                                       pattern = '_',
                                       simplify = TRUE)[, 1],
        alteration = stri_split_fixed(variable,
                                      pattern = '_',
                                      simplify = TRUE)[, 2])

  ## significant effects

  sub_genet$significant <- sub_genet$test %>%
    map(filter, regulation == 'enriched') %>%
    map(~.x$variable) %>%
    map(unique)

  ## common significant effects

  sub_genet$common_significant <- sub_genet$significant %>%
    shared_features(m = 2) %>%
    as.character

  ## global testing results to be displayed in the plots

  sub_genet$stack_test <- sub_genet$test %>%
    map(filter, !duplicated(variable)) %>%
    map(mutate,
        alteration_short = car::recode(alteration,
                                       "'mutation' = 'mut';
                                       'deletion' = 'del';
                                       'amplification' = 'amp'")) %>%
    map(arrange, p_value)

# Stack plots for single genetic features --------

  insert_msg('Stack plots')

  ## plots: top 5 significant effects

  for(i in names(sub_genet$data)) {

    sub_genet$plots[[i]] <-
      list(variable = sub_genet$stack_test[[i]][1:5, ]$variable,
           plot_title = sub_genet$stack_test[[i]][1:5, ]$gene_symbol %>%
             html_italic %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           plot_subtitle = sub_genet$stack_test[[i]][1:5, ]$plot_cap) %>%
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
      set_names(sub_genet$stack_test[[i]][1:5, ]$variable)

  }

  ## styling

  sub_genet$plots[c("imvigor", "bcan")] <-
    sub_genet$plots[c("imvigor", "bcan")] %>%
    map(map,
        ~.x +
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
         sub_genet$stack_test %>%
           map(mutate,
               eff_size = global_eff_size,
               significance = global_significance),
         merge_stat_test) %>%
    compress(names_to = 'cohort') %>%
    format_genetics %>%
    mutate(Alteration = stri_split_fixed(`Gene symbol`,
                                         pattern = '_',
                                         simplify = TRUE)[, 2],
           `Gene symbol` = stri_split_fixed(`Gene symbol`,
                                            pattern = '_',
                                            simplify = TRUE)[, 1],
           n_total = stri_extract(`Cohort (n samples)`,
                                  regex = "\\d+"),
           `MIBC consensus class (n samples)` = paste0(consensusClass,
                                                       " (", n_total, ")"),
           `Cohort` = stri_replace(`Cohort (n samples)`,
                                   regex = "\\s{1}\\(\\d+\\)",
                                   replacement = "")) %>%
    arrange(`Cohort`, `Gene symbol`, consensusClass) %>%
    select(`Cohort`,
           `Gene symbol`,
           Alteration,
           `MIBC consensus class (n samples)`,
           `Alteration frequency, n samples (percentage)`,
           `FDR p value`,
           `Effect size, Cramer's V`)

# Oncoplot representation of alterations of interest -------

  insert_msg('Oncoplots')

  ## these are: FGFR3 mutations and amplifications of FGF3, FGF4, and FGF19
  ## the full representation is possible only for the TCGA cohort

  sub_genet$onco_variables <-
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

  sub_genet$onco_data <- sub_genet$data %>%
    map(select, consensusClass, any_of(sub_genet$onco_variables)) %>%
    map(mutate,
        consensusClass = fct_recode(consensusClass,
                                    `Stroma\nrich` = 'Stroma-rich'))

  sub_genet$onco_variables <- sub_genet$onco_data %>%
    map(names) %>%
    map(intersect, sub_genet$onco_variables)

  ## Y axis labels: significant effects are highlighted

  sub_genet$axis_labs <-
    map2(sub_genet$stack_test,
         sub_genet$significant,
         ~mutate(.x,
                 axis_lab = paste(html_italic(gene_symbol),
                                  alteration_short),
                 axis_lab = ifelse(variable %in% .y,
                                   html_bold(axis_lab), axis_lab))) %>%
    map(~set_names(.x$axis_lab, .x$variable))

  ## heat maps

  sub_genet$onco_plots <-
    list(data = sub_genet$onco_data,
         variables = sub_genet$onco_variables,
         plot_title = globals$cohort_labs[names(sub_genet$onco_data)]) %>%
    pmap(plot_bionco,
         split_fct = 'consensusClass',
         hide_x_axis_text = TRUE,
         cust_theme = globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.line = element_blank(),
                 panel.grid.major.x = element_blank(),
                 strip.text.x = element_text(hjust = 0, angle = 90)),
         color_scale = c('gray85', 'orangered3'),
         one_plot = FALSE,
         x_lab = 'cancer sample') %>%
    map(~.x$main) %>%
    map2(., sub_genet$axis_labs,
         ~.x +
           labs(subtitle = stri_replace(.x$labels$subtitle,
                                        fixed = '\n',
                                        replacement = '-')) +
           theme(axis.text.y = element_markdown()) +
           scale_y_discrete(labels = .y))

# Caching --------

  insert_msg('Caching')

  #sub_genet$data <- NULL
  #sub_genet$hm_data <- NULL
  #sub_genet$onco_data <- NULL
  #sub_genet$stack_test <- NULL

  sub_genet <- compact(sub_genet)

  save(sub_genet, file = './cache/sub_genet.RData')

# END ------

  rm(i)

  plan('sequential')

  insert_tail()
