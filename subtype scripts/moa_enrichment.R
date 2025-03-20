# Enrichment analysis of mechanisms of actions (MOA) and molecular targets
# for drugs that differed in predicted resistance between MIBC consensus
# classes in at least two cohorts.
#
# The enrichment analysis is done separately for the common sensitive and
# resistant effects and employs random sampling test.
# significant enrichment is considered for OR >= 1.44 and pFDR < 0.05.


  insert_head()

# container --------

  sub_moa <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## the drugs of interest
  ## dictionaries of drugs assigned to MOA and molecular targets
  ## removal of entries with single drugs

  sub_moa$drugs <- sub_drugs$common_significant %>%
    unlist(recursive = FALSE)

  sub_moa$dictionary[c('moa', 'targets')] <-
    drugs[c("moa_dictionary", "target_dictinoary")]

  sub_moa$dictionary$targets <- sub_moa$dictionary$targets %>%
    map(function(x) if(length(x) == 1) NULL else x) %>%
    compact

  ## vectors with all drugs present in the dictionaries

  sub_moa$all <- sub_moa$dictionary %>%
    map(reduce, union)

  ## levels of the consensus classes

  sub_moa$class_levels <- sub_drugs$common_significant %>%
    names

  ## plot titles, suffixes, and fill colors

  sub_moa$plot_titles <-
    map2_chr(sub_drugs$common_significant %>%
               names %>%
               map(rep, 2) %>%
               unlist,
             sub_drugs$common_significant[[1]] %>%
               names %>%
               rep(length(sub_drugs$common_significant)),
             ~paste(.x, 'vs LumP,', .y)) %>%
    set_names(names(sub_moa$drugs))

  sub_moa$plot_suffixes <-
    c(moa = 'MOA',
      targets = 'molecular targets')

  sub_moa$plot_colors <-
    c(moa = 'aquamarine4',
      targets = 'lightskyblue3')

# The enrichment analysis ---------

  insert_msg('Enrichment analysis')

  set.seed(12345)

  for(i in names(sub_moa$dictionary)) {

    sub_moa$test[[i]] <- sub_moa$drugs %>%
      map(~f_enrichment(.x[.x %in% sub_moa$all[[i]]],
                        dict = sub_moa$dictionary[[i]],
                        all = sub_moa$all[[i]],
                        type = 'fisher',
                        adj_method = 'BH')) %>%
      map(as.data.frame) %>%
      map(rownames_to_column, 'variable') %>%
      map(as_tibble)

  }

  ## formatting the enrichment testing results

  sub_moa$test <- sub_moa$test %>%
    map(map, re_adjust) %>%
    map(map,
        mutate,
        log_or = log(or),
        regulation = ifelse(p_adjusted >= 0.05, 'ns',
                            ifelse(log_or > log(1.44), 'enriched',
                                   ifelse(log_or <= -log(1.44),
                                          'depleted', 'ns'))),
        regulation = factor(regulation, c('enriched', 'depleted', 'ns')))

# significant enrichment ---------

  insert_msg('Significant enrichment')

  sub_moa$significant <- sub_moa$test %>%
    map(map, filter, regulation == 'enriched') %>%
    map(map, ~.x$variable)

# volcano plots -------

  insert_msg('volcano plots')

  for(i in names(sub_moa$test)) {

    sub_moa$volcano_plots[[i]] <-
      list(data = sub_moa$test[[i]],
           plot_title = paste(sub_moa$plot_titles,
                              sub_moa$plot_suffixes,
                              sep = ', ')) %>%
      pmap(plot_volcano,
           regulation_variable = 'log_or',
           p_variable = 'p_adjusted',
           regulation_level = log(1.44),
           x_lab = 'enrichment over all drugs, log OR',
           y_lab = expression('-log'[10] * ' pFDR'),
           cust_theme = globals$common_theme +
             theme(plot.title = element_markdown(),
                   legend.title = element_markdown()),
           label_variable = 'variable',
           top_significant = 20,
           txt_size = 2.5,
           txt_face = 'plain',
           label_type = 'text') %>%
      map(tag2subtitle) %>%
      map(~.x +
            scale_fill_manual(values = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray60'),
                              labels = c(upregulated = 'enriched',
                                         downregulated = 'depleted',
                                         ns = 'ns'),
                              name = paste('enrichment vs all drugs')))


  }

# Dot plots ----------

  insert_msg('Dot plots')

  ## plotting data

  sub_moa$dot_data <- sub_moa$test %>%
    map(map, filter, regulation == 'enriched') %>%
    map(compress, names_to = 'resistance') %>%
    map(mutate,
        consensusClass = stri_split_fixed(resistance,
                                          pattern = '.',
                                          simplify = TRUE)[, 1],
        resistance = stri_split_fixed(resistance,
                                      pattern = '.',
                                      simplify = TRUE)[, 2],
        consensusClass = factor(consensusClass, sub_moa$class_levels),
        resistance = factor(resistance,
                            c('resistant', 'sensitive')),
        p_adjusted = ifelse(p_adjusted == 0,
                            min(p_adjusted[p_adjusted != 0]),
                            p_adjusted)) %>%
    map(select, consensusClass, variable, resistance, p_adjusted, or) %>%
    map(blast, consensusClass)

  ## the plots: OR in the X axis, point size codes for the FDR-corrected
  ## p values

  for(i in names(sub_moa$dot_data)) {

    sub_moa$dot_plots[[i]] <-
      list(x = sub_moa$dot_data[[i]],
           y = sub_moa$dot_data[[i]] %>%
             names %>%
             paste('vs LumP,',
                   sub_moa$plot_suffixes[[i]])) %>%
      pmap(function(x, y) x %>%
             ggplot(aes(x = or,
                        y = reorder(variable, or),
                        size = -log10(p_adjusted),
                        fill = resistance)) +
             facet_grid(resistance ~ .,
                        scales = 'free',
                        space = 'free') +
             geom_point(shape = 21) +
             scale_size_area(max_size = 4.5,
                             labels = function(x) signif(10^-x, 2),
                             name = 'pFDR') +
             scale_y_discrete(labels = function(x) stri_replace(x,
                                                                fixed = ' ',
                                                                replacement = '\n')) +
             scale_fill_manual(values = c(resistant = 'firebrick',
                                          sensitive = 'steelblue'),
                               name = 'resistance vs LumP') +
             globals$common_theme +
             theme(axis.title.y = element_blank(),
                   plot.title = element_markdown()) +
             labs(title = y,
                  x = 'enrichment over all drugs, OR'))

  }

# END ------

  sub_moa$drugs <- NULL
  sub_moa$dictionary <- NULL
  sub_moa$all <- NULL
  sub_moa$dot_data <- NULL

  sub_moa <- compact(sub_moa)

  insert_tail()
