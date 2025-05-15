# Enrichment analysis of mechanisms of actions (MOA) and molecular targets
# for drugs that differed in predicted resistance between FGFR3 WT and FGFRF3
# mutated cancers in at least two cohorts.
#
# The enrichment analysis is done separately for the common sensitive and
# resistant effects and employs random sampling test.
# significant enrichment is considered for OR >= 1.44 and pFDR < 0.05.


  insert_head()

# container --------

  fgfr_moa <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## the drugs of interest
  ## dictionaries of drugs assigned to MOA and molecular targets
  ## removal of entries with single drugs

  fgfr_moa$drugs <- fgfr_drugs$common_significant

  fgfr_moa$dictionary[c('moa', 'targets')] <-
    drugs[c("moa_dictionary", "target_dictinoary")]

  fgfr_moa$dictionary$targets <- fgfr_moa$dictionary$targets %>%
    map(function(x) if(length(x) == 1) NULL else x) %>%
    compact

  ## vectors with all drugs present in the dictionaries

  fgfr_moa$all <- fgfr_moa$dictionary %>%
    map(reduce, union)

# The enrichment analysis ---------

  insert_msg('Enrichment analysis')

  set.seed(12345)

  for(i in names(fgfr_moa$dictionary)) {

    fgfr_moa$test[[i]] <- fgfr_moa$drugs %>%
      map(~f_enrichment(.x[.x %in% fgfr_moa$all[[i]]],
                        dict = fgfr_moa$dictionary[[i]],
                        all = fgfr_moa$all[[i]],
                        type = 'fisher',
                        adj_method = 'BH')) %>%
      map(as.data.frame) %>%
      map(rownames_to_column, 'variable') %>%
      map(as_tibble)

  }

  ## formatting the enrichment testing results

  fgfr_moa$test <- fgfr_moa$test %>%
    map(map, re_adjust) %>%
    map(map, p_formatter, text = TRUE) %>%
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

  fgfr_moa$significant <- fgfr_moa$test %>%
    map(map, filter, regulation == 'enriched') %>%
    map(map, ~.x$variable)

# volcano plots -------

  insert_msg('volcano plots')

  for(i in names(fgfr_moa$test)) {

    fgfr_moa$volcano_plots[[i]] <-
      list(data = fgfr_moa$test[[i]],
           plot_title = fgfr_moa$volcano_plots[[i]] %>%
             names %>%
             stri_capitalize_first %>%
             paste('in', html_italic('FGFR3'), 'mutated vs WT')) %>%
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

  fgfr_moa$dot_data <- fgfr_moa$test %>%
    map(map, filter, regulation == 'enriched') %>%
    map(compress, names_to = 'resistance') %>%
    map(mutate,
        resistance = factor(resistance, c('resistant', 'sensitive')),
        p_adjusted = ifelse(p_adjusted == 0,
                            min(p_adjusted[p_adjusted != 0]),
                            p_adjusted)) %>%
    map(select, variable, resistance, p_adjusted, or)

  ## the plots: OR in the X axis, point size codes for the FDR-corrected
  ## p values

  fgfr_moa$dot_plots <-
    list(x = fgfr_moa$dot_data,
         y = paste('Predicted drug response',
                   html_italic('FGFR3'),
                   'mutated vs WT,',
                   c('MOA', 'molecular targets'))) %>%
    pmap(function(x, y, z) x %>%
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
                             name = 'resistance vs WT')  +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 plot.title = element_markdown()) +
           labs(title = y,
                x = 'enrichment over all drugs, OR'))

# END ------

  fgfr_moa$drugs <- NULL
  fgfr_moa$dictionary <- NULL
  fgfr_moa$all <- NULL
  fgfr_moa$dot_data <- NULL

  fgfr_moa <- compact(fgfr_moa)

  insert_tail()
