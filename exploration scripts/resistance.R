# Exploratory analysis for resistance to pan-FGFRi in DepMap cell lines

  insert_tail()

# container -------

  expl_res <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## resistance data and FGFR3 mutation status
  ## doxorubicin serves as a positive sensitivity control (one per data set!)

  expl_res$data <- depmap[c("gdsc", "prism")]

  expl_res$data$gdsc <- expl_res$data$gdsc %>%
    select(-GDSC1_133_Doxorubicin, -GDSC1_1421_FGFR_0939)

  expl_res$variables <- expl_res$data %>%
    map(select, -sample_id) %>%
    map(names)

  expl_res$lexicon <- expl_res$variables %>%
    map(~filter(depmap$lexicon,
                variable %in% .x))

  ## GDSC: conversion to identity scale

  expl_res$data$gdsc[expl_res$variables$gdsc] <-
    expl_res$data$gdsc[expl_res$variables$gdsc] %>%
    map_dfc(~2^.x)

  ## FGFR3 mutation status

  expl_res$FGFR3_mutation <-
    depmap$mutation[, c("sample_id", "FGFR3_mut")] %>%
    transmute(sample_id = sample_id,
              FGFR3_mutation = fct_recode(FGFR3_mut,
                                       WT = 'no',
                                       mutated = 'yes'))

  expl_res$data <- expl_res$data %>%
    map(left_join,
        expl_res$FGFR3_mutation,
        by = 'sample_id') %>%
    map(relocate, sample_id, FGFR3_mutation)

  ## identification of the sensitive cell lines:
  ## for GDSC - by IC50 below the maximal concentration in the assay
  ## for PRISM - by AUC < 1

  expl_res$sensitive$gdsc <- expl_res$data$gdsc

  expl_res$sensitive$gdsc[expl_res$lexicon$gdsc$variable] <-
    map2_dfc(expl_res$sensitive$gdsc[expl_res$lexicon$gdsc$variable],
             expl_res$lexicon$gdsc$max_conc,
             ~cut(.x,
                  c(-Inf, .y, Inf),
                  c('sensitive', 'resistant'),
                  right = FALSE))

  expl_res$sensitive$prism <- expl_res$data$prism

  expl_res$sensitive$prism[expl_res$variables$prism] <-
    expl_res$sensitive$prism[expl_res$variables$prism] %>%
    map_dfc(cut,
            c(-Inf, 1, Inf),
            c('sensitive', 'resistant'),
            right = FALSE)

  ## edits of the variable lexicon

  expl_res$lexicon <- expl_res$lexicon %>%
    reduce(rbind) %>%
    mutate(data_set = stri_extract(variable, regex = 'GDSC1|GDSC2|PRISM'))

  ## long-format data to be used for plots, appending with the cell line
  ## names for the drug-sensitive cell lines

  for(i in c('data', 'sensitive')) {

    expl_res$long_data[[i]] <-
      map2(expl_res[[i]],
           expl_res$variables,
           ~pivot_longer(.x,
                         cols = all_of(.y),
                         names_to = 'variable',
                         values_to = 'resistance'))

  }

  expl_res$long_data <- expl_res$long_data %>%
    transpose %>%
    map(reduce,
        left_join,
        by = c('sample_id', 'FGFR3_mutation', 'variable')) %>%
    map(set_names,
        c('sample_id', 'FGFR3_mutation',
          'variable', 'resistance', 'status')) %>%
    map(mutate,
        name_lab = ifelse(status == 'sensitive',
                          exchange(sample_id,
                                   depmap$cell_lines,
                                   key = 'sample_id',
                                   value = 'cell_line_name'),
                          NA),
        data_set = exchange(variable,
                            expl_res$lexicon,
                            value = 'data_set'),
        label = exchange(variable,
                         expl_res$lexicon),
        label = stri_replace(label, fixed = ', ', replacement = '\n'))

# Box plots -------

  insert_msg('Box plots')

  expl_res$box_plots <-
    list(x = expl_res$long_data,
         v = c(10, 1),
         u = c('resistance<br>vs max. concentration',
               'resistance<br>vs DMSO'),
         w = c('log10', 'identity'),
         y = c('resistance, IC50, µM',
               'resistance vs DMSO, AUC'),
         z = c('GDSC drug screening',
               'PRISM drug screening')) %>%
    pmap(function(x, v, u, w, y, z) x %>%
           ggplot(aes(x = resistance,
                      y = reorder(label,
                                  -resistance,
                                  FUN = function(x) median(x, na.rm = TRUE)),
                      color = status)) +
           facet_grid(data_set ~ .,
                      scales = 'free',
                      space = 'free') +
           geom_vline(xintercept = v,
                      linetype = 'dashed') +
           geom_boxplot(outlier.color = NA,
                        alpha = 0.25,
                        color = 'black',
                        linewidth = 0.5,
                        fill = 'cornsilk') +
           geom_point(aes(shape = FGFR3_mutation,
                          size = FGFR3_mutation,
                          alpha = FGFR3_mutation,
                          fill = status),
                      position = position_jitter(0, 0.1)) +
           #geom_text_repel(aes(label = name_lab),
            #               size = 2.3,
             #              show.legend = FALSE) +
           scale_shape_manual(values = c(WT = 8,
                                         mutated = 23),
                              name = html_italic('FGFR3')) +
           scale_size_manual(values = c(WT = 2,
                                        mutated = 3),
                             name = html_italic('FGFR3')) +
           scale_alpha_manual(values = c(WT = 0.5,
                                         mutated = 1),
                              name = html_italic('FGFR3')) +
           scale_color_manual(values = c(sensitive = 'steelblue',
                                         resistant = 'firebrick'),
                              name = '') +
           scale_fill_manual(values = c(sensitive = 'steelblue',
                                        resistant = 'firebrick'),
                             name = '') +
           scale_x_continuous(trans = w) +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 legend.title = element_markdown()) +
           labs(title = z,
                x = y))

# END --------

  expl_res <- expl_res[c("sensitive", "box_plots")]

  insert_tail()
