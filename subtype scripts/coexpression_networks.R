# Co-expression networks of FGF-, FGFR-, and FGFBP-coding genes found to be
# differentially regulated between consensus molecular classes in at least
# two of the TCGA, IMvigor, and BCAN cohorts.

  insert_head()

# container -----

  sub_exnet <- list()

# analysis data --------

  insert_msg('Analysis data')

  ## variables of interest and their regulation status
  ## in the consensus subtypes, protein product classification

  sub_exnet$attr_table <- sub_dge$common_significant %>%
    map(map,
        ~tibble(variable = .x)) %>%
    map2(., names(.),
         ~compress(.x, names_to = .y)) %>%
    reduce(full_join, by = 'variable')

  sub_exnet$attr_table[-1] <- sub_exnet$attr_table[-1] %>%
    map_dfc(~ifelse(is.na(.x), 'ns', .x)) %>%
    map_dfc(factor, c('upregulated', 'downregulated', 'ns'))

  sub_exnet$attr_table <- sub_exnet$attr_table %>%
    mutate(protein_type = ifelse(variable %in% globals$ligands,
                                 'ligand',
                                 ifelse(variable %in% globals$receptors,
                                        'receptor', 'binding protein')),
           protein_type = factor(protein_type,
                                 c('receptor', 'ligand', 'binding protein')))

  sub_exnet$sub_levels <-
    names(sub_exnet$attr_table)[!names(sub_exnet$attr_table) %in% c('protein_type', 'variable')]

  sub_exnet$variables <- sub_exnet$attr_table$variable

  ## log2 expression data

  sub_exnet$data <- list(tcga = tcga,
                         imvigor = imvigor,
                         bcan = bcan) %>%
    map(~.x$expression[c('sample_id', sub_exnet$variables)]) %>%
    map(column_to_rownames, 'sample_id')

# Graph objects ---------

  insert_msg('Graph objects')

  sub_exnet$graph_obj <- sub_exnet$data %>%
    map(as_iGraph,
        fun = cor,
        cutoff = 0.3,
        method = 'spearman',
        weighted = TRUE,
        diag = FALSE)

# Community stats --------

  insert_msg('Community stats')

  sub_exnet$stats <- sub_exnet$graph_obj %>%
    map(summary) %>%
    map(select, -index)

# Setting the node attributes: nod stats and regulation status --------

  insert_msg('Setting the vertex attributes')

  sub_exnet$graph_obj <-
    map2(sub_exnet$graph_obj,
         sub_exnet$stats,
         set_vertex_attributes) %>%
    map(set_vertex_attributes,
        sub_exnet$attr_table)

# Plots ---------

  insert_msg('Plots')

  for(i in names(sub_exnet$graph_obj))

    sub_exnet$plots[[i]] <-
      list(vertex_fill_variable = sub_exnet$sub_levels,
           plot_title = sub_exnet$sub_levels %>%
             paste(globals$cohort_labs[[i]], sep = ' vs LumP, ')) %>%
      pmap(plot,
           x = sub_exnet$graph_obj[[i]],
           layout = layout.fruchterman.reingold,
           weighting_order = 1,
           label_vertices = TRUE,
           vertex_shape_variable = 'protein_type',
           vertex_color_variable = NULL,
           vertex_txt_face = 'italic',
           vertex_size = 3,
           vertex_color = 'black',
           cust_theme = globals$net_theme,
           seed = 12345) %>%
      map(~.x +
            guides(fill = guide_legend(override.aes = list(shape = 21,
                                                           size = 3))) +
            scale_linewidth_continuous(limits = c(0.3, 1),
                                       range = c(0.25, 1),
                                       name = '\u03C1') +
            scale_alpha_continuous(limits = c(0.3, 1),
                                   range = c(0.25, 1),
                                   name = '\u03C1') +
            scale_shape_manual(values = c(21, 22, 23),
                               name = 'Protein type',
                               drop = FALSE) +
            scale_fill_manual(values = c(upregulated = 'firebrick',
                                         downregulated = 'steelblue',
                                         ns = 'gray75'),
                              name = 'regulation\nvs LumP')) %>%
      set_names(sub_exnet$sub_levels)

# END -------

  rm(i)

  sub_exnet <- sub_exnet[c("graph_obj", "stats", "plots")]

  insert_tail()
