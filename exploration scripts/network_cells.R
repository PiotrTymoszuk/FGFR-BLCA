# correlation and network analysis of gene expression variables in the
# GDSC cell line data set.

  insert_head()

# container ------

  expl_celnet <- list()

# analysis globals --------

  insert_msg('analysis globals')

  expl_celnet$variables <- depmap$lexicon %>%
    filter(type == 'expression') %>%
    .$variable %>%
    intersect(expl_cells$mod_variables)

  expl_celnet$data <- depmap$expression %>%
    select(all_of(expl_celnet$variables)) %>%
    filter(complete.cases(.))

  ## attribute data frame

  expl_celnet$attr_table <-
    globals[c("ligands", "receptors", "binding_proteins")] %>%
    map(~tibble(gene_symbol = .x)) %>%
    compress(names_to = 'protein_type') %>%
    mutate(protein_type = stri_replace(protein_type,
                                       regex = 's$',
                                       replacement = ''),
           protein_type = stri_replace(protein_type,
                                       fixed = '_',
                                       replacement = ' '),
           protein_type = factor(protein_type,
                                 c('receptor', 'ligand', 'binding protein'))) %>%
    filter(gene_symbol %in% expl_celnet$variables)

# graph object and communities --------

  insert_msg('Graph object and communities')

  expl_celnet$graph_obj <- expl_celnet$data %>%
    as_iGraph(fun = cor,
              cutoff = 0.3,
              method = 'spearman',
              weighted = TRUE,
              diag = FALSE) %>%
    prune_degree(cutoff = 0)

  set.seed(1234)

  expl_celnet$communities <- expl_celnet$graph_obj %>%
    cluster_leiden(objective_function = 'modularity',
                   n_iterations = 100,
                   resolution_parameter = 0.5) %>%
    comm_recode(c('FGFR2/3' = '1',
                  'FGFR1' = '2'))

  expl_celnet$graph_obj <- add_communities(expl_celnet$graph_obj,
                                           expl_celnet$communities)

# Calculation of vertex importance stats in the communities --------

  insert_msg('Calculation of vertex importance stats for the communities')

  expl_celnet$stats <- c('FGFR2/3' = 'FGFR2/3',
                         'FGFR1' = 'FGFR1') %>%
    map(~select_vertices(expl_celnet$graph_obj, community_id == .x)) %>%
    map(summary) %>%
    map(select, -index)

# appending the graph objects with importance stats and protein types -------

  insert_msg('Setting the attributes')

  expl_celnet$graph_obj <- expl_celnet$graph_obj %>%
    set_vertex_attributes(expl_celnet$attr_table) %>%
    set_vertex_attributes(reduce(expl_celnet$stats, rbind))

# Plots --------

  insert_msg('Plots')

  ## community assignment

  expl_celnet$plots$communities <- expl_celnet$graph_obj %>%
    plot(layout = layout.fruchterman.reingold,
         weighting_order = 1,
         label_vertices = TRUE,
         vertex_fill_variable = 'community_id',
         vertex_color_variable = NULL,
         vertex_shape_variable = 'protein_type',
         vertex_txt_face = 'italic',
         vertex_size = 3,
         vertex_color = 'black',
         cust_theme = globals$net_theme,
         seed = 12345,
         plot_title = 'Gene communities, DepMap cell lines') +
    scale_linewidth_continuous(limits = c(0.3, 1),
                               range = c(0.25, 1),
                               name = '\u03C1') +
    scale_alpha_continuous(limits = c(0.3, 1),
                           range = c(0.05, 0.35),
                           name = '\u03C1') +
    scale_fill_manual(values = c('aquamarine3',
                                 'plum4'),
                      name = 'Gene community') +
    scale_color_manual(values = c('aquamarine3',
                                  'plum4'),
                      name = 'Gene community') +
    scale_shape_manual(values = c(21, 22, 23),
                       name = 'Protein type',
                       drop = FALSE) +
    guides(shape = guide_legend(override.aes = list(fill = 'cornsilk2')),
           fill = guide_legend(override.aes = list(shape = 21)))

  ## hub scores within the community

  expl_celnet$plots$hubs <- expl_celnet$graph_obj %>%
    plot(layout = layout.fruchterman.reingold,
         weighting_order = 1,
         label_vertices = TRUE,
         vertex_fill_variable = 'hub_score',
         vertex_color_variable = NULL,
         vertex_shape_variable = 'protein_type',
         vertex_txt_face = 'italic',
         vertex_size = 3,
         vertex_color = 'black',
         cust_theme = globals$net_theme,
         seed = 12345,
         plot_title = 'Hub genes, DepMap cell lines') +
    scale_linewidth_continuous(limits = c(0.3, 1),
                               range = c(0.25, 1),
                               name = '\u03C1') +
    scale_alpha_continuous(limits = c(0.3, 1),
                           range = c(0.05, 0.35),
                           name = '\u03C1') +
    scale_fill_gradient2(low = 'steelblue',
                         mid = 'black',
                         high = 'firebrick',
                         limits = c(0, 1),
                         midpoint = 0.5,
                         name = 'Hub score<br>within community') +
    scale_shape_manual(values = c(21, 22, 23),
                       name = 'Protein type',
                       drop = FALSE) +
    guides(shape = guide_legend(override.aes = list(fill = 'cornsilk2')),
           fill = guide_legend(override.aes = list(shape = 21)))


# END -------

  expl_celnet <- expl_celnet[c("graph_obj", "stats", "plots")]

  insert_tail()
