# Networks of mutations, deletions, and amplifications.
#
# The affinity matrices are constructed with Jaccard similarity coefficient J,
# and the cutoff of J >= 0.2.


  insert_head()

# container -------

  expl_alnet <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## attribute table:

  expl_alnet$attr_table <-
    globals[c("ligands", "receptors", "binding_proteins")] %>%
    map(~tibble(variable = .x)) %>%
    compress(names_to = 'protein_type') %>%
    mutate(protein_type = stri_replace(protein_type,
                                       regex = 's$',
                                       replacement = ''),
           protein_type = stri_replace(protein_type,
                                       fixed = '_',
                                       replacement = ' '),
           protein_type = factor(protein_type,
                                 c('receptor', 'ligand', 'binding protein')))

  expl_alnet$variables <- expl_alnet$attr_table$variable

  expl_alnet$attr_table <-
    list(x = c('_mut', '_amp', '_del'),
         y = c('mutation', 'amplification', 'deletion')) %>%
    pmap_dfr(function(x, y) expl_alnet$attr_table %>%
               mutate(gene_symbol = variable,
                      variable = paste0(variable, x),
                      alteration = y))

  ## analysis data: binary-coded presence of alterations

  for(i in list(c('mutation', '_mut'),
                c('amplification', '_amp'),
                c('deletion', '_del'))) {

    expl_alnet$data[[i[[1]]]] <- globals$cohort_expr %>%
      eval %>%
      map(~.x[[i[[1]]]]) %>%
      compact %>%
      map(column_to_rownames, 'sample_id') %>%
      map(select, any_of(expl_alnet$variables)) %>%
      map(~set_names(.x,
                     paste0(names(.x), i[[2]])))

  }

  expl_alnet$data <- expl_alnet$data %>%
    transpose %>%
    map(compact) %>%
    map(reduce, cbind)

  ## cohort-specific variable vectors

  expl_alnet$variables <- expl_alnet$data %>%
    map(names)

# Jaccard similarity and construction of graphs --------

  insert_msg('Jaccard simlarity and graphs')

  expl_alnet$simil_mtx <- expl_alnet$data %>%
    map(t) %>%
    map(calculate_dist, method = 'jaccard') %>%
    map(~1 - .x)

  expl_alnet$graph_obj <- expl_alnet$simil_mtx %>%
    map(as_iGraph,
        input_type = 'similarity',
        weighted = TRUE,
        diag = FALSE,
        cutoff = 0.2) %>%
    map(prune_degree, cutoff = 0)

# Vertex importance stats -----

  insert_msg('Vertex importance stats')

  expl_alnet$stats <- expl_alnet$graph_obj %>%
    map(summary) %>%
    map(select, -index)

# Setting the graph attributes: vertex stats and protein type -------

  insert_msg('Graph attributes')

  expl_alnet$graph_obj <-
    map2(expl_alnet$graph_obj,
         expl_alnet$stats,
         set_vertex_attributes) %>%
    map(set_vertex_attributes,
        expl_alnet$attr_table)

# Network plots --------

  insert_msg('Network plots')

  expl_alnet$plots <-
    list(x = expl_alnet$graph_obj,
         vertex_shape_variable = 'protein_type',
         plot_title = globals$cohort_labs[names(expl_alnet$graph_obj)]) %>%
    pmap(plot,
         layout = layout.fruchterman.reingold,
         weighting_order = 1,
         label_vertices = TRUE,
         vertex_label_variable = 'gene_symbol',
         vertex_fill_variable = 'alteration',
         vertex_color_variable = NULL,
         vertex_txt_face = 'italic',
         vertex_size = 3,
         vertex_color = 'black',
         cust_theme = globals$net_theme,
         seed = 12345) %>%
    map(~.x +
          guides(fill = guide_legend(override.aes = list(shape = 21,
                                                         size = 3))) +
          scale_linewidth_continuous(limits = c(0.2, 1),
                                     range = c(0.25, 1),
                                     name = 'J') +
          scale_alpha_continuous(limits = c(0.2, 1),
                                 range = c(0.25, 1),
                                 name = 'J') +
          scale_shape_manual(values = c(21, 22, 23),
                             name = 'Protein type',
                             drop = FALSE) +
          scale_fill_manual(values = c(mutation = globals$mut_status_color[["mutated"]],
                                       amplification = globals$amp_status_color[["amplified"]],
                                       deletion = globals$del_status_color[["deleted"]]),
                            name = 'Alteration type',
                            drop = FALSE))

# END -------

  expl_alnet <- expl_alnet[c("simil_mtx", "graph_obj", "stats", "plots")]

  insert_tail()
