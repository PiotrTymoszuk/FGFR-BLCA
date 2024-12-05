# Co-occurrence of the top most frequent somatic mutations, amplifications
# and deletions investigated by a network analysis.
#
# The graphs are created with Jaccard similarities as affinity matrices;
# cutoff of J >= 0.2.

  insert_head()

# container -------

  lca_ocnet <- list()

# analysis data -----

  insert_msg('Analysis globals')

  ## variables and lexicon

  lca_ocnet$variables <- lca_globals$variables

  lca_ocnet$lexicon <- lca_globals$lexicon

  lca_ocnet$data <- lca_globals$data %>%
    map(column_to_rownames, 'sample_id')

# Jaccard similarity and graph objects --------

  insert_msg('Similarity matrices and graph objects')

  lca_ocnet$simil_mtx <- lca_ocnet$data %>%
    map(t) %>%
    map(calculate_dist, method = 'jaccard') %>%
    map(~1 - .x)

  lca_ocnet$graph_obj <- lca_ocnet$simil_mtx %>%
    map(as_iGraph,
        input_type = 'similarity',
        weighted = TRUE,
        diag = FALSE,
        cutoff = 0.2) %>%
    map(prune_degree, cutoff = 0)

# Vertex importance stats -----

  insert_msg('Vertex importance stats')

  lca_ocnet$stats <- lca_ocnet$graph_obj %>%
    map(summary) %>%
    map(select, -index)

# Setting the graph attributes: vertex stats and protein type -------

  insert_msg('Graph attributes')

  lca_ocnet$graph_obj <-
    map2(lca_ocnet$graph_obj,
         lca_ocnet$stats,
         set_vertex_attributes) %>%
    map(set_vertex_attributes,
        lca_ocnet$lexicon)

# Network plots --------

  insert_msg('Network plots')

  lca_ocnet$plots <-
    list(x = lca_ocnet$graph_obj,
         plot_title = globals$cohort_labs[names(lca_ocnet$graph_obj)]) %>%
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
          scale_fill_manual(values = c(mutation = globals$mut_status_color[["mutated"]],
                                       amplification = globals$amp_status_color[["amplified"]],
                                       deletion = globals$del_status_color[["deleted"]]),
                            name = 'Alteration type',
                            drop = FALSE))

# END -------

  lca_ocnet <- lca_ocnet[c("simil_mtx", "graph_obj", "stats", "plots")]

  insert_tail()
