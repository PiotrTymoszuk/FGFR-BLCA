# Networks of co-expression of FGF ligands and FGFR receptors.
#
# Spearman's correlation-weighted graph objects are created for gene - gene
# correlations of rho >= 0.3.
# Identification of the hubs: top hub score as a metric of overall correlation
# strength.

  insert_head()

# container -------

  expl_exnet <- list()

# analysis variables and data -------

  insert_msg('Analysis variables and data')

  ## attribute table: indicates if a gene codes for a ligand
  ## or a receptor

  expl_exnet$attr_table <-
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
                                 c('receptor', 'ligand', 'binding protein')))

  expl_exnet$variables <- expl_exnet$attr_table$gene_symbol

  ## analysis data and a cohort-specific list of variables

  expl_exnet$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, any_of(expl_exnet$variables))

  expl_exnet$variables <- expl_exnet$data %>%
    map(names) %>%
    map(intersect, expl_exnet$variables)

# Graph objects -------

  insert_msg('Graph objects')

  expl_exnet$graph_obj <- expl_exnet$data %>%
    map(as_iGraph,
        fun = cor,
        cutoff = 0.3,
        method = 'spearman',
        weighted = TRUE,
        diag = FALSE) %>%
    map(prune_degree, cutoff = 0)

# Node importance stats --------

  insert_msg('Node importance stats')

  expl_exnet$stats <- expl_exnet$graph_obj %>%
    map(summary) %>%
    map(select, -index)

# appending the graphs with attributes and node importance stats -------

  insert_msg('Setting the node importance')

  expl_exnet$graph_obj <-
    map2(expl_exnet$graph_obj,
         expl_exnet$stats,
         set_vertex_attributes)

  expl_exnet$graph_obj <- expl_exnet$graph_obj %>%
    map(set_vertex_attributes,
        expl_exnet$attr_table)

# Plots ---------

  insert_msg('Plots')

  expl_exnet$plots <-
    list(x = expl_exnet$graph_obj,
         vertex_shape_variable = 'protein_type',
         plot_title = globals$cohort_labs[names(expl_exnet$graph_obj)]) %>%
    pmap(plot,
         layout = layout.fruchterman.reingold,
         weighting_order = 1,
         label_vertices = TRUE,
         vertex_fill_variable = 'hub_score',
         vertex_color_variable = NULL,
         vertex_txt_face = 'italic',
         vertex_size = 3,
         vertex_color = 'black',
         cust_theme = globals$net_theme,
         seed = 12345) %>%
    map(~.x +
          scale_linewidth_continuous(limits = c(0.3, 1),
                                     range = c(0.25, 1),
                                     name = '\u03C1') +
          scale_alpha_continuous(limits = c(0.3, 1),
                                 range = c(0.05, 0.35),
                                 name = '\u03C1') +
          scale_shape_manual(values = c(21, 22, 23),
                             name = 'Protein type',
                             drop = FALSE) +
          scale_fill_gradient2(low = 'steelblue',
                               mid = 'black',
                               high = 'firebrick',
                               limits = c(0, 1),
                               midpoint = 0.5,
                               name = 'Hub score'))

# END -------

  expl_exnet <- expl_exnet[c("graph_obj", "stats", "plots")]

  insert_tail()
