# elements of figure 1 of the short correspondence.
# the figure is subsequently stitched by hand

  insert_head()

# container -------

  fig <- list()

# Figure 1A: oncoplot, FGFR and FGF alterations, TCGA cohort -------

  insert_msg('Figure 1A: genetic alterations, FGFR/FGF genes, oncoplot')

  ## genetic alterations of interest in consensus classes of MIBC
  ## not assigned and NE-like cancers are hidden

  fig$figure_1a <-
    as_figure(sub_genet$onco_plots$tcga +
                labs(subtitle = sub_genet$onco_plots$tcga$labels$subtitle %>%
                       stri_replace(regex = ',\\s{1}NE.*',
                                    replacement = '')) +
                theme(legend.position = 'bottom',
                      legend.text = globals$figure_text,
                      legend.title = element_blank(),
                      axis.title.x = globals$figure_text,
                      axis.text.y = element_markdown(size = 10),
                      plot.subtitle = globals$figure_text,
                      plot.title = element_text(size = 10, face = 'bold'),
                      strip.text.x = globals$figure_text),
              label = 'figure_1a_fgfr_fgf_mutations_oncoplot_tcga',
              ref_name = 'figure_1a',
              caption = paste('Alterations of FGFR-, FGF-, and FGFBP-coding',
                              'genes in MIBC consensus molecular classes of',
                              'urothelial cancer.'),
              w = 180,
              h = 110)

# Figure 1B: mutation density and mutation hotspots, FGFR3 ------

  insert_msg('Figure 1B: mutation density and mutation hotspots, FGFR3')

  ## density plot and detailed plot

  fig$figure_1b <-
    plot_grid(expl_hot$domain_density_plots_paper$FGFR3,
              expl_hot$domain_percentage_plots_paper$protein_domain$FGFR3,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_fgfr$domain_plots$FGFR3 +
                           globals$figure_theme +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(margin =
                                                              ggplot2::margin(r = 10)))),
              nrow = 2,
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_1b_fgfr3_mutation_hotspots',
              ref_name = 'figure_1b',
              caption = 'Mutational hot spots of the FGFR3 gene.',
              w = 190,
              h = 160)

# Figure 1C: mRNA and protein expression of the genes of interest --------

  insert_msg('Figure 1C: mRNA and protein expression')

  ## left: heat map with the mRNA expression Z scores
  ## right: box plots with IHVC scores

  fig$figure_1c <-
    list(x = list(expl_expr$hm_plot, expl_expr$plots$hpa),
         y = c('mRNA expression', 'Protein expression'),
         z = c('italic', 'plain'),
         v = c(45, 0)) %>%
    pmap(function(x, y, z, v) x +
           guides(x = guide_axis(angle = v)) +
           labs(title = y) +
           globals$figure_theme +
           theme(plot.subtitle = element_blank(),
                 axis.text.y = element_text(face = z),
                 axis.title.y = element_blank()))

  fig$figure_1c[[1]] <- fig$figure_1c[[1]] +
    theme(axis.title.x = element_blank())

  fig$figure_1c <-
    plot_grid(fig$figure_1c[[1]] +
                theme(legend.position = 'none'),
              plot_grid(fig$figure_1c[[2]] +
                          theme(legend.position = 'none'),
                        get_legend(fig$figure_1c[[1]] +
                                     theme(legend.position = 'bottom')),
                        nrow = 2,
                        rel_heights = c(0.73, 0.27)),
              ncol = 2,
              rel_widths = c(1, 1.15))

  fig$figure_1c <- fig$figure_1c %>%
    as_figure(label = 'figure_1c_mRNA_protein_expression',
              ref_name = 'figure_1c',
              caption = paste('Expression of genes coding FGFR, FGF, and FGFBP',
                              'at mRNA and protein level.'),
              w = 180,
              h = 183)

# Figure 1D: co-expression networks ---------

  insert_msg('Figure 1D: co-expression networks')

  fig$figure_1d <- expl_exnet$paper_plots %>%
    map(~.x +
          theme(plot.subtitle = element_blank(),
                legend.position = 'none')) %>%
    c(list(get_legend(expl_exnet$paper_plots[[1]] +
                        theme(legend.box = 'horizontal',
                              legend.text = element_text(margin = ggplot2::margin(r = 10)))))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_1d_fgfr_fgf_coexpression_networks',
              ref_name = 'figure_1d',
              caption = paste('Co-expression networks of FGFR-, FGF-,',
                              'and FGFBP-coding genes.'),
              w = 190,
              h = 190)

# Figure 1E: differential gene expression in consensus molecular subsets -------

  insert_msg('Figure 1E: differential gene expression, consensus classes')

  fig$figure_1e <-
    as_figure(sub_dge$hm_plots$tcga +
                theme(plot.title = element_text(size = 10),
                      plot.subtitle = globals$figure_text,
                      legend.title = globals$figure_text,
                      legend.text = globals$figure_text,
                      axis.title.x = globals$figure_text,
                      axis.text.y = element_markdown(size = 10),
                      strip.text.x = globals$figure_text),
              label = 'figure_1e_differential_expression_consensus_subsets',
              ref_name = 'figure_1e',
              caption = paste('Differential expression of FGFR-, FGF-, and',
                              'FGFR-coding genes in the consensus molecular',
                              'classes of urothelial cancers.'),
              w = 140,
              h = 115)

# Figure 1F: SHAP importance -----------

  insert_msg('Figure 1F: SHAP importances')

  ## recommended: to show importance for just one model
  ## hence, we're preparing two figure panels

  fig$figure_1f_elnet <- ml_splots$paper_plots$elnet$panels %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(ml_splots$paper_plots$elnet$legend,
              nrow = 2,
              rel_heights = c(0.92, 0.1))

  fig$figure_1f_ranger <- ml_splots$paper_plots$ranger$panels %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(ml_splots$paper_plots$ranger$legend,
              nrow = 2,
              rel_heights = c(0.92, 0.1))

  ## figure objects

  fig[c("figure_1f_elnet", "figure_1f_ranger")] <-
    fig[c("figure_1f_elnet", "figure_1f_ranger")] %>%
    list(label = c('figure_1f_shap_elastic_net_model',
                   'figure_1f_shap_random_forest_model'),
         ref_name = names(.),
         caption = paste('Importance of explanatory variables in the',
                         c('Elastic Net', 'Random Forest'),
                         'model of MIBC molecular consensus classes of',
                         'urothelial cancers assessed by the',
                         'SHAP algorithm.')) %>%
    pmap(as_figure,
         w = 190,
         h = 170)

# Figure 1G:  GDSC and PRISM drug screening ------

  insert_msg('Figure 1G: drug screening results')

  fig$figure_1g <-
    list(x = expl_res$box_plots[c("gdsc", "prism")],
         y = list(element_rect(), element_blank()),
         z = list(element_text(), element_blank())) %>%
    pmap(function(x, y, z) x +
           scale_y_discrete(labels = function(x) {

             x %>%
               stri_replace(regex = '\\n.*', replacement = '') %>%
               stri_capitalize_first

           }) +
           globals$figure_theme +
           theme(legend.title = element_markdown(),
                 strip.background.y = y,
                 strip.text.y = z,
                 axis.title.y = element_blank()))

  fig$figure_1g <-
    plot_grid(fig$figure_1g[[1]] +
                theme(legend.position = 'none'),
              plot_grid(fig$figure_1g[[2]] +
                          theme(legend.position = 'none'),
                        get_legend(fig$figure_1g[[1]] +
                                     theme(legend.position = 'bottom',
                                           legend.box = 'vertical')),
                        nrow = 2,
                        rel_heights = c(0.6, 0.4)),
              ncol = 2,
              rel_widths = c(1.1, 1)) %>%
    as_figure(label = 'figure_1g_fgfr_pan_inhibitors_reistance_cell_lines',
              ref_name = 'resistance',
              caption = paste('Resistance and sensitivity to pan-FGFR',
                              'inhibitors of DepMap cancer urothelial',
                              'cell lines.'),
              w = 190,
              h = 120)

# Figure 1H: resistance and sensitivity, cell lines WT and FGFR3 mutant ------

  insert_msg('Figure 1H: cell lines, FGFR3 mutation/WT')

  fig$figure_1h <-
    as_figure(ggdraw(),
              label = 'figure_1h_in_vitro_cell_lines_placeholder',
              ref_name = 'figure_1h',
              caption = paste('Confluence of cultured FGFR1-4 wild-type',
                              '(UM-UC-3, 5637 and RT112) and FGFR3 mutated',
                              '(UM-UC-14 and UM-UC-6) UC cell lines after',
                              'treatment with different concentrations of',
                              'erdafitinib.'),
              w = 180,
              h = 180)

# Saving figures on the disc -------

  insert_msg('Saving figure panels on the disc')

  ## as PDF and PNG files

  fig %>%
    walk(pickle,
         format = 'pdf',
         path = './report/paper figure 1',
         device = cairo_pdf)

  fig %>%
    walk(pickle,
         format = 'png',
         path = './report/paper PNG figure 1',
         dpi = 800)

# END -----

  insert_tail()
