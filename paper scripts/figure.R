# Elements of Figure 1 of the short correspondence.
# the figure is subsequently stitched by hand.

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
              caption = paste('Alterations of FGFR/FGF/FGFBP-coding',
                              'genes in consensus molecular classes of MIBC.'),
              w = 150,
              h = 110)

# Figure 1B: co-expression networks ---------

  insert_msg('Figure 1B: co-expression networks')

  ## showing only the TCGA cohort, the rest for the supplements

  fig$figure_1b <-
    plot_grid(expl_exnet$paper_plots$tcga +
                theme(legend.position = 'none'),
              ggdraw(),
              get_legend(expl_exnet$paper_plots$tcga),
              ncol = 3,
              rel_widths = c(1, 0.07, 0.35)) %>%
    as_figure(label = 'figure_1b_fgfr_fgf_coexpression_networks',
              ref_name = 'figure_1b',
              caption = paste('Co-expression networks of',
                              'FGFR/FGF/FGFBP-coding genes.'),
              w = 150,
              h = 110)

# Figure 1C: differential gene expression in consensus molecular subsets -------

  insert_msg('Figure 1E: differential gene expression, consensus classes')

  ## showing the TCGA; the res to be shown in supplements

  fig$figure_1c <-
    as_figure(sub_dge$hm_plots$tcga +
                theme(plot.title = element_text(size = 10),
                      plot.subtitle = globals$figure_text,
                      legend.title = globals$figure_text,
                      legend.text = globals$figure_text,
                      axis.title.x = globals$figure_text,
                      axis.text.y = element_markdown(size = 10),
                      strip.text.x = globals$figure_text,
                      legend.position = 'bottom',
                      plot.margin = ggplot2::margin(t = 3, unit = 'mm')),
              label = 'figure_1c_differential_expression_consensus_subsets',
              ref_name = 'figure_1c',
              caption = paste('Differential expression of FGFR/FGF/FGFR-coding',
                              'genes in consensus molecular classes of MIBC.'),
              w = 140,
              h = 150)

# Figure 1D: SHAP importance, Elastic Net model, FGFR/FGF/FGFB independent variables -----------

  insert_msg('Figure 1D: SHAP importances, Elastic Net')

  ## SHAP importance for Elastic Net,
  ## importance for the Random Forest model goes into the supplements

  fig$figure_1d <- ml_splots$paper_plots$elnet$panels %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(ml_splots$paper_plots$elnet$legend,
              nrow = 2,
              rel_heights = c(0.87, 0.13)) %>%
    as_figure(label = 'figure_1d_elastic_net_model_mibc_classes_shap',
              ref_name = 'figure_1d',
              caption = paste('SHAP variable importance in the Elastic Net',
                              'model of consensus molecular classes of MIBC',
                              'with FGFR/FGF/FGFBP-coding genes as explanatory',
                              'factors.'),
              w = 160,
              h = 150)

# Figure 1E:  PRISM drug screening ------

  insert_msg('Figure 1E: PRISM drug screening results')

  ## violin plots, according to the rule 21 of the figure/table guidelines
  ## Results of GDSC screening go into Supplementary Material

  fig$figure_1e <- expl_res$violin_plots$prism +
    scale_y_discrete(labels = drug_labeller) +
    globals$figure_theme +
    theme(legend.title = element_markdown(),
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical")

  fig$figure_1e <- fig$figure_1e %>%
    as_figure(label = 'figure_1e_fgfr_pan_inhibitors_reistance_cell_lines',
              ref_name = 'resistance',
              caption = paste('Sensitivity to pan-FGFR',
                              'inhibitors of DepMap urothelial cancer',
                              'cell lines.'),
              w = 100,
              h = 100)

# Figure 1F: resistance and sensitivity, cell lines WT and FGFR3 mutant ------

  insert_msg('Figure 1F: cell lines, FGFR3 mutation/WT')

  fig$figure_1f <-
    plot_grid(iv_plots$line_plot +
                globals$figure_theme +
                theme(legend.position = "none",
                      strip.text.x = element_markdown())) %>%
    as_figure(label = 'figure_1f_in_vitro_cell_lines_erdafitinib',
              ref_name = 'figure_1f',
              caption = paste('Confluence of cultured FGFR1-4 wild-type',
                              '(UM-UC-3, 5637 and RT112) and FGFR3 mutated',
                              '(UM-UC-14 and UM-UC-6) UC cell lines after',
                              'treatment with different concentrations of',
                              'erdafitinib.'),
              w = 200,
              h = 100)

# The complete Figure 1 ---------

  insert_msg('The complete Figure 1')

  ## upper panels A and B

  fig$figure_1_complete$upper <-
    plot_grid(fig$figure_1a$plot,
              fig$figure_1b$plot,
              ncol = 2,
              rel_widths = c(fig$figure_1a$w, fig$figure_1b$w),
              labels = c('A', 'B'),
              label_size = 14)

  ## middle panels C and D

  fig$figure_1_complete$middle <-
    plot_grid(fig$figure_1c$plot,
              fig$figure_1d$plot,
              ncol = 2,
              rel_widths = c(fig$figure_1c$w, fig$figure_1d$w),
              labels = c('C', 'D'),
              label_size = 14)

  ## bottom panels E and F

  fig$figure_1_complete$bottom <-
    plot_grid(fig$figure_1e$plot,
              fig$figure_1f$plot,
              ncol = 2,
              rel_widths = c(fig$figure_1e$w, fig$figure_1f$w),
              labels = c('E', 'F'),
              label_size = 14)

  ## the entire figure

  fig$figure_1_complete <-
    plot_grid(fig$figure_1_complete$upper,
              fig$figure_1_complete$middle,
              fig$figure_1_complete$bottom,
              nrow = 3,
              rel_heights = c(fig$figure_1a$h,
                              fig$figure_1c$h,
                              fig$figure_1e$h)) %>%
    as_figure(label = 'figure_1_complete',
              ref_name = 'figure_1_complete',
              caption = paste('Distinct FGFR signaling activation mechanisms',
                              'specific for consensus molecular classes of',
                              'MIBC.'),
              w = fig[c("figure_1a", "figure_1b")] %>%
                map_dbl(~.x$w) %>%
                sum,
              h = fig[c("figure_1a", "figure_1c", "figure_1e")] %>%
                map_dbl(~.x$h) %>%
                sum)

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
