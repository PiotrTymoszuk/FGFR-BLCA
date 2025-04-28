# Supplementary figures

  insert_head()

# container -------

  suppl_figs <- list()

# oncoplots for FGFR/FGF/FGF alterations in consensus classes --------

  insert_msg('Oncoplots, genetic alterations of FGFR/FGF/FGFBP genes')

  suppl_figs$sub_genet <- sub_genet$onco_plots[c("imvigor", "bcan")] %>%
    map(~.x +
          theme(plot.title.position = 'plot',
                legend.position = 'none',
                strip.text.y = element_blank(),
                strip.background.y = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(sub_genet$onco_plots[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(size = 7),
                                 legend.title = element_text(size = 7))),
              nrow = 2,
              rel_heights = c(0.85, 0.15)) %>%
    as_figure(label = 'consensus_subtypes_fgf_fgfr_gene_alterations',
              ref_name = 'sub_genet',
              caption = paste('Alterations of FGFR-, FGF-, and FGFBP-coding',
                              'genes in consensus molecular classes of MIBC.'),
              w = 180,
              h = 90)

# classification of FGFR3 mutation by protein domain -------

  insert_msg('Classification of FGFR3 mutations by protein domain')

  suppl_figs$fgfr3_domains <-
    expl_fgfr$stack_plots[c("protein_domain.genie",
                            "protein_domain.msk",
                            "protein_domain.tcga",
                            "protein_domain.bcan")] %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_fgfr$stack_plots$protein_domain.genie +
                           theme(legend.position = 'bottom',
                                 legend.text = element_text(margin = ggplot2::margin(r = 5,
                                                                                     l = 5))) +
                           guides(fill = guide_legend(nrow = 2))),
              nrow = 2,
              rel_heights = c(0.75, 0.25))

  suppl_figs$fgfr3_domains <- suppl_figs$fgfr3_domains %>%
    as_figure(label = 'fgfr3_mutation_domains',
              ref_name = 'fgfr3_domains',
              caption = paste('Frequency of somatic mutations in FGFR1/2/3/4',
                              'genes split by protein domain'),
              w = 180,
              h = 140)

# Mutation density and mutation hotspots, FGFR3 ------

  insert_msg('Mtation density and mutation hotspots, FGFR3')

  ## density plot and detailed plot

  suppl_figs$fgfr3_hot <-
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
    as_figure(label = 'fgfr3_mutation_hotspots',
              ref_name = 'fgfr3_hot',
              caption = 'Mutational hot spots of the FGFR3 gene.',
              w = 190,
              h = 160)

# mRNA and protein expression of the genes of interest --------

  insert_msg('mRNA and protein expression')

  ## left: heat map with the mRNA expression Z scores
  ## right: box plots with IHVC scores

  suppl_figs$expression <-
    list(x = list(expl_expr$hm_plot, expl_expr$plots$hpa),
         y = c('mRNA expression', 'Protein expression'),
         z = c('italic', 'plain'),
         v = c(45, 0)) %>%
    pmap(function(x, y, z, v) x +
           guides(x = guide_axis(angle = v)) +
           labs(title = y) +
           #globals$figure_theme +
           theme(plot.subtitle = element_blank(),
                 axis.text.y = element_text(face = z),
                 axis.title.y = element_blank()))

  suppl_figs$expression[[1]] <- suppl_figs$expression[[1]] +
    theme(axis.title.x = element_blank())

  suppl_figs$expression <-
    plot_grid(suppl_figs$expression[[1]] +
                theme(legend.position = 'none'),
              plot_grid(suppl_figs$expression[[2]] +
                          theme(legend.position = 'none'),
                        get_legend(suppl_figs$expression[[1]] +
                                     theme(legend.position = 'bottom')),
                        nrow = 2,
                        rel_heights = c(0.73, 0.27)),
              ncol = 2,
              rel_widths = c(1, 1.15))

  suppl_figs$expression <- suppl_figs$expression %>%
    as_figure(label = 'expression_mRNA_protein_expression',
              ref_name = 'expression',
              caption = paste('Expression of genes coding FGFR, FGF, and FGFBP',
                              'at mRNA and protein level.'),
              w = 180,
              h = 183)

# Representative staining for FGFR/FGF/FGFBP -------

  insert_msg('Representative IHC stainings')

  suppl_figs[c('ihc_receptors', 'ihc_ligands', 'ihc_fgfbp')] <-
    list(expl_ihc$panels[c("FGFR1", "FGFR2", "FGFR3", "FGFR4")],
         expl_ihc$panels[c("FGF7", "FGF10")],
         expl_ihc$panels[c("SDC2", "GPC4", "CD44", "FGFBP1")]) %>%
    map(map, eval) %>%
    map(~plot_grid(plotlist = .,
                   ncol = 1))
  suppl_figs[c('ihc_receptors', 'ihc_ligands', 'ihc_fgfbp')] <-
    suppl_figs[c('ihc_receptors', 'ihc_ligands', 'ihc_fgfbp')] %>%
    list(x = .,
         label = names(.),
         ref_name = names(.),
         caption = paste('Representative immunohistochemistry staining images',
                         'for',
                         c('FGFR proteins', 'FGF proteins', 'FGFBP proteins'),
                         'retrieved from the Human Protein Atlas.'),
         h = c(160, 80, 160)) %>%
    pmap(as_figure,
         w = 120)

# Co-expression networks, IMvigor and BCAN ---------

  insert_msg('Co-expression networks, Imvigor and BCAN')

  suppl_figs$bulk_networks <-
    expl_exnet$paper_plots[c("imvigor", "bcan")] %>%
    map(~.x +
          theme(plot.subtitle = element_blank(),
                legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(ggdraw(),
              get_legend(expl_exnet$paper_plots[[1]] +
                           theme(legend.box = 'horizontal',
                                 legend.text = element_text(margin = ggplot2::margin(r = 10)))),
              nrow = 3,
              rel_heights = c(0.7, 0.05, 0.3)) %>%
    as_figure(label = 'fgfr_fgf_coexpression_networks',
              ref_name = 'bulk_networks',
              caption = paste('Co-expression networks of FGFR-, FGF-,',
                              'and FGFBP-coding genes in the IMvigor',
                              'and BCAN cohorts.'),
              w = 190,
              h = 140)

# Co-expression network, cell lines --------

  insert_msg('Co-expression networks, cancer tissue')

  suppl_figs$cell_networks <-
    plot_grid(expl_celnet$plots$communities +
                guides(alpha = 'none',
                       size = 'none',
                       color = 'none',
                       linewidth = 'none') +
                theme(plot.subtitle = element_blank(),
                      legend.position = 'right'),
              expl_celnet$plots$hubs +
                guides(shape = 'none') +
                theme(plot.subtitle = element_blank(),
                      legend.position = 'right'),
              nrow = 2,
              align = 'hv',
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'fgfr_fgf_coexpression_networks_cell_lines',
              ref_name = 'cell_networks',
              caption = paste('Co-expression network of',
                              'FGF-, FGFR-, and FGFBP-coding genes in DepMap',
                              'urothelial cancer cell lines'),
              w = 180,
              h = 230)

# Expression of FGFR/FGF/FGFBP genes in consensus classes of MIBC -------

  insert_msg('Expression of FGFR/FGF/FGFBP genes in the consensus classes')

  suppl_figs$sub_dge <- sub_dge$hm_plots[c("imvigor", "bcan")] %>%
    map(~.x +
          guides(y = guide_axis(n.dodge = 1)) +
          theme(legend.position = 'none',
                plot.margin = ggplot2::margin(r = 2, l = 2,
                                              t = 0, b = 1,
                                              unit = 'mm'))) %>%
    plot_grid(plotlist = .,
              nrow = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(sub_dge$hm_plots[[1]] +
                           theme(legend.position = 'right')),
              ncol = 2,
              rel_widths = c(0.8, 0.2)) %>%
    as_figure(label = 'molecular_subtypes_expression_FGF_FGFR_FGFBP',
              ref_name = 'sub_fgfr_dge',
              caption = paste('Differential expression of FGF-, FGFR-,',
                              'and FGFR-coding genes in the consensus molecular',
                              'classes of MIBC.'),
              w = 180,
              h = 190)

# Expression of FGFR-related gene signatures --------

  insert_msg('Expression of FGFR-related gene signatures in consensus classes')

  suppl_figs$sub_sig <- sub_sig$hm_plots[c("tcga", "imvigor", "bcan")] %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title.position = 'plot',
                axis.text.y = element_markdown(size = 6),,
                axis.title.y = element_text(angle = 90)) +
          labs(y = 'FGFR-related gene signature')) %>%
    c(list(get_legend(sub_sig$hm_plots[[1]]))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'consensus_classes_fgfr_gene_signatures',
              ref_name = 'sub_sig',
              caption = paste('Scores of FGFR-related gene signatures in',
                              'the consensus molecular classes of MIBC.'),
              w = 190,
              h = 210)

# Predicted resistance to pan-FGFR inhibitors, consensus subsets --------

  insert_msg('Predicted response to FGFRi in consensus molecular classes')

  ## upper panel: heat maps

  suppl_figs$sub_drugs$upper <- sub_drugs$hm_plots$fgfr_drugs %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.title.position = 'plot',
                axis.text.y = element_markdown(size = 6))) %>%
    c(list(get_legend(sub_drugs$hm_plots$fgfr_drugs[[1]]))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: box plots for erdafitinib

  suppl_figs$sub_drugs$bottom <- sub_drugs$fgfr_plots %>%
    map(~.x$erdafitinib) %>%
    map2(., names(.),
         ~.x +
           labs(title = paste(.x$labels$title,
                              globals$cohort_labs[.y],
                              sep = ', ')) +
           guides(x = guide_axis(angle = 45)) +
           theme(legend.position = 'none',
                 axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  suppl_figs$sub_drugs <-
    plot_grid(suppl_figs$sub_drugs$upper,
              suppl_figs$sub_drugs$bottom,
              nrow = 2,
              rel_heights = c(2.3, 1),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'consensus_classes_predicted_fgfr_inhibitors',
              ref_name = 'sub_drugs',
              caption = paste('Predicted response to pan-FGFR inhibitors in',
                              'the consensus molecular classes of MIBC.'),
              w = 190,
              h = 230)

# Predicted response to pan-FGFR inhibitors, FGFR3 WT and mutant cancers -------

  insert_msg('Predicted response to pan-FGFR inhibitors, FGFR3 WT and mutant')

  ## upper panel: heat maps

  suppl_figs$fgfr_drugs <- fgfr_drugs$hm_plots$fgfr_drugs %>%
    map(~.x +
          theme(legend.position = 'none',
                strip.text.x = element_text(hjust = 0, angle = 90))) %>%
    c(list(get_legend(fgfr_drugs$hm_plots$fgfr_drugs[[1]]))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'predicted_pan_fgfr_inhibitors',
              ref_name = 'fgfr_drugs',
              caption = paste('Predicted response to pan-FGFR inhibitors in',
                              'the urothelial cancers with and without FGFR3',
                              'mutations.'),
              w = 180,
              h = 190)

# Evaluation of the machine learning models --------

  insert_msg('Evaluation of the machine learning models')

  ## upper panel: global performance stats,
  ## no ensemble model, no training subset

  suppl_figs$ml_eval$upper <- mibc_comp$stat_plot

  suppl_figs$ml_eval$upper$data <-
    suppl_figs$ml_eval$upper$data %>%
    filter(algorithm != 'ensemble',
           cohort != 'tcga_train')

  suppl_figs$ml_eval$upper <-
    plot_grid(suppl_figs$ml_eval$upper,
              ncol = 2,
              rel_widths = c(0.9, 0.1))

  ## bottom panel: Brier squares

  suppl_figs$ml_eval$bottom <-
    mibc_comp$brier_square_plots[c("elnet", "ranger")] %>%
    map(~.x[c('tcga_test', 'imvigor', 'bcan')]) %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(mibc_comp$brier_square_plots[[1]][[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## the entire figure

  suppl_figs$ml_eval <-
    plot_grid(suppl_figs$ml_eval$upper,
              suppl_figs$ml_eval$bottom,
              nrow = 2,
              rel_heights = c(1, 2),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'consensus_classes_machine_learning_evaluation',
              ref_name = 'ml_eval',
              caption = paste('Evaluation of predictions of the consensus',
                              'molecular classes of MIBC made by machine',
                              'learning trained with the genuine consensus',
                              'class-defining genes and the',
                              'FGFR/FGF/FGFBP-coding genes.'),
              w = 190,
              h = 230)

# Machine learning models: confusion matrices -------

  insert_msg('Machine leaning models, confusion matrices')

  suppl_figs$confusion <-
    list(mibc_genes = mibc_eval,
         fgfr_genes = ml_eval) %>%
    map(~.x$confusion_plots[c('elnet', 'ranger')]) %>%
    map(map, ~.x[c('tcga_test', 'imvigor', 'bcan')]) %>%
    map(unlist, recursive = FALSE)

  for(i in names(suppl_figs$confusion)) {

    suppl_figs$confusion[[i]] <- suppl_figs$confusion[[i]] %>%
      map(~.x +
            theme(legend.position = 'none',
                  plot.title.position = 'plot') +
            labs(title = .x$labels$title %>%
                   paste(mibc_comp$model_labels[[i]], sep = ', ') %>%
                   stri_replace(fixed = ' BLCA, test',
                                replacement = ', test') %>%
                   stri_replace(fixed = 'Random Forest',
                                replacement = 'RF')))

  }

  suppl_figs$confusion <- suppl_figs$confusion %>%
    unlist(recursive = FALSE) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', '',
                         '', '', '',
                         'B', '', '',
                         '', '', ''),
              label_size = 10) %>%
    as_figure(label = 'consensus_classes_machine_learning_confusion',
              ref_name = 'confusion',
              caption = paste('Confusion matrices of predictions of the consensus',
                              'molecular classes of MIBC made by machine',
                              'learning trained with the genuine consensus',
                              'class-definingg genes and the',
                              'FGFR/FGF/FGFBP-coding genes.'),
              w = 190,
              h = 230)

# SHAP importance, Random Forest model, FGFR/FGF/FGFB independent variables -----------

  insert_msg('Figure 1D: SHAP importances, Random Forest')

  suppl_figs$rf_shap <- ml_splots$plots$ranger$panels %>%
    map(~.x + guides(x = guide_axis(angle = 45))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(ml_splots$plots$ranger$legend,
              nrow = 2,
              rel_heights = c(0.87, 0.13)) %>%
    as_figure(label = 'random_forest_model_mibc_classes_shap',
              ref_name = 'rf_shap',
              caption = paste('SHAP variable importance in the',
                              'Random Forest model of consensus molecular classes',
                              'of MIBC with FGFR/FGF/FGFBP-coding genes as',
                              'explanatory factors.'),
              w = 180,
              h = 180)

# Variable importance, SHAP, full predictor models ---------

  insert_msg('Variable importance, full predictor sets')

  suppl_figs$shap_elnet <- mibc_splots$plots$elnet$panels %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(mibc_splots$plots$elnet$legend,
              nrow = 2,
              rel_heights = c(0.92, 0.1))

  suppl_figs$shap_ranger <- mibc_splots$plots$ranger$panels %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(mibc_splots$plots$ranger$legend,
              nrow = 2,
              rel_heights = c(0.92, 0.1))

  ## figure objects

  suppl_figs[c("shap_elnet", "shap_ranger")] <-
    suppl_figs[c("shap_elnet", "shap_ranger")] %>%
    list(label = c('shap_elastic_net_model_mibc_defining_genes',
                   'shap_random_forest_model_mibc_defining_genes'),
         ref_name = names(.),
         caption = paste('Importance of explanatory variables in the',
                         c('Elastic Net', 'Random Forest'),
                         'model of MIBC molecular consensus classes of',
                         'urothelial cancers with the full set of',
                         'class-defining genes by Kamoun et al. assessed by',
                         'the SHAP algorithm.')) %>%
    pmap(as_figure,
         w = 190,
         h = 170)


# Variable importance, SHAP, FGFR/FGF/FGFBP genes, full predictor models ------

  insert_msg('SHAP importance, FGFR/FGF/FGFBP genes, full predictor models')

  suppl_figs$shap_fgfr <- mibc_splots$mean_plots %>%
    map(~.x +
          labs(title = .x$labels$title %>%
                 paste(mibc_comp$model_labels["mibc_genes"],
                       sep = ', ')) +
          theme(legend.position = 'none',
                plot.subtitle = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(mibc_splots$mean_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.85, 0.15)) %>%
    as_figure(label = 'shap_mibc_defining_genes_FGFR_FGF_FGFBP',
              ref_name = 'shap_fgfr',
              caption = paste('Importance of FGFR/FGF/FGFBP-coding genes',
                              'for predictions of the MIBC consensus classes',
                              'by machine learning models employing the full',
                              'set of the class-defining genes by Kamoun et al.'),
              w = 190,
              h = 110)

# Response to pan-FGFR inhibitors, GDSC ---------

  insert_msg('Response to pan-FGFRi, GDSC')

  suppl_figs$pani_response <- expl_res$box_plots$gdsc +
    theme(legend.title = element_markdown(),
          axis.title.y = element_blank())

  suppl_figs$pani_response <- suppl_figs$pani_response %>%
    as_figure(label = 'fgfr_pan_inhibitors_reistance_cell_lines_gdsc',
              ref_name = 'pani_response',
              caption = paste('Sensitivity to pan-FGFR',
                              'inhibitors of DepMap urothelial cancer',
                              'cell lines in the GDSC1/2 drug screening',
                              'experiments.'),
              w = 135,
              h = 100)

# In vitro erdafitinib treatments --------

  insert_msg('In vitro erdafitinib treatment')

  ## images of the confluence plots provided by the study team

  suppl_figs$in_vitro <-
    paste0('./report/aux files/',
          c('UMUC3 erdafitinib.PNG',
            '5637 erdafitinib.PNG',
            'RT112 erdafitinib.PNG')) %>%
    map(~ggdraw() + draw_image(.x)) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'erdafitinib_confluence_in_vitro',
              ref_name = 'in_vivo',
              caption = paste('Dose response curves of cultured FGFR1-4',
                              'wild-type cell lines (UMUC3, 5637 and RT112)',
                              'after treatment with increasing concentrations',
                              'of erdafitinib.'),
              w = 190,
              h = 160)

# Saving the supplementary figures on the disc ------

  insert_msg('Saving the supplementary figures')

  suppl_figs <- suppl_figs %>%
    number_figures(prefix = 'Supplementary Figure S')

  suppl_figs[!names(suppl_figs) %in% c('ihc_receptors',
                                       'ihc_ligands',
                                       'ihc_fgfbp')] %>%
    walk(pickle,
         path = './report/supplementary figures',
         format = 'pdf',
         device = cairo_pdf)

  suppl_figs[c("ihc_receptors",
               "ihc_ligands",
               "ihc_fgfbp")] %>%
    walk(pickle,
         path = './report/supplementary figures',
         format = 'png',
         dpi = 600)

# END ------

  insert_tail()
