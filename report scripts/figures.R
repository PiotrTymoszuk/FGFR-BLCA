# Figures for the analysis report

  insert_head()

# container ------

  rep_figs <- list()

# General frequency of mutations, amplification and deletions  --------

  insert_msg('General frequency of genetic alterations')

  ## facetting by the gene (?) and common X scale

  rep_figs$general_freq <- expl_genet$plots %>%
    map(~.x +
         # facet_grid(variable ~ .,
          #           scales = 'free',
           #          space = 'free') +
          scale_x_continuous(limits = expl_genet$plots %>%
                               map(~.x$data$percent) %>%
                               reduce(c) %>%
                               max %>%
                               c(0, .)) +
          theme(plot.subtitle = element_blank(),
                legend.position = 'none',
                strip.background.y = element_blank(),
                strip.text.y = element_blank()))

  rep_figs$general_freq <-
    rep_figs$general_freq[c("deletion", "amplification")] %>%
    plot_grid(plotlist = .,
              nrow = 2,
              rel_heights = c(3.5, 8)) %>%
    plot_grid(rep_figs$general_freq$mutation,
              .,
              ncol = 2) %>%
    plot_grid(get_legend(expl_genet$plots$mutation +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  rep_figs$general_freq <- rep_figs$general_freq %>%
    as_figure(label = 'general_frequency_genetic_alterations',
              ref_name = 'general_freq',
              caption = paste('The most frequent somatic mutations and gene',
                              'copy alterations in urothelial cancer.'),
              w = 180,
              h = 180)

# Frequency of alterations of FGF and FGFR genes ------

  insert_msg('Frequency of alterations of FGF and FGFR genes')

  rep_figs$fgfr_freq <- expl_genet$fgf_plots %>%
    map(~.x +
          scale_x_continuous(limits = expl_genet$fgf_plots %>%
                               map(~.x$data$percent) %>%
                               reduce(c) %>%
                               max %>%
                               c(0, .)) +
          theme(plot.subtitle = element_blank(),
                legend.position = 'none')) %>%
    c(list(legend = get_legend(expl_genet$fgf_plots)))

  rep_figs$fgfr_freq <-
    rep_figs$fgfr_freq[c("mutation", "amplification", "deletion", "legend")] %>%
    plot_grid(plotlist = .,
              ncol = 2,
              rel_heights = c(10, 4.5)) %>%
    as_figure(label = 'fgfr_fgf_alteration_frequency',
              ref_name = 'fgfr_freq',
              caption = paste('Frequency of somatic mutations and copy number',
                              'alterations of FGF- and FGFR-coding genes.'),
              w = 180,
              h = 180)

# Classification of mutations of FGFR genes ------

  insert_msg('Mutation classification')

  rep_figs[c('mutation_type',
             'variant_type',
             'protein_domain_split')] <-
    list(expl_fgfr$stack_plots[c("mutation_type.genie", "mutation_type.tcga")],
         expl_fgfr$stack_plots[c("variant_type.genie", "variant_type.tcga")],
         expl_fgfr$stack_plots[c("protein_domain.genie", "protein_domain.tcga")]) %>%
    map(function(x) plot_grid(x[[1]] + theme(legend.position = 'none'),
                              x[[2]] + theme(legend.position = 'none'),
                              ncol = 2) %>%
          plot_grid(get_legend(x[[1]] + theme(legend.position = 'bottom')),
                    nrow = 2,
                    rel_heights = c(2, 1)))

  rep_figs[c('mutation_type',
             'variant_type',
             'protein_domain_split')] <-
    rep_figs[c('mutation_type',
               'variant_type',
               'protein_domain_split')] %>%
    list(x = .,
         label = c('fgfr_mutation_types',
                   'fgfr_variant_types',
                   'fgfr_protein_domain_split'),
         ref_name = names(.),
         caption = paste('Frequency of somatic mutations in FGFR1/2/3/4 genes',
                         'split by',
                         c('mutation type.',
                           'variant type.',
                           'protein domain.'))) %>%
    pmap(as_figure,
         w = 180,
         h = 120)

# Location of mutations in FGFR genes ------

  insert_msg('location of mutations in FGFR genes')

  rep_figs[c('fgfr1_position',
             'fgfr2_position',
             'fgfr3_position',
             'fgfr4_position')] <-
    paste0('FGFR', 1:4) %>%
    map(function(gene) map(expl_fgfr$position_plots, ~.x[[gene]])) %>%
    map2(., list(c(4, 1),
                 c(4, 2),
                 c(9, 5),
                 c(3, 1)),
         function(plots, dims) plots %>%
           plot_grid(plotlist = .,
                     nrow = 2,
                     rel_heights = dims + 2))

  ## appending with the domain legends

  rep_figs[c('fgfr1_position',
             'fgfr2_position',
             'fgfr3_position',
             'fgfr4_position')] <-
    map2(rep_figs[c('fgfr1_position',
                    'fgfr2_position',
                    'fgfr3_position',
                    'fgfr4_position')],
         expl_fgfr$domain_plots,
         ~plot_grid(.x,
                    get_legend(.y),
                    ncol = 2,
                    rel_widths = c(0.8, 0.2)))

  ## the figures

  rep_figs[c('fgfr1_position',
             'fgfr2_position',
             'fgfr3_position',
             'fgfr4_position')] <-
    rep_figs[c('fgfr1_position',
               'fgfr2_position',
               'fgfr3_position',
               'fgfr4_position')] %>%
    list(x = .,
         label = paste0(c('fgfr1', 'fgfr2', 'fgfr3', 'fgfr4'),
                        '_mutations_protein_position'),
         ref_name = names(.),
         caption = paste('Positions of', c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'),
                         'protein affected by somatic mutations.'),
         h = c(160, 170, 230, 140)) %>%
    pmap(as_figure,
         w = 180)

# Co-occurrence of the most common gene alterations --------

  insert_msg('Co-occurrence of the most common genetic features')

  ## upper panel: MDS scatter plots

  rep_figs$genetic_overlap$top <- lca_occur$mds_plots %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv') %>%
    plot_grid(get_legend(lca_occur$mds_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.83, 0.17))

  ## bottom panel: heat maps of contingency tables

  rep_figs$genetic_overlap$bottom <- lca_occur$hm_plots %>%
    map(~.x[c('FGFR3_mutation|RB1_mutation',
              'FGFR3_mutation|CCND1_amplification',
              'FGFR3_mutation|CDKN2A_deletion',
              'TP53_mutation|FGFR3_mutation')]) %>%
    map(map,
        ~.x +
          theme(legend.position = 'none',
                plot.title = element_blank())) %>%
    map(~plot_grid(plotlist = .,
                   ncol = 2,
                   align = 'hv',
                   axis = 'tblr'))

  rep_figs$genetic_overlap$bottom <-
    plot_grid(rep_figs$genetic_overlap$bottom[[1]],
              rep_figs$genetic_overlap$bottom[[2]],
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = globals$cohort_labs[c("genie", "tcga")],
              label_size = 8)

  ## the entire figure

  rep_figs$genetic_overlap <-
    plot_grid(rep_figs$genetic_overlap$top,
              ggdraw(),
              rep_figs$genetic_overlap$bottom,
              nrow = 3,
              rel_heights = c(1.1, 0.1, 1),
              labels = c('A', 'B', ''),
              label_size = 10) %>%
    as_figure(label = 'overlap_most_common_gene_alterations',
              ref_name = 'genetic_overlap',
              caption = paste('Co-occurrence of the most common genetic',
                              'alterations in urothelial cancer.'),
              w = 180,
              h = 210)

# Co-occurrence of FGF and FGFR gene alterations -------

  insert_msg('Co-occurrence of FGF and FGFR gene alterations')

  ## upper panel: MDS scatter plots

  rep_figs$fgfr_overlap$top <- expl_occur$mds_plots %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv') %>%
    plot_grid(get_legend(expl_occur$mds_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.83, 0.17))

  ## bottom panel: heat maps of contingency tables

  rep_figs$fgfr_overlap$bottom <- expl_occur$hm_plots %>%
    map(~.x[c('FGF3_amplification|FGF4_amplification',
             'FGF3_amplification|FGF19_amplification',
             'FGF4_amplification|FGF19_amplification',
             'FGFR3_mutation|FGF3_amplification')]) %>%
    map(map,
        ~.x +
          theme(legend.position = 'none',
                plot.title = element_blank())) %>%
    map(~plot_grid(plotlist = .,
                   ncol = 2,
                   align = 'hv',
                   axis = 'tblr'))

  rep_figs$fgfr_overlap$bottom <-
    plot_grid(rep_figs$fgfr_overlap$bottom[[1]],
              rep_figs$fgfr_overlap$bottom[[2]],
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = globals$cohort_labs[c("genie", "tcga")],
              label_size = 8)

  ## the entire figure

  rep_figs$fgfr_overlap <-
    plot_grid(rep_figs$fgfr_overlap$top,
              ggdraw(),
              rep_figs$fgfr_overlap$bottom,
              nrow = 3,
              rel_heights = c(1.1, 0.1, 1),
              labels = c('A', 'B', ''),
              label_size = 10) %>%
    as_figure(label = 'overlap_fgfr_fgf_gene_alterations',
              ref_name = 'fgfr_overlap',
              caption = paste('Co-occurrence of genetic alterations of',
                              'FGF and FGFR genes in urothelial cancer.'),
              w = 180,
              h = 210)

# Correlation of FGF and FGFR gene expression --------

  insert_msg('Correlation of FGF and FGFR gene expression')

  rep_figs$fgfr_correlation <- expl_corr$bubble_plots %>%
    map(~.x +
          guides(x = guide_axis(angle = 45)) +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              nrow = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(expl_corr$bubble_plots[[1]]),
              ncol = 2,
              rel_widths = c(0.85, 0.15)) %>%
    as_figure(label = 'fgfr_fgf_expression_correlation',
              ref_name = 'fgfr_correlation',
              caption = paste('Correlation of mRNA levels of FGF- and',
                              'FGFR-coding genes.'),
              w = 180,
              h = 220)

# Differential gene expression in tumors with and without FGFR3 mutations ------

  insert_msg('FGFR3 mutations: differential gene expression')

  ## upper panel: Volcano plots

  rep_figs$fgfr_dge$top <- expl_dge$volcano_plots$FGFR3_mut %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: box plots

  rep_figs$fgfr_dge$bottom <- expl_dge$box_plots %>%
    map(~.x[c('FGFR1', 'FGFR3', 'FGF2', 'FGF5', 'FGF7')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  rep_figs$fgfr_dge <- plot_grid(rep_figs$fgfr_dge$top,
                                 rep_figs$fgfr_dge$bottom,
                                 nrow = 2,
                                 rel_heights = c(1.2, 3),
                                 labels = LETTERS,
                                 label_size = 10) %>%
    as_figure(label = 'fgfr3_mutation_gene_expression',
              ref_name = 'fgfr_dge',
              caption = paste('Expression of FGF- and FGFR-coding genes',
                              'in urothelial cancers stratified by presence',
                              'of FGFR3 mutations.'),
              w = 180,
              h = 230)

# Total mutation burden in cancers with and without FGFR3 mutations -------

  insert_msg('Total mutation burden and FGFR3 mtations')

  rep_figs$fgfr_tmb <-
    list(tmb = map(fgfr_tmb$plots,
                   ~.x$tmb_per_mb),
         list(ggdraw()),
         mutations = map(fgfr_tmb$plots,
                         ~.x$mutations),
         list(ggdraw()),
         deletions = map(fgfr_tmb$plots[c("genie", "tcga")],
                         ~.x$deletions),
         amplifications = map(fgfr_tmb$plots[c("genie", "tcga")],
                              ~.x$amplifications)) %>%
    unlist(recursive = FALSE) %>%
    map(~.x + theme(plot.title.position = 'plot')) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', '', '',
                         'B', '', '', '',
                         'C', '', 'D', ''),
              label_size = 10) %>%
    as_figure(label = 'fgfr3_mutations_total_mutation_burden',
              ref_name = 'fgfr_tmb',
              caption = paste('Total mutation burden and counts of copy number',
                              'alterations in FGFR3 WT and',
                              'FGFR3-mutated cancers.'),
              w = 180,
              h = 210)

# Overall survival with cancers with and without FGFR3 mutations --------

  insert_msg('Overall survival in cancers with and without FGFR3 mutations')

  rep_figs$fgfr_os <-
    fgfr_os$km_plots[c("genie", "imvigor", "tcga", "tcga_wo_t4")] %>%
    map2(., c('approximate overall survival, days',
              rep('overall survival, days', 3)),
         ~.x +
           labs(x = .y) +
           theme(legend.position = 'bottom')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'fgfr3_mutation_overall_survival',
              ref_name = 'fgfr_os',
              caption = paste('Overall survival of patients with FGFR3 WT',
                              'and FGFR3-mutated cancers.'),
              w = 180,
              h = 180)

# Disease-specific and relapse-free survival with WT and FGFR3 mutated cancers ------

  insert_msg('DSS and RFS in cancers with and without FGFR3 mutations')

  rep_figs$fgfr_rfs <-
    fgfr_rfs$km_plots[c("tss", "tss_wo_t4", "rfs", "rfs_wo_t4")] %>%
    map2(., rep(c(globals$cohort_labs["tcga"],
                  paste(globals$cohort_labs["tcga"],
                        'without pT4 cancers')), 2),
         ~.x +
           labs(title = .y) +
           theme(legend.position = 'bottom')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', 'B', ''),
              label_size = 10) %>%
    as_figure(label = 'fgfr3_disease_specific_relapse_free_survival',
              ref_name = 'fgfr_rfs',
              caption = paste('Disease-specific and relapse-free survival of',
                              'patients with FGFR3 WT and FGFR3 mutated',
                              'cancers.'),
              w = 180,
              h = 180)

# Differential gene expression in cancers with and without 11q13 amplification -------

  insert_msg('Differential gene expression, 11q13 amplification')

  rep_figs$amp11q13 <-
    fgfr_ampl$plots[c("FGFR1", "FGF3", "FGF4", "FGF19")] %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 4,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = '11q13_amplification_fgf_fgfr_expression',
              ref_name = 'amp11q13',
              caption = paste('Expression of FGF- and FGFR-coding genes in',
                              'cancers with and without amplification of the',
                              '11q13 chromosome region.'),
              w = 180,
              h = 80)

# Molecular subtypes and gene alteration: heat maps --------

  insert_msg('Consensus subtypes and top gene alteration heat map')

  ## upper panel: distribution of the consensus classes

  rep_figs$fgfr_subtypes$top <-
    plot_grid(expl_sub$plot +
                coord_flip() +
                theme(legend.position = 'none',
                      axis.title.y = element_blank()),
              get_legend(expl_sub$plot +
                           theme(legend.position = 'bottom')),
              ncol = 2)

  ## bottom panel: heat maps

  rep_figs$fgfr_subtypes$bottom <- sub_genet$hm_plots %>%
    map(~.x +
          theme(legend.position = 'none',
                strip.text.y = element_blank(),
                strip.background.y = element_blank(),
                axis.text.y = element_markdown(size = 7),
                strip.text.x = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              nrow = 2,
              rel_heights = c(1, 0.5),
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(sub_genet$hm_plots[[1]]),
              ncol = 2,
              rel_widths = c(0.85, 0.15))

  ## the entire figure

  rep_figs$fgfr_subtypes <-
    plot_grid(rep_figs$fgfr_subtypes$top,
              rep_figs$fgfr_subtypes$bottom,
              nrow = 2,
              rel_heights = c(1, 3.5),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'consensus_subtypes_fgf_fgfr_gene_alterations',
              ref_name = 'fgfr_subtypes',
              caption = paste('Alterations of FGF- and FGFR-coding genes',
                              'in consensus molecular classes of urothelial',
                              'cancers.'),
              w = 180,
              h = 230)

# Molecular subtypes and gene alterations: stack plots ------

  insert_msg('Consensus molecular subtypes and FGF/FGRF gene alterations')

  ## top panel: FGFR3

  rep_figs$fgfr_subdetails$top <- sub_genet$plots %>%
    map(~.x$FGFR3_mutation) %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.text.x = element_text(size = 7),
                axis.title.x = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(sub_genet$plots[[1]]$FGFR3_mutation +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## bottom panel: FGF3/4/19 amplification

  rep_figs$fgfr_subdetails$bottom <-
    sub_genet$plots$tcga[c("FGF3_amplification",
                           "FGF4_amplification",
                           "FGF19_amplification")] %>%
    map(~.x +
          theme(legend.position = 'none',
                axis.text.x = element_text(size = 7))) %>%
    c(list(get_legend(sub_genet$plots$tcga$FGF3_amplification))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr')

  ## the complete figure

  rep_figs$fgfr_subdetails <-
    plot_grid(rep_figs$fgfr_subdetails$top,
              rep_figs$fgfr_subdetails$bottom,
              nrow = 2,
              rel_heights = c(1.1, 2),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'consensus_subtypes_fgf_fgfr_alteration_details',
              ref_name = 'fgfr_subdetails',
              caption = paste('Frequency of somatic mutations of FGFFR3',
                              'and amplification of FGF3/4/19 in consensus',
                              'molecular classes of urothelial cancer.'),
              w = 180,
              h = 210)

# Molecular subtypes and expression of FGF and FGFR genes -------

  insert_msg('Molecular subtypes and expression of FGF- and FGFR-coding genes')

  rep_figs[c('sub_fgfr_expression',
             'sub_fgf_expression1',
             'sub_fgf_expression2')] <-
    list(map(sub_dge$plots,
             ~.x[c('FGFR1', 'FGFR3')]),
         map(sub_dge$plots,
             ~.x[c('FGF2', 'FGF5', 'FGF7')]),
         map(sub_dge$plots,
             ~.x[c('FGF10', 'FGF18')])) %>%
    map(transpose) %>%
    map(unlist, recursive = FALSE) %>%
    map(map,
        ~.x + theme(legend.position = 'none',
                    axis.text.x = element_text(size = 7),
                    axis.title.x = element_blank())) %>%
    map(~plot_grid(plotlist = .x,
                   ncol = 2,
                   align = 'hv',
                   axis = 'tblr'))

  rep_figs[c('sub_fgfr_expression',
             'sub_fgf_expression1',
             'sub_fgf_expression2')] <-
    rep_figs[c('sub_fgfr_expression',
               'sub_fgf_expression1',
               'sub_fgf_expression2')] %>%
    list(x = .,
         label = c('molecular_subtypes_expression_FGFR',
                   'molecular_subtypes_expression_FGF',
                   'molecular_subtypes_expression_FGF'),
         ref_name = names(.),
         caption = paste('Differential expression of',
                         c('FGFR-coding', 'FGF-coding', 'FGF-coding'),
                         'genes in the consensus molecular classes of urothelial',
                         'cancers.'),
         h = c(120, 180, 120)) %>%
    pmap(as_figure,
         w = 180)

# Genetic subtypes: development -------

  insert_msg('Genetic subtypes: development')

  ## upper panel: log-likelihood and BIC

  rep_figs$lca_devel$top <- lca_dev$tune_object$plots[c("llik", "bic")] %>%
    map(~.x +
          geom_vline(xintercept = 6, linetype = 'dashed') +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv') %>%
    plot_grid(get_legend(lca_dev$tune_object$plots$bic),
              ncol = 2,
              rel_widths = c(0.9, 0.1))

  ## bottom panel: posterior probabilities

  rep_figs$lca_devel$bottom <- lca_pred$posterior_plots %>%
    map(~.x +
          scale_y_continuous(limits = c(0, 1)) +
          theme(legend.position = 'none',
                strip.text.x = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              nrow = 2,
              align = 'hv',
              axis = 'tblr')

  ## the entire figure

  rep_figs$lca_devel <- plot_grid(rep_figs$lca_devel$top,
                                  rep_figs$lca_devel$bottom,
                                  nrow = 2,
                                  rel_heights = c(1, 2),
                                  labels = LETTERS,
                                  label_size = 10) %>%
    as_figure(label = 'genetic_subsets_development',
              ref_name = 'lca_devel',
              caption = paste('Development of genetic subsets of urothelial',
                              'carcinoma in the GENIE BLCA training cohort.',
                              'Prediction of the genetic subset assignment',
                              'for the TCGA BLCA samples.'),
              w = 180,
              h = 210)

# Genetic subtypes: class distributions ------

  insert_msg('Genetic subtypes, class distribution')

  rep_figs$lca_size <-
    plot_grid(lca_eval$n_plot +
                theme(legend.position = 'none'),
              get_legend(lca_eval$n_plot +
                           theme(legend.position = 'bottom',
                                 legend.title = element_blank()))) %>%
    as_figure(label = 'genetic_subset_size',
              ref_name = 'lca_size',
              caption = paste('Distribution of sizes of the genetic subsets',
                              'in the training GENIE BLCA cohort and the',
                              'test TCGA BLCA collective.'),
              w = 180,
              h = 60)

# Genetic subsets: most common genetic alterations -------

  insert_msg('Genetic clusters: top genetic alterations')

  rep_figs$lca_genet <- lca_eval$hm_plots %>%
    map(~.x +
          theme(legend.position = 'none',
                strip.text = element_text(size = 7),
                axis.text.y = element_markdown(size = 7))) %>%
    plot_grid(plotlist = .,
              nrow = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(lca_eval$hm_plots[[1]]),
              ncol = 2,
              rel_widths = c(0.87, 0.13)) %>%
    as_figure(label = 'genetic_clusters_top_genetic_alterations',
              ref_name = 'lca_genet',
              caption = paste('Distribution of the most frequent somatic',
                              'mutations and copy number alterations in the',
                              'genetic subsets of urothelial cancers.'),
              w = 180,
              h = 230)

# Genetic subsets: key genetic features -------

  insert_msg('Genetic subsets, key genetic features')

  rep_figs[c('lca_details_mutations', 'lca_details_cna')] <-
    list(c('TP53_mutation',
           'RB1_mutation',
           'FGFR3_mutation',
           'ERBB2_mutation'),
         c('CCND1_amplification',
           'E2F3_amplification',
           'MDM2_amplification',
           'CDKN2A_deletion')) %>%
    map(function(ft) lca_eval$plots %>%
          map(~.x[ft]) %>%
          transpose) %>%
    map(unlist, recursive = FALSE) %>%
    map(map,
        ~.x +
          guides(x = guide_axis(angle = 45)) +
          theme(legend.position = 'none',
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 7))) %>%
    map(~plot_grid(plotlist = .x,
                   ncol = 2,
                   align = 'hv',
                   axis = 'tblr',
                   labels = c('A', '',
                              'B', '',
                              'C', '',
                              'D', ''),
                   label_size = 10))

  ## appending with the legends

  rep_figs$lca_details_mutations <-
    rep_figs$lca_details_mutations %>%
    plot_grid(get_legend(lca_eval$plots[[1]]$TP53_mutation),
              ncol = 2,
              rel_widths = c(0.85, 0.15))

  rep_figs$lca_details_cna <-
    rep_figs$lca_details_cna %>%
    plot_grid(.,
              plot_grid(ggdraw(),
                        get_legend(lca_eval$plots[[1]]$MDM2_amplification),
                        ggdraw(),
                        get_legend(lca_eval$plots[[1]]$CDKN2A_deletion),
                        nrow = 4),
              ncol = 2,
              rel_widths = c(0.85, 0.15))

  ## entire figures

  rep_figs[c('lca_details_mutations', 'lca_details_cna')] <-
    rep_figs[c('lca_details_mutations', 'lca_details_cna')] %>%
    list(x = .,
         label = c('genetic_subsets_somatic_mutation_details',
                   'genetic_subsets_copy_number_variant_details'),
         ref_name = names(.),
         caption = paste('Frequency of selected',
                         c('somatic mutations', 'copy number alterations'),
                         'in the genetic subsets of urothelial cancers.')) %>%
    pmap(as_figure,
         w = 180,
         h = 230)

# Genetic subsets: mutation burdens -------

  insert_msg('Genetic subset: mutation burden')

  rep_figs$lca_tmb <- lca_tmb$plots %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map2(., c('identity', rep('sqrt', 7)),
         ~.x +
           scale_y_continuous(trans = .y) +
           guides(x = guide_axis(angle = 45)) +
           theme(legend.position = 'none',
                 axis.text.x = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', 'B', '', 'C', '', 'D', ''),
              label_size = 10) %>%
    as_figure(label = 'genetic_subsets_mutation_burden',
              ref_name = 'lca_tmb',
              caption = paste('Total mutation burden, and counts of somatic',
                              'mutations and copy number alterations in the',
                              'genetic subsets of urothelial cancers.'),
              w = 180,
              h = 230)

# Genetic subsets: differential gene expression -------

  insert_msg('Genetic subsets: differential gene expression')

  rep_figs[c('lca_dge_fgfr', 'lca_dge_fgf')] <-
    list(lca_dge$plots[c("FGFR1", "FGFR2", "FGFR3")],
         lca_dge$plots[c('FGF2', 'FGF5', 'FGF7', 'FGF17')]) %>%
    map(map,
        ~.x +
          guides(x = guide_axis(angle = 45)) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(size = 7),
                legend.position = 'none')) %>%
    map(~plot_grid(plotlist = .,
                   ncol = 2,
                   align = 'hv',
                   axis = 'tblr'))

  rep_figs[c('lca_dge_fgfr', 'lca_dge_fgf')] <-
    rep_figs[c('lca_dge_fgfr', 'lca_dge_fgf')] %>%
    list(x = .,
         label = c('genetic_subsets_fgfr_receptor_expression',
                   'genetic_subsets_fgf_ligand_expression'),
         ref_name = names(.),
         caption = paste('Differential expression of',
                         c('FGFR-coding', 'FGF-coding'),
                         'genes in the genetic subsets of',
                         'urothelial cancers.')) %>%
    pmap(as_figure,
         w = 180,
         h = 140)

# Genetic subsets: clinical features --------

  insert_msg('Genetic subsets: clinical features')

  rep_figs$lca_clinic <- lca_clinic$plots %>%
    map(~.x[c('age', 'sex')]) %>%
    transpose %>%
    unlist(recursive = FALSE) %>%
    map(~.x +
          guides(x = guide_axis(angle = 45)) +
          theme(legend.position = 'none',
                axis.title.x = element_blank(),
                axis.text.x = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv') %>%
    plot_grid(plot_grid(ggdraw(),
                        get_legend(lca_clinic$plots[[1]]$sex),
                        nrow = 2),
              ncol = 2,
              rel_widths = c(0.85, 0.15)) %>%
    as_figure(label = 'genetic_subsets_age_gender',
              ref_name = 'lca_clinic',
              caption = 'Age and gender in the genetic subsets.',
              w = 180,
              h = 140)

# Genetic subsets: overall survival ---------

  insert_msg('Genetic subsets: overall survival')

  rep_figs$lca_os <- lca_os$km_plots %>%
    map2(., c('approximate overall survival, days',
             rep('overall survival, days', 2)),
        ~.x +
          labs(x = .y) +
          theme(legend.position = 'right',
                legend.key.height = unit(11, 'mm'))) %>%
    plot_grid(plotlist = .,
              nrow = 3) %>%
    as_figure(label = 'genetic_subsets_overall_survival',
              ref_name = 'lca_os',
              caption = paste('Overall survival in the genetic subsets',
                              'of urothelial cancers.'),
              w = 140,
              h = 230)

# Saving on the disc ------

  insert_msg('Saving figures on the disc')

  rep_figs %>%
    number_figures %>%
    walk(pickle,
         path = './report/figures',
         format = 'pdf',
         device = cairo_pdf)

# END ----

  insert_tail()
