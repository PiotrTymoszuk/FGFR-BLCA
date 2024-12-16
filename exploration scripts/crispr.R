# Exploratory analysis for the results of CHRONOS CRISPR screening for the
# FGFR1-, FGF-, and FGFBP-coding genes of interest.
#
# Cell lines with significant gene effects are defined by gene effects
# beyond the 95th percentile of gene effects of non-expressed genes.

  insert_head()

# container -------

  expl_crispr <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## effects of genes of interest

  expl_crispr$data <- depmap$crispr %>%
    filter(complete.cases(.))

  expl_crispr$variables <- names(expl_crispr$data)[-1]

  ## background effects of non-expressed genes,
  ## and their 95% percentile range

  expl_crispr$background <-
    tibble(background = as.numeric(as.matrix(depmap$crispr_bcg[-1]))) %>%
    filter(complete.cases(.))

  expl_crispr$background_ci <-
    quantile(expl_crispr$background$background, c(0.025, 0.975))

  ## identification of significant gene effects

  expl_crispr$effects <- expl_crispr$data

  expl_crispr$effects[-1] <- expl_crispr$effects[-1] %>%
    map_dfc(cut,
            c(-Inf, expl_crispr$background_ci, Inf),
            c('negative', 'ns', 'positive')) %>%
    as.data.frame %>%
    set_rownames(rownames(expl_crispr$data))

  ## FGFR3 mutation status

  expl_crispr$FGFR3_mutation <-
    depmap$mutation[, c("sample_id", "FGFR3_mut")] %>%
    mutate(FGFR3_mutation = fct_recode(FGFR3_mut,
                                       WT = 'no',
                                       mutated = 'yes'))

  ## plotting data in a long format used for plotting,
  ## appended with cell line names for the significant effects,
  ## the background gene effects, and FGFR3 mutation status

  expl_crispr$long_data <-
    map2(expl_crispr[c("data", "effects")],
         c('effect', 'effect_class'),
         ~pivot_longer(.x,
                       cols = all_of(expl_crispr$variables),
                       names_to = 'variable',
                       values_to = .y)) %>%
    reduce(left_join, by = c('variable', 'sample_id')) %>%
    mutate(name_lab = ifelse(effect_class != 'ns',
                             exchange(sample_id,
                                      depmap$cell_lines,
                                      key = 'sample_id',
                                      value = 'cell_line_name'),
                             NA),
           gene_symbol = stri_replace(variable,
                                      regex = '_crispr$',
                                      replacement = ''),
           axis_lab = html_italic(gene_symbol),
           protein_type = ifelse(gene_symbol %in% globals$receptors,
                                 'receptor',
                                 ifelse(gene_symbol %in% globals$ligands,
                                        'ligand', 'binding_protein')),
           protein_type = factor(protein_type,
                                 c('receptor', 'ligand', 'binding_protein')))

  expl_crispr$long_data <-
    rbind(expl_crispr$long_data,
          tibble(sample_id = 'background',
                 variable = 'background',
                 effect = expl_crispr$background$background,
                 effect_class = factor('ns', levels(expl_crispr$effects[[1]])),
                 name_lab = NA,
                 gene_symbol = 'background',
                 axis_lab = 'background',
                 protein_type = 'background'))

  expl_crispr$long_data <-
    left_join(expl_crispr$long_data,
              expl_crispr$FGFR3_mutation,
              by = 'sample_id')

# Box plots of gene effects --------

  insert_msg('Box plots of gene effects')

  ## separate box plot panels for the receptors,
  ## ligands, and binding protein genes.
  ## Setting factor levels to make the background appear always at the
  ## bottom of the plot

  expl_crispr$box_data <-
    list(receptors = c('receptor', 'background'),
         ligands = c('ligand', 'background'),
         binding_proteins = c('binding_protein', 'background')) %>%
    map(~filter(expl_crispr$long_data, protein_type %in% .x)) %>%
    map2(., globals[c("receptors", "ligands", "binding_proteins")],
         ~mutate(.x,
                 axis_lab = factor(axis_lab,
                                   levels = c('background',
                                              rev(html_italic(.y))))))

  expl_crispr$box_plots <-
    list(x = expl_crispr$box_data,
         y = c('FGFR', 'FGF', 'FGFBP') %>%
           paste0('-coding genes')) %>%
    pmap(function(x, y) x %>%
           ggplot(aes(x = effect,
                      y = axis_lab,
                      fill = protein_type,
                      alpha = protein_type,
                      color = effect_class)) +
           facet_grid(protein_type ~ .,
                      scales = 'free',
                      space = 'free',
                      labeller = as_labeller(function(x) stri_replace(x,
                                                                      fixed = '_',
                                                                      replacement = '\n'))) +
           geom_vline(xintercept = expl_crispr$background_ci[[1]],
                      linetype = 'dashed') +
           geom_vline(xintercept = expl_crispr$background_ci[[2]],
                      linetype = 'dashed') +
           geom_boxplot(outlier.color = NA,
                        alpha = 0.25,
                        color = 'black',
                        linewidth = 0.5) +
           geom_point(aes(shape = FGFR3_mutation),
                      size = 2,
                      position = position_jitter(0, 0.1)) +
           geom_text_repel(aes(label = name_lab),
                           size = 2.3,
                           show.legend = FALSE) +
           scale_shape_manual(values = c(WT = 16,
                                         mutated = 4),
                              name = html_italic('FGFR3')) +
           scale_alpha_manual(values = c(background = 0.0,
                                         receptor = 1,
                                         ligand = 1,
                                         binding_protein = 1)) +
           scale_fill_manual(values = globals$protein_colors) +
           scale_color_manual(values = c(negative = 'steelblue',
                                         ns = 'gray70',
                                         positive = 'firebrick'),
                              name = 'Significant<br>gene effect') +
           scale_x_continuous(limits = range(expl_crispr$long_data$effect)) +
           guides(fill = 'none',
                  alpha = 'none') +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_markdown(),
                 strip.background = element_blank(),
                 strip.text = element_blank(),
                 legend.title = element_markdown()) +
           labs(title = y,
                subtitle = 'CHRONOS, CRISPR screening',
                x = 'gene effect, CHRONOS score'))

# END -------

  expl_crispr <-
    expl_crispr[c("background_ci", "effects", "box_plots")]

  insert_tail()
