# Project globals: graphics, variable names, labels and units

  insert_head()

# container ------

  globals <- list()

# graphic globals ------

  globals$common_text <- element_text(size = 8,
                                      face = 'plain',
                                      color = 'black')

  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = 'mm')

  ## common theme for plots

  globals$common_theme <- theme_classic() +
    theme(axis.text = globals$common_text,
          axis.title = globals$common_text,
          plot.title = element_text(size = 8,
                                    face = 'bold'),
          plot.subtitle = globals$common_text,
          plot.tag = element_text(size = 8,
                                  face = 'plain',
                                  color = 'black',
                                  hjust = 0,
                                  vjust = 1),
          plot.tag.position = 'bottom',
          legend.text = globals$common_text,
          legend.title = globals$common_text,
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'gray95',
                                          color = 'gray80'),
          plot.margin = globals$common_margin,
          panel.grid.major = element_line(color = 'gray90'))

  ## common theme for networks

  globals$net_theme <- theme_void() +
    theme(plot.title = element_text(size = 8,
                                    face = 'bold'),
          plot.subtitle = globals$common_text,
          legend.text = globals$common_text,
          legend.title = globals$common_text,
          plot.margin = globals$common_margin)

# cohort labs and colors ------

  insert_msg('Color labs and colors')

  globals$cohort_labs <- c(genie = 'GENIE BLCA',
                           msk = 'MSK IMPACT',
                           tcga = 'TCGA BLCA',
                           imvigor = 'IMvigor',
                           bcan = 'BCAN',
                           hpa = 'HPA')

  globals$cohort_colors <- c(genie = 'indianred3',
                             msk = 'darkolivegreen4',
                             tcga = 'steelblue',
                             imvigor = 'gray60',
                             bcan = 'plum4')

  ## a lexicon data frame for convenient annotation with `exchange()`

  globals$cohort_lexicon <-
    globals[c("cohort_labs", "cohort_colors")] %>%
    map2(., c('label', 'color'),
         ~compress(.x, names_to = 'cohort', values_to = .y)) %>%
    reduce(left_join, by = 'cohort')

  ## a ready-to-use expression for a list with the study data set
  ## this expression includes only the studies passing the gene number
  ## criterion

  globals$analysis_cohorts <-
    globals$cohort_labs[names(globals$cohort_labs) != 'hpa'] %>%
    names

  globals$cohort_expr <- globals$analysis_cohorts %>%
    map_chr(~paste(.x, .x, sep = ' = ')) %>%
    paste(collapse = ', ') %>%
    paste0('list(', ., ')') %>%
    parse_expr

  ## common units

  globals$cohort_axis_labs <-
    c(tcga = 'log<sub>2</sub> gene count',
      imvigor = 'log<sub>2</sub> gene count',
      bcan = 'Z-score',
      hpa = 'IHC score')

# Genes of interest -------

  insert_msg('Genes of interest')

  ## gene symbols

  globals$receptors <- c(paste0('FGFR', 1:4),
                         'FGFRL1', 'TGFBR3')

  globals$ligands <- paste0('FGF', c(1:14, 17:23))

  globals$binding_proteins <- c(paste0('FGFBP', 1:3),
                                'KL', 'KLB', 'PTX3',
                                'SDC1', 'SDC2', 'SDC4',
                                'DCN', 'CD44', 'HSPG2',
                                'GPC1', 'GPC3', 'GPC4',
                                'TNFAIP6', 'FIBP', 'ANOS1')

  ## protein lengths

  globals$receptor_lengths <-
    c(FGFR1 = 822,
      FGFR2 = 821,
      FGFR3 = 806,
      FGFR4 = 802)

# Consensus subtype colors -------

  globals$sub_colors <-
    c('LumP' = 'steelblue',
      'LumU' =  'goldenrod',
      'LumNS' = 'gray40',
      'Stroma-rich' = 'orangered3',
      'Ba/Sq' = 'darkolivegreen4',
      'NE-like' = 'plum4')

  globals$sub_labels <-
    c('LumP' = 'luminal papillary',
      'LumU' = 'luminal unstable',
      'LumNS' = 'luminal non-specified',
      'Stroma-rich' = 'stroma-rich',
      'Ba/Sq' = 'basal/squamous',
      'NE-like' = 'neuroendocrine-like',
      'not assigned' = 'not assigned')

# clinical variables ------

  insert_msg('Clinical variables')

  globals$clinic_lexicon <-
    tibble(variable = c('age', 'sex',
                        'tissue', 'invasiveness',
                        'pt_stage', 'pn_stage', 'pm_stage',
                        'tmb_per_mb',
                        'bi_response', 'death'),
           label = c('Age', 'Gender',
                     'Tissue', 'Invasiveness',
                     'pT stage', 'pN stage', 'pM stage',
                     'Mutation burden',
                     'Best overall response', 'Mortality')) %>%
    mutate(format = ifelse(variable %in% c('age', 'tmb_per_mb'),
                           'numeric', 'factor'),
           unit = ifelse(variable == 'age', 'years',
                         ifelse(variable == 'tmb_per_mb',
                                'mutations/MB', NA)))

# Mutation and variant types, and domains -------

  globals$mutation_colors <-
    c('missense' = 'indianred3',
      'nonsense' = 'aquamarine3',
      'stop loss' = 'brown2',
      'in-frame deletion' = 'steelblue2',
      'in-frame insertion' = 'firebrick2',
      'frame shift deletion' = 'steelblue4',
      'frame shift insertion' = 'firebrick4',
      'splice region' = 'cornsilk',
      'splice site' = 'bisque3',
      'intron' = 'darkolivegreen',
      "3'UTR" = 'orangered',
      'silent' = 'gray70')

  globals$variant_colors <-
    c('SNP' = 'orangered3',
      'DNP' = 'bisque2',
      'ONP' = 'bisque4',
      'deletion' = 'steelblue',
      'insertion' = 'firebrick')

  globals$domain_colors <-
    c('Ig-like C2-type 1' = 'steelblue1',
      'Ig-like C2-type 2' = 'steelblue3',
      'Ig-like C2-type 3' = 'steelblue4',
      'Protein kinase' = 'orangered3',
      'not assigned' = 'gray70')

# mutation, amplification and deletion status colors -------

  globals$mut_status_color <- c('WT' = 'gray70',
                                'mutated' = 'plum4')

  globals$del_status_color <- c('WT' = 'gray70',
                                'deleted' = 'steelblue')

  globals$amp_status_color <- c('WT' = 'gray70',
                                'amplified' = 'firebrick')

# genetic cluster colors ------

  globals$genet_colors <-
    c('mutRB1' = 'bisque3',
      'hyperMut' = 'plum4',
      'mutFGFR3' = 'aquamarine3',
      'oligoMut' = 'orangered3',
      'ampMDM2' = 'firebrick',
      'del9p21' = 'steelblue')

# END -------

  insert_tail()
