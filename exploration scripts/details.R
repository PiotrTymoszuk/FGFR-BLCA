# Details of mutation type and mutation location for FGFR1 - 4.
# The analysis is done for the GENIE and TCGA cohorts.

  insert_head()

# container -----

  expl_fgfr <- list()

# analysis data -------

  insert_msg('Analysis data')

  ## detailed data and numbers of complete cases

  expl_fgfr$data <- list(genie = genie,
                         msk = msk,
                         tcga = tcga) %>%
    map(~.x$mutation_detail)

  expl_fgfr$n_totals <- list(genie = genie,
                             msk = msk,
                             tcga = tcga) %>%
    map(~.x$mutation) %>%
    map_dbl(nrow)

  ## some minimal wrangling for common mutation classification
  ## and consistent protein positions

  expl_fgfr$data$genie <- expl_fgfr$data$genie %>%
    transmute(sample_id = sampleId,
              gene_symbol = hugoGeneSymbol,
              dna_start_position = as.numeric(startPosition),
              dna_end_position = as.numeric(endPosition),
              dna_reference = referenceAllele,
              dna_variant = variantAllele,
              protein_change = proteinChange,
              protein_start_position = proteinPosStart,
              protein_end_position = proteinPosEnd,
              mutation_type = car::recode(mutationType,
                                          "'Frame_Shift_Del' = 'frame shift deletion';
                                          'Frame_Shift_Ins' = 'frame shift insertion';
                                          'In_Frame_Del' = 'in-frame deletion';
                                          'In_Frame_Ins' = 'in-frame insertion';
                                          'Missense_Mutation' = 'missense';
                                          'Nonsense_Mutation' = 'nonsense';
                                          'Nonstop_Mutation' = 'stop loss';
                                          'Splice_Region' = 'splice region';
                                          'Splice_Site' = 'splice site'"),
              variant_type = car::recode(variantType,
                                         "'DEL' = 'deletion';
                                         'INS' = 'insertion'"))

  expl_fgfr$data[c('msk', 'tcga')] <- expl_fgfr$data[c('msk', 'tcga')] %>%
    map(transmute,
        sample_id = Tumor_Sample_Barcode,
        gene_symbol = Hugo_Symbol,
        dna_start_position = Start_Position,
        dna_end_position = End_Position,
        dna_reference = Reference_Allele,
        dna_variant = Tumor_Seq_Allele2,
        protein_change = stri_replace(HGVSp_Short,
                                      regex = '^p\\.',
                                      replacement = ''),
        protein_start_position = Protein_position,
        mutation_type = car::recode(Variant_Classification ,
                                    "'Frame_Shift_Del' = 'frame shift deletion';
                                    'Frame_Shift_Ins' = 'frame shift insertion';
                                    'In_Frame_Ins' = 'in-frame insertion';
                                    'Intron' = 'intron';
                                    'Missense_Mutation' = 'missense';
                                    'Nonsense_Mutation' = 'nonsense';
                                    'Silent' = 'silent';
                                    'Splice_Region' = 'splice region';
                                    'Splice_Site' = 'splice site'"),
        variant_type = car::recode(Variant_Type,
                                   "'DEL' = 'deletion'; 'INS' = 'insertion'"))

  ## common factor levels

  expl_fgfr$data <- expl_fgfr$data %>%
    map(mutate,
        mutation_type = factor(mutation_type,
                               c('missense',
                                 'nonsense',
                                 'stop loss',
                                 'in-frame deletion',
                                 'in-frame insertion',
                                 'frame shift deletion',
                                 'frame shift insertion',
                                 'splice region',
                                 'splice site',
                                 'intron',
                                 "3'UTR")),
        variant_type = factor(variant_type,
                              c('SNP', 'DNP', 'ONP',
                                'deletion', 'insertion')))

# Domain data -------

  insert_msg('Domain data')

  expl_fgfr$domain_lst <- paste0('FGFR', 1:4) %>%
    paste0('./data/domains/', ., '.json') %>%
    map(read_json, simplifyVector = TRUE) %>%
    map(~.x$features) %>%
    map(filter, type == 'Domain') %>%
    map(~cbind(type = .x[['type']],
              .x[['location']][['start']],
              .x[['location']][['end']],
              desciption = .x[['description']])) %>%
    map(set_names,
        c('type', 'start_position', 'start_type',
          'end_position', 'end_type', 'description')) %>%
    set_names(paste0('FGFR', 1:4)) %>%
    map(as_tibble)

  ## appending the mutation data

  expl_fgfr$data <- expl_fgfr$data %>%
    map(append_domain,
        domain_lst = expl_fgfr$domain_lst) %>%
    map(mutate,
        protein_domain = factor(protein_domain,
                                c('Ig-like C2-type 1',
                                  'Ig-like C2-type 2',
                                  'Ig-like C2-type 3',
                                  'Protein kinase',
                                  'not assigned')))

# tallying the mutation types and domains -------

  insert_msg('Frequencies of mutation and variant types, and domains')

  for(i in c('mutation_type', 'variant_type', 'protein_domain')) {

    expl_fgfr$stats[[i]] <- expl_fgfr$data %>%
      map(blast, gene_symbol) %>%
      map(map, count, .data[[i]]) %>%
      map(map, arrange, desc(.data[[i]])) %>%
      map(map, mutate, plot_pos = cumsum(n) - 0.5 * n) %>%
      map(compress, names_to = 'gene_symbol') %>%
      map2(., expl_fgfr$n_totals,
           ~mutate(.x,
                   n_total = .y,
                   percent = n/n_total * 100,
                   plot_pos = plot_pos/n_total * 100))

  }

# Stack plots with mutation and variant types -------

  insert_msg('Stack plots with mutation and variant types')

  expl_fgfr$stack_plots <-
    list(data = unlist(expl_fgfr$stats, recursive = FALSE),
         fill_var = c(rep('mutation_type', 3),
                      rep('variant_type', 3),
                      rep('protein_domain', 3)),
         plot_title = c(paste('Mutation type,',
                              globals$cohort_labs[c("genie", "msk", "tcga")]),
                        paste('Variant type,',
                              globals$cohort_labs[c("genie", "msk", "tcga")]),
                        paste('Protein domain,',
                              globals$cohort_labs[c("genie", "msk", "tcga")])),
         plot_subtitle = paste('total samples: n =',
                               rep(expl_fgfr$n_totals, 3)),
         fill_title = '',
         palette = globals[c(rep("mutation_colors", 3),
                             rep("variant_colors", 3),
                             rep("domain_colors", 3))],
         reverse_fill = c(rep(TRUE, 6), rep(FALSE, 3))) %>%
    pmap(plot_stack,
         x_var = 'percent',
         y_var = 'gene_symbol',
         y_limits = paste0('FGFR', c(1, 4, 2, 3))) %>%
    map(~.x +
          theme(axis.text.y = element_text(face = 'italic')))

# Position plots ---------

  insert_msg('Position plots')

  expl_fgfr$position_data <- expl_fgfr$data %>%
    map(blast, gene_symbol)

  for(i in names(expl_fgfr$position_data)) {

    expl_fgfr$position_plots[[i]] <-
      list(data = expl_fgfr$position_data[[i]],
           plot_title = expl_fgfr$position_data[[i]] %>%
             names %>%
             html_italic %>%
             paste(globals$cohort_labs[[i]], sep = ', '),
           protein_length = globals$receptor_lengths) %>%
      pmap(position_plot,
           fill_var = 'protein_domain',
           palette = globals$domain_colors) %>%
      map(~.x +
            theme(plot.title = element_markdown()) +
            labs(fill = 'Protein domain'))

  }

# Domain schemes -------

  insert_msg('Domain schemes')

  expl_fgfr$domain_plots <-
    list(domain_tbl = expl_fgfr$domain_lst,
         protein_name = names(expl_fgfr$domain_lst),
         protein_length = globals$receptor_lengths) %>%
    pmap(plot_domains)

  ## adding the domain schemes to the position plots

  expl_fgfr$position_plots$genie <-
    list(position_plot = expl_fgfr$position_plots$genie,
         domain_tbl = expl_fgfr$domain_lst,
         protein_name = names(expl_fgfr$domain_lst),
         protein_length = globals$receptor_lengths,
         rel_heights = list(c(2, 4 + 0.5),
                            c(2, 4 + 0.5),
                            c(2, 9 + 0.5),
                            c(2, 3 + 0.5))) %>%
    pmap(append_domain_scheme,
         align = 'v',
         rm_position_legend = TRUE,
         rm_domain_legend = TRUE)

  expl_fgfr$position_plots$msk <-
    list(position_plot = expl_fgfr$position_plots$msk,
         domain_tbl = expl_fgfr$domain_lst,
         protein_name = names(expl_fgfr$domain_lst),
         protein_length = globals$receptor_lengths,
         rel_heights = list(c(2, 3 + 0.5),
                            c(2, 2 + 0.5),
                            c(2, 6 + 0.5),
                            c(2, 1 + 0.5))) %>%
    pmap(append_domain_scheme,
         alig = 'v',
         rm_position_legend = TRUE,
         rm_domain_legend = TRUE)

  expl_fgfr$position_plots$tcga <-
    list(position_plot = expl_fgfr$position_plots$tcga,
         domain_tbl = expl_fgfr$domain_lst,
         protein_name = names(expl_fgfr$domain_lst),
         protein_length = globals$receptor_lengths,
         rel_heights = list(c(2, 1 + 0.5),
                            c(2, 2 + 0.5),
                            c(2, 5 + 0.5),
                            c(2, 1 + 0.5))) %>%
    pmap(append_domain_scheme,
         alig = 'v',
         rm_position_legend = TRUE,
         rm_domain_legend = TRUE)

# Result tables --------

  insert_msg('Result tables')

  ## frequency of mutations split by type, variant type and protein domain

  expl_fgfr$result_tbl <- expl_fgfr$stats %>%
    map(compress, names_to = 'cohort') %>%
    map(mutate, cohort = globals$cohort_labs[cohort]) %>%
    map2(., names(.),
         ~select(.x, cohort, gene_symbol, all_of(.y), percent, n, n_total)) %>%
    map2(., names(.),
         ~arrange(.x, cohort, gene_symbol, .data[[.y]])) %>%
    map2(., c('Mutation type',
              'Variant type',
              'Protein domain'),
         ~set_names(.x,
                    c('Cohort', 'Gene symbol', .y,
                      'Percentage of samples',
                      'Samples with mutation, N',
                      'Total samples, N')))

# END ------

  rm(i)

  expl_fgfr$position_data <- NULL

  expl_fgfr <- compact(expl_fgfr)

  insert_tail()
