# Globals used by multiple scripts

  insert_head()

# container -------

  expl_globals <- list()

# Details of FGFR gene mutations --------

  insert_msg('Detailed mutation data')

  ## detailed data and numbers of complete cases

  expl_globals$data <- globals$cohort_expr %>%
    eval %>%
    map(~.x$mutation_detail) %>%
    compact

  expl_globals$n_totals <- globals$cohort_expr %>%
    eval %>%
    map(~.x$mutation) %>%
    map_dbl(nrow)

  expl_globals$n_totals <- expl_globals$n_totals[names(expl_globals$data)]

  ## some minimal wrangling for common mutation classification
  ## and consistent protein positions

  expl_globals$data$genie <- expl_globals$data$genie %>%
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

  expl_globals$data[c('msk', 'tcga', 'bcan')] <-
    expl_globals$data[c('msk', 'tcga', 'bcan')] %>%
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
        mutation_type = ifelse(mutation_type == "3'Flank",
                               "3'UTR", mutation_type),
        variant_type = car::recode(Variant_Type,
                                   "'DEL' = 'deletion'; 'INS' = 'insertion'"))

  ## common factor levels

  expl_globals$data <- expl_globals$data %>%
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
                                 "3'UTR",
                                 'silent')),
        variant_type = factor(variant_type,
                              c('SNP', 'DNP', 'ONP',
                                'deletion', 'insertion')))

  # Domain data -------

  insert_msg('Domain data')

  expl_globals$domain_lst <- paste0('FGFR', 1:4) %>%
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

  expl_globals$data <- expl_globals$data %>%
    map(append_domain,
        domain_lst = expl_globals$domain_lst) %>%
    map(mutate,
        protein_domain = factor(protein_domain,
                                c('Ig-like C2-type 1',
                                  'Ig-like C2-type 2',
                                  'Ig-like C2-type 3',
                                  'Protein kinase',
                                  'not assigned')))

# END -------

  insert_tail()
