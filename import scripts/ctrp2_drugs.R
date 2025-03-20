# Import, wrangling and basic exploratory analysis of durg sensitivity data
# of the CTRP2 project.

  insert_head()

# container --------

  ctrp2_drugs <- list()

# drug sensitivity: AUC -------

  insert_msg('Drug sensitivity, AUC')

  ctrp2_drugs$auc <-
    read_csv('./data/CTRP2/Drug_sensitivity_AUC_(CTD^2).csv')

  names(ctrp2_drugs$auc)[1] <- 'sample_id'

# drug lexicon -------

  insert_msg('Drug lexicon')

  ctrp2_drugs$lexicon <-
    tibble(variable = names(ctrp2_drugs$auc)[-1]) %>%
    mutate(drug_id = stri_extract(variable,
                                  regex = '\\(CTRP.*\\)$'),
           drug_id = stri_replace_all(drug_id,
                                      regex = '\\(|\\)',
                                      replacement = ''),
           drug_id = stri_replace(drug_id,
                                  regex = '^CTRP:',
                                  replacement = ''))

  ## appending the lexicon with drug targets and MOA

  ctrp2_drugs$appendix <- read_tsv('./data/CTRP2/v20.meta.per_compound.txt')

  ctrp2_drugs$appendix <- ctrp2_drugs$appendix %>%
    transmute(drug_id = as.character(master_cpd_id),
              drug_name = cpd_name,
              label = drug_name,
              targets = stri_split_fixed(gene_symbol_of_protein_target,
                                         pattern = ';'),
              moa = target_or_activity_of_compound,
              broad_cpd_id = broad_cpd_id)

  ctrp2_drugs$lexicon <-
    left_join(ctrp2_drugs$lexicon,
              ctrp2_drugs$appendix,
              by = 'drug_id')

# cell line details --------

  insert_msg('Cell line details')

  ctrp2_drugs$samples <-
    read_csv('./data/CTRP2/Model.csv')

  ctrp2_drugs$samples <- ctrp2_drugs$samples %>%
    transmute(sample_id = ModelID,
              cell_line = StrippedCellLineName,
              tcga_entity = OncotreeCode,
              tissue = OncotreeLineage,
              tissue_subtype = OncotreeSubtype) %>%
    filter(sample_id %in% ctrp2_drugs$auc$sample_id)

# Exploratory analysis: counts of tissue types, entities and subtypes -------

  insert_msg('Exploratory analysis: counts')

  ctrp2_drugs[c('entity_counts',
                'tissue_counts',
                'tissue_subtype_counts')] <-
    c('tcga_entity',
      'tissue',
      'tissue_subtype') %>%
    map(~count(ctrp2_drugs$samples, .data[[.x]]))

# Caching the results -------

  insert_msg('Caching the results')

  ctrp2_drugs$appendix <- NULL

  ctrp2_drugs <- compact(ctrp2_drugs)

  save(ctrp2_drugs, file = './data/CTRP2/ctrp2_drugs.RData')

# END -------

  insert_tail()
