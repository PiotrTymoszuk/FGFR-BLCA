# Import and clearing of drug response data from the PRISM project.
# The response metric is AUC, i.e. area under the dose-response curve, provided
# by the study authors.
# The data was obtained from the depMap portal:
# https://depmap.org/portal/download/all/?releasename=CTRP+CTD%5E2&filename=CTRPv2.0_2015_ctd2_ExpandedDataset.zip

  insert_head()

# container -------

  prism_drugs <- list()

# drug sensitivity: AUC and drug lexicon -------

  insert_msg('Drug sensitivity, AUC and drug lexicon')

  ## raw AUCs

  prism_drugs$auc <-
    read_csv('./data/PRISM/secondary-screen-dose-response-curve-parameters.csv')

  ## drug lexicon skeleton

  prism_drugs$lexicon <- prism_drugs$auc %>%
    transmute(variable = make.names(name),
              drug_name = name,
              label = name,
              moa = moa,
              targets = stri_split_fixed(target, pattern = ', '))

# Selection of reliable drug AUC, wrangling -------

  insert_msg('Selection of reliable AUC, wide format')

  ## a liberal selection of reliable drug response curves with r2 >= 0.26
  ## averaging AUC in the screens

  prism_drugs$auc <- prism_drugs$auc %>%
    filter(r2 >= 0.26) %>%
    mutate(sample_id = depmap_id) %>%
    select(sample_id, screen_id, name, auc) %>%
    group_by(sample_id, name) %>%
    summarise(auc = mean(auc, na.rm = TRUE)) %>%
    ungroup

  ## conversion to a wide format: drug names in the columns,
  ## cell lines in the rows

  prism_drugs$auc <- prism_drugs$auc %>%
    mutate(name = make.names(name)) %>%
    pivot_wider(id_cols = 'sample_id',
                names_from = 'name',
                values_from = 'auc')

  ## update of the lexicon

  prism_drugs$lexicon <- prism_drugs$lexicon %>%
    filter(!duplicated(variable)) %>%
    filter(variable %in% names(prism_drugs$auc))

# cell line details --------

  insert_msg('Cell line details')

  prism_drugs$samples <-
    read_csv('./data/PRISM/Model.csv')

  prism_drugs$samples <- prism_drugs$samples %>%
    transmute(sample_id = ModelID,
              cell_line = StrippedCellLineName,
              tcga_entity = OncotreeCode,
              tissue = OncotreeLineage,
              tissue_subtype = OncotreeSubtype) %>%
    filter(sample_id %in% prism_drugs$auc$sample_id)

  prism_drugs$auc <- prism_drugs$auc %>%
    filter(sample_id %in% prism_drugs$samples$sample_id)

# Exploratory analysis: counts of tissue types, entities and subtypes -------

  insert_msg('Exploratory analysis: counts')

  prism_drugs[c('entity_counts',
                'tissue_counts',
                'tissue_subtype_counts')] <-
    c('tcga_entity',
      'tissue',
      'tissue_subtype') %>%
    map(~count(prism_drugs$samples, .data[[.x]]))

# Caching the results -------

  insert_msg('Caching the results')

  prism_drugs$appendix <- NULL

  prism_drugs <- compact(prism_drugs)

  save(prism_drugs, file = './data/PRISM/prism_drugs.RData')

# END -------

  insert_tail()
