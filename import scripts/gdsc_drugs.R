# Wrangling and basic exploratory analysis for GDSC1 and GDSC2 drug
# response data sets.

  insert_head()

# container -------

  gdsc_drugs <- list()

# drug lexicons --------

  insert_msg('Drug lexicons')

  ## variable naming scheme:
  ## compound name_identifier, drug targets stored as a list

  gdsc_drugs$lexicons <-
    c(gdsc1 = './data/GDSC/compounds_gdsc1.tab',
      gdsc2 = './data/GDSC/compounds_gdsc2.tab') %>%
    map(read_tsv) %>%
    map(mutate,
        variable = paste(drug_name, drug_id, sep = '_'),
        targets = ifelse(targets == '-', NA, targets),
        pathway_name = ifelse(pathway_name == '-', NA, pathway_name),
        targets = stri_split_fixed(targets,
                                   pattern = ', ',
                                   simplify = FALSE))


# IC50 data --------

  insert_msg('Experiment design and IC50')

  gdsc_drugs$design <- c('gdsc1' = './data/GDSC/GDSC1_IC50.csv',
                         'gdsc2' = './data/GDSC/GDSC2_IC50.csv') %>%
    map(read_csv)

  ## clearing

  gdsc_drugs$design <- gdsc_drugs$design %>%
    map(transmute,
        variable = paste(`Drug Name`, `Drug ID`, sep = '_'),
        drug_id = `Drug ID`,
        drug_name = `Drug Name`,
        sample_id = paste0('cosmic_', `Cosmic ID`),
        cosmic_id = `Cosmic ID`,
        cell_line = `Cell Line Name`,
        tcga_entity = `TCGA Classification`,
        tissue = Tissue,
        tissue_subtype = `Tissue Sub-type`,
        log_ic50 = IC50,
        auc = AUC,
        rmse = RMSE)

# Wide-format data sets with lof IC50 and AUC --------

  insert_msg('Wide-format IC50 and AUC data sets')

  for(i in c('log_ic50', 'auc')) {

    gdsc_drugs[[i]] <- gdsc_drugs$design %>%
      map(select, sample_id, variable, all_of(i)) %>%
      map(pivot_wider,
          id_cols = 'sample_id',
          names_from = 'variable',
          values_from = all_of(i))

  }

# Exploratory analysis: cancer entities and tissues -------

  insert_msg('Explorative analysis: tissues and cancer entities')

  ## counts of cancer entities according to classification of TCGA,
  ## tissue types and tissue subtypes

  gdsc_drugs$entity_counts <- gdsc_drugs$design %>%
    map(filter,
        !duplicated(sample_id)) %>%
    map(count, tcga_entity) %>%
    reduce(full_join, by = 'tcga_entity') %>%
    set_names(c('tcga_entity', 'n_gdsc1', 'n_gdsc2'))

  gdsc_drugs$tissue_counts <- gdsc_drugs$design %>%
    map(filter,
        !duplicated(sample_id)) %>%
    map(count, tissue) %>%
    reduce(full_join, by = 'tissue') %>%
    set_names(c('tissue', 'n_gdsc1', 'n_gdsc2'))

  gdsc_drugs$tissue_subtype_counts <- gdsc_drugs$design %>%
    map(filter,
        !duplicated(sample_id)) %>%
    map(count, tissue_subtype) %>%
    reduce(full_join, by = 'tissue_subtype') %>%
    set_names(c('tissue_subtype', 'n_gdsc1', 'n_gdsc2'))

# Exploratory analysis: common drugs and cell lines --------

  insert_msg('Common drugs and cell lines')

  ## drugs and cell lines shared by both GDSC experiments

  gdsc_drugs$common_drugs <- gdsc_drugs$log_ic50 %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  gdsc_drugs$common_samples <- gdsc_drugs$log_ic50 %>%
    map(~.x[['sample_id']]) %>%
    reduce(intersect)

# Exploratory analysis: systematic differences between experiments -------

  insert_msg('Systematic differences between experiments')

  ## systematic differences in log IC50 between the experiments
  ## investigated in common cell lines by paired T test

  ## data for testing

  gdsc_drugs$differences <- gdsc_drugs$log_ic50 %>%
    map(filter, sample_id %in% gdsc_drugs$common_samples) %>%
    map(select, sample_id, all_of(gdsc_drugs$common_drugs)) %>%
    compress(names_to = 'experiment') %>%
    mutate(experiment = factor(experiment, c('gdsc2', 'gdsc1'))) %>%
    arrange(sample_id, experiment)

  ## tests

  gdsc_drugs$differences <- gdsc_drugs$differences %>%
    test_two_groups(variables = gdsc_drugs$common_drugs,
                    split_fct = 'experiment',
                    type = 't',
                    adj_method = 'BH',
                    paired = TRUE) %>%
    mutate(regulation = ifelse(p_adjusted >= 0.05, 'ns',
                               ifelse(estimate > 0, 'upregulated',
                                      ifelse(estimate < 0,
                                             'downregulated', 'ns'))),
           regulation = factor(regulation,
                               c('upregulated', 'downregulated', 'ns')))

  ## volcano plot: log-difference and pFDR

  gdsc_drugs$difference_plot <-
    plot_volcano(gdsc_drugs$differences,
                 regulation_variable = 'estimate',
                 p_variable = 'p_adjusted',
                 signif_level = 0.05,
                 regulation_level = 0,
                 x_lab = '\u0394 logIC50, GDSC2 versus GDSC1',
                 y_lab = expression('-log'[10] * ' pFDR'),
                 plot_title = 'Systematic differences in drug sensitivity',
                 plot_subtitle = paste0('common compounds: n = ',
                                        nrow(gdsc_drugs$differences),
                                        '\nupregulated: n = ',
                                        table(gdsc_drugs$differences$regulation)[1],
                                        ', downregulated: n = ',
                                        table(gdsc_drugs$differences$regulation)[2]),
                 cust_theme = globals$common_theme +
                   theme(plot.tag = element_blank()),
                 label_variable = 'response',
                 top_regulated = 20,
                 label_type = 'text',
                 txt_size = 2.5)

# Caching the analysis results ------

  insert_msg('Caching the analysis results')

  save(gdsc_drugs, file = './data/GDSC/gdsc_drugs.RData')

# END ------

  rm(i)

  insert_tail()
