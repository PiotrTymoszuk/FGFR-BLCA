# Import of the cell line data: baseline expression of the FGFR-, FGF-, and
# FGFBP-coding genes, somatic mutations of FGFRs, and IC50 of FGFR inhibitors
# tested in the GDSC screening experiments.

  insert_head()

# tools --------

  library(figur)

# container -------

  depmap <- list()

# Cell line lexicon -------

  insert_msg('Cell line lexicon')

  depmap$cell_lines <- read_csv('./data/cell lines/Model.csv')

  ## clinical and demographic information

  depmap$cell_lines <- depmap$cell_lines %>%
    transmute(sample_id = ModelID,
              cosmic_id = COSMICID,
              cell_line_name = CellLineName,
              oncotree = DepmapModelType,
              lineage = OncotreeLineage,
              age = Age,
              sex = Sex,
              metastasis = ifelse(PrimaryOrMetastasis == 'Primary',
                                  'no', 'yes'),
              metastasis = factor(metastasis, c('no', 'yes')),
              patient_treatment = tolower(PatientTreatmentStatus),
              patient_treatment = factor(patient_treatment,
                                         c('pre-treatment', 'post-treatment')),
              systemic_chemotherapy = ifelse(stri_detect(PatientTreatmentType,
                                                         regex = 'chemo|Chemo'),
                                             'yes', 'no'),
              systemic_chemotherapy = factor(systemic_chemotherapy,
                                             c('no', 'yes')),
              pt_stage = paste0('T', Stage),
              pt_stage = factor(pt_stage, c('Ta', 'Tcis', paste0('T', 0:4))),
              pt_stage = droplevels(pt_stage))

  ## filtering for urothelial cancer cell lines

  depmap$cell_lines <- depmap$cell_lines %>%
    filter(lineage == 'Bladder/Urinary Tract')

# Clinical information lexicon for the cell lines --------

  insert_msg('Clinical information lexicon for the cell lines')

  depmap$clinic_lexicon <-
    tibble(variable = c('oncotree', 'age', 'sex',
                        'metastasis',
                        'patient_treatment',
                        'systemic_chemotherapy',
                        'pt_stage'),
           label = c('Oncotree', 'Age', 'Gender',
                     'Metastsis as a source',
                     'Patient treatment status',
                     'Systemic chemotherapy',
                     'pT stage')) %>%
    mutate(format = ifelse(variable %in% c('age'),
                           'numeric', 'factor'),
           unit = ifelse(variable == 'age', 'years', NA))

# Expression of FGFR-, FGF-, and FGFBP-coding genes --------

  insert_msg('Gene expression')

  ## provided as log2 TPM

  depmap$raw_expression <-
    read_csv('./data/cell lines/OmicsExpressionProteinCodingGenesTPMLogp1.csv')

  names(depmap$raw_expression)[1] <- 'sample_id'

  depmap$raw_expression <- depmap$raw_expression %>%
    filter(sample_id %in% depmap$cell_lines$sample_id)

  depmap$annotation <-
    tibble(variable = names(depmap$raw_expression)[-1]) %>%
    mutate(gene_symbol = stri_extract(variable, regex = '\\w+'),
           entrez_id = stri_extract(variable, regex = '\\(\\d+\\)'),
           entrez_id = stri_extract(entrez_id, regex = '\\d+')) %>%
    filter(gene_symbol %in% reduce(globals[c("receptors",
                                             "ligands",
                                             "binding_proteins")],
                                   union))

  depmap$expression <- depmap$raw_expression %>%
    select(sample_id, all_of(depmap$annotation$variable)) %>%
    set_names(c('sample_id', depmap$annotation$gene_symbol))

  names(depmap$raw_expression)[-1] <-
    stri_split_fixed(names(depmap$raw_expression)[-1],
                     pattern = ' ',
                     simplify = TRUE)[, 1]

# Somatic mutations of FGFRs --------

  insert_msg('Somatic mutations')

  depmap$mutation_detail <-
    read_csv('./data/cell lines/OmicsSomaticMutations.csv')

  depmap$mutation_detail <- depmap$mutation_detail %>%
    filter(!is.na(ProteinChange),
           ProteinChange != '',
           HugoSymbol %in% c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'))

  depmap$mutation <- extract_mutations(depmap$mutation_detail,
                                       sample_id = 'ModelID',
                                       gene_id = 'HugoSymbol')

  depmap$mutation <-
    full_join(depmap$expression[, "sample_id"],
              depmap$mutation,
              by = 'sample_id')

  depmap$mutation[c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4')] <-
    depmap$mutation[c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4')] %>%
    map_dfc(~ifelse(is.na(.x), 0, .x)) %>%
    map_dfc(~ifelse(.x == 0, 'no', 'yes')) %>%
    map_dfc(factor, c('no', 'yes'))

  depmap$mutation <- depmap$mutation %>%
    filter(sample_id %in% depmap$cell_lines$sample_id) %>%
    set_names(c('sample_id', paste0(names(depmap$mutation)[-1], '_mut')))

# CRISPR screening  ---------

  insert_msg('CRISPR screening')

  depmap$raw_crispr <- read_csv('./data/cell lines/CRISPRGeneEffect.csv')

  names(depmap$raw_crispr)[1] <- 'sample_id'

  names(depmap$raw_crispr) <-
    stri_split_fixed(names(depmap$raw_crispr), pattern = ' ', simplify = TRUE)[, 1]

  ## selecting the genes and cell lines of interest

  depmap$raw_crispr <- depmap$raw_crispr %>%
    filter(sample_id %in% depmap$cell_lines$sample_id)

  depmap$crispr <- depmap$raw_crispr %>%
    select(sample_id,
           any_of(reduce(globals[c("receptors", "ligands", "binding_proteins")],
                         union)))

  names(depmap$crispr)[-1] <-
    paste0(names(depmap$crispr)[-1], '_crispr')

# RNAi screening ---------

  insert_msg('RNAi gene effects')

  depmap$raw_rnai <-
    read_csv('./data/cell lines/RNAi_(Achilles+DRIVE+Marcotte,_DEMETER2)_subsetted.csv')

  names(depmap$raw_rnai)[1] <- 'sample_id'

  ## selecting cell lines and genes of interest

  depmap$raw_rnai <- depmap$raw_rnai %>%
    filter(sample_id %in% depmap$cell_lines$sample_id)

  depmap$rnai <- depmap$raw_rnai %>%
    select(sample_id,
           any_of(reduce(globals[c("receptors", "ligands", "binding_proteins")],
                         union)))

# Backgrounds of gene effects: gene effects for non-expressed genes --------

  insert_msg('Backgrounds of gene effects')

  ## non-expressed genes: 95th percentile of log2 gene counts is equal to 0

  depmap$non_expressed <- depmap$raw_expression[-1] %>%
    colQuantiles(0.95)

  depmap$non_expressed <-
    depmap$non_expressed[depmap$non_expressed == 0, , drop = FALSE]

  depmap$non_expressed <- rownames(depmap$non_expressed)

  ## the corresponding gene effects

  depmap[c("crispr_bcg", "rnai_bcg")] <-
    depmap[c("raw_crispr", "raw_rnai")] %>%
    map(select, sample_id, any_of(depmap$non_expressed))

# Drug resistance, GDSC -------

  insert_msg('Drug resistance')

  ## drugs of interest (source: https://www.cancerrxgene.org/compounds)
  ## doxorubicin serves as a positive control

  depmap$gdsc_interest <-
    tibble(drug_name = c('AZD4547',
                         'FGFR_0939',
                         'FGFR_3831',
                         'PD173074',
                         'Doxorubicin'),
           specificity = c('FGRF1, FGFR2, FGFR',
                           'FGFR4',
                           'FGFR1, FGFR2, FGFR3, FGFR4',
                           'FGFR1, FGFR2, FGFR3',
                           NA),
           max_conc = c(10, 2.5, 10, 10, 5))

  depmap$gdsc_lexicon <-
    paste0('./data/cell lines/', c('GDSC1.csv', 'GDSC2.csv')) %>%
    map_dfr(read_csv)

  depmap$gdsc_lexicon <- depmap$gdsc_lexicon %>%
    transmute(cosmic_id = `Cosmic ID`,
              drug_id = `Drug ID`,
              drug_name = `Drug Name`,
              data_set = `Dataset Version`,
              variable = paste(data_set, drug_id, drug_name, sep = '_'),
              log2_ic50 = log2(exp(IC50))) %>%
    filter(drug_name %in% depmap$gdsc_interest$drug_name,
           cosmic_id %in% depmap$cell_lines$cosmic_id) %>%
    left_join(depmap$cell_lines[, c("sample_id", "cell_line_name", "cosmic_id")],
              by = 'cosmic_id')

  ## drug log2 IC50 in a wide format

  depmap$gdsc <- depmap$gdsc_lexicon %>%
    select(sample_id, variable, log2_ic50) %>%
    pivot_wider(id_cols = 'sample_id',
                names_from = 'variable',
                values_from = 'log2_ic50')

  depmap$gdsc_lexicon <- depmap$gdsc_lexicon %>%
    select(-cosmic_id, -cell_line_name, -sample_id, -log2_ic50) %>%
    filter(!duplicated(variable)) %>%
    left_join(depmap$gdsc_interest, by = 'drug_name')

# Drug resistance, PRISM ---------

  insert_msg('Drug resistance, PRISM')

  depmap$prism_lexicon <-
    read_csv('./data/cell lines/secondary-screen-dose-response-curve-parameters.csv')

  depmap$prism_lexicon <- depmap$prism_lexicon %>%
    transmute(sample_id = depmap_id,
              data_set = 'PRISM',
              moa = moa,
              drug_id = broad_id,
              drug_name = name,
              variable = paste(data_set, drug_id, drug_name, sep = '_'),
              log2_ic50 = log2(ic50),
              auc = auc) %>%
    filter(sample_id %in% depmap$cell_lines$sample_id) %>%
    filter((moa %in% c('FGFR inhibitor', 'FGFR antagonist')) |
             drug_name == 'doxorubicin')

  ## drug AUC in a wide format

  depmap$prism <- depmap$prism_lexicon %>%
    select(sample_id, variable, auc) %>%
    pivot_wider(id_cols = 'sample_id',
                names_from = 'variable',
                values_from = 'auc')

  depmap$prism_lexicon <- depmap$prism_lexicon %>%
    filter(!duplicated(variable)) %>%
    select(-sample_id, -log2_ic50, -auc) %>%
    mutate(specificity = NA,
           max_conc = NA)

# Common variable lexicon ----------

  insert_msg('Common variable lexicon')

  depmap$lexicon[c('gdsc', 'prism')] <-
    depmap[c('gdsc_lexicon', 'prism_lexicon')] %>%
    list(x = .,
         y = c('GDSC', 'PRISM'),
         z = c('log2 IC50, µM', 'AUC')) %>%
    pmap(function(x, y, z) x %>%
           transmute(x,
                     variable = variable,
                     label = paste(drug_name, drug_id, sep = ', ID: '),
                     label_html = label,
                     specificity = specificity,
                     max_conc = max_conc,
                     data_set = y,
                     unit = z,
                     type = 'resistance',
                     modeling = 'response'))

  depmap$lexicon$crispr <-
    tibble(variable = names(depmap$crispr)[-1],
           label = names(depmap$crispr)[-1] %>%
             stri_replace(regex = '_crispr$',
                          replacement = '') %>%
             paste('gene effect'),
           label_html = names(depmap$crispr)[-1] %>%
             stri_replace(regex = '_crispr$',
                          replacement = '') %>%
             html_italic %>%
             paste('gene effect'),
           unit = NA,
           type = 'crispr',
           modeling = 'response',
           data_set = 'CHRONOS')

  depmap$lexicon$expression <- depmap$annotation %>%
    transmute(variable = gene_symbol,
              label = paste(gene_symbol, 'expression'),
              label_html = paste(html_italic(gene_symbol), 'expression'),
              specificity = NA,
              unit = 'log2 TPM',
              type = 'expression',
              modeling = 'explanatory',
              data_set = NA)

  depmap$lexicon$mutation <-
    tibble(variable = names(depmap$mutation)[-1],
           label = names(depmap$mutation)[-1] %>%
             stri_extract(regex = '^FGFR\\d{1}') %>%
             paste('mutation'),
           label_html = names(depmap$mutation)[-1] %>%
             stri_extract(regex = '^FGFR\\d{1}') %>%
             html_italic %>%
             paste('mutation'),
           specificity = NA,
           unit = NA,
           type = 'mutation',
           modeling = 'explanatory',
           data_set = NA)

  depmap$lexicon <- depmap$lexicon %>%
    reduce(full_rbind) %>%
    mutate(type = factor(type,
                         c('resistance', 'crispr', 'expression', 'mutation')),
           modeling = factor(modeling,
                             c('response', 'explanatory')))

  ## axis and table labels

  depmap$lexicon <- depmap$lexicon %>%
    mutate(axis_lab = ifelse(is.na(unit),
                             label,
                             paste(label, unit, sep = ', ')),
           table_lab = ifelse(is.na(data_set),
                              axis_lab,
                              paste(data_set, axis_lab, sep = ', ')))

# caching -------

  insert_msg('Caching')

  depmap <-
    depmap[c("cell_lines", "clinic_lexicon",
             "raw_expression", "expression", "annotation",
             "mutation",
             "crispr", "crispr_bcg",
             "rnai", "rnai_bcg",
             "gdsc_lexicon", "gdsc",
             "prism_lexicon", "prism",
             "lexicon")]

  save(depmap, file = './data/depmap.RData')

# END -------

  insert_tail()
