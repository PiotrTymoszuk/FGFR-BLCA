# Import of the BCAN study data: clinical information, 0/1-coded gene-level
# mutation information, and Z-scores of log2 expression counts.

  insert_head()

# container -------

  bcan <- list()

# clinical information -------

  insert_msg('Clinical information')

  bcan$clinic$patient <-
    read_tsv('./data/BCAN/data_clinical_patient.txt', skip = 4)

  ## patient's data

  bcan$clinic$patient <- bcan$clinic$patient %>%
    transmute(
              ## demography

              patient_id = PATIENT_ID,
              sex = factor(tolower(SEX), c('female', 'male')),
              race = factor(RACE),
              race = fct_recode(race, Other = 'Unknown'),
              age = as.numeric(AGE_AT_DIAGNOSIS),
              smoking = factor(SMOKING_STATUS),
              smoking = fct_recode(smoking,
                                   `current or previous` = 'Current Smoker',
                                   `current or previous` = 'Former Smoker',
                                   never = 'Never Smoker'),
              smoking = factor(smoking, c('never', 'current or previous')),
              baseline_ecog = factor(BASELINE_ECOG),

              ## anatomy

              primary_location = stri_extract_first(PRIMARY_TUMOR_LOCATION,
                                                    regex = 'Bladder|Ureter|Urethra|(Renal\\s{1}pelvis)'),
              primary_location = factor(tolower(primary_location)),
              primary_location = fct_relevel(primary_location, 'bladder'),
              tissue = factor(tolower(SAMPLE_SITE)),
              tissue = fct_relevel(tissue, 'bladder'),
              tissue_type = factor(tolower(PRIMARY_MET)),
              tissue_type = fct_recode(tissue_type,
                                       tumor = 'primary',
                                       metastasis = 'met'),
              tissue_type = fct_relevel(tissue_type, 'tumor', 'metastasis'),

              ## surgery

              cystectomy = ifelse(stri_detect(PRIMARY_SURGERY,
                                              regex = '(c|C)yst(o|e)'),
                                  'yes', 'no'),
              cystectomy = factor(cystectomy, c('no', 'yes')),
              turbt = ifelse(stri_detect(PRIMARY_SURGERY, regex = 'TURBT'),
                             'yes', 'no'),
              turbt = factor(turbt, c('no', 'yes')),

              ## stages are restricted to the main ones

              pt_stage = stri_extract(PRIMARY_SURGERY_T_STAGE, regex = '(T\\d{1})|Ta|Tis'),
              pt_stage = factor(pt_stage,
                                c('T0', 'Ta', 'Tis', 'T1', 'T2', 'T3', 'T4')),
              invasiveness = ifelse(pt_stage %in% c('T2', 'T3', 'T4'),
                                    'muscle invasive', 'non-muscle invasive'),
              invasiveness = factor(invasiveness,
                                    c('muscle invasive', 'non-muscle invasive')),
              pn_stage = stri_extract(PRIMARY_SURGERY_N_STAGE, regex = 'N\\d{1}'),
              pn_stage = factor(pn_stage, paste0('N', 0:3)),
              pm_stage = stri_extract(PRIMARY_SURGERY_M_STAGE, regex = 'M\\d{1}'),
              pm_stage = factor(pm_stage, c('M0', 'M1')),

              ## therapy

              neoadjuvant = factor(tolower(NEOADJUVANT_THERAPY),
                                   c('no', 'yes')),
              radiation = factor(tolower(PRIMARY_RADIATION_THERAPY),
                                 c('no', 'yes')),
              adjuvant_therapy = factor(tolower(ADJUVANT_THERAPY),
                                        c('no', 'yes')),
              meta_radiation = factor(tolower(METASTATIC_RADIATION_THERAPY),
                                      c('no', 'yes')),
              systemic_chemotherapy = ifelse(is.na(CHEMOTHERAPY),
                                             'no', tolower(CHEMOTHERAPY)),
              systemic_chemotherapy = factor(systemic_chemotherapy,
                                             c('no', 'yes')),
              chemotherapy_bor = factor(BEST_RESPONSE_CHEMOTHERAPY),
              chemotherapy_bor = fct_recode(chemotherapy_bor,
                                            CR = 'Complete Response',
                                            PR = 'Partial Response',
                                            PD = 'Progressive Disease',
                                            SD = 'Stable Disease'),
              chemotherapy_bor = factor(chemotherapy_bor,
                                        c('CR', 'PR', 'PD', 'SD')),
              chemotherapy_biresponse = fct_recode(chemotherapy_bor,
                                                   `CR/PR` = 'CR',
                                                   `CR/PR` = 'PR',
                                                   `SD/PD` = 'PD',
                                                   `SD/PD` = 'SD'),
              immunotherapy = ifelse(`ADC+IMMUNOTHERAPY` == 'Yes' |
                                       `CHEMOTHERAPY+IMMUNOTHERAPY` == 'Yes',
                                     'yes', 'no'),
              immunotherapy = ifelse(is.na(immunotherapy), 'no', immunotherapy),

              ## survival: chemo and immuno

              death = ifelse(stri_detect(SURVIVAL_STATUS,
                                         regex = '^Dea'),
                             1, 0),
              tumor_death = ifelse(stri_detect(DEATH_REASON,
                                               regex = 'disease'),
                                   1, 0),
              os_days = 30.437 * as.numeric(SURVIVAL_TIME),
              tss_days = os_days,
              chemo_os_days = 30.437 * as.numeric(CHEMO_SURVIVAL),
              immuno_os_days = 30.437 * as.numeric(IMMUNO_SURVIVAL))

  ## sample data: essentially only sample ID

  bcan$clinic$sample <-
    read_tsv('./data/BCAN/data_clinical_sample.txt', skip = 4) %>%
    set_names(c('sample_id', 'patient_id', 'oncotree'))

  bcan$clinic <- bcan$clinic %>%
    reduce(inner_join, by = 'patient_id') %>%
    relocate(sample_id, patient_id)

# Expression Z-scores -------

  insert_msg('Expression Z-scores')

  bcan$raw_expression <-
    read_tsv('./data/BCAN/data_RNA_Seq_v2_mRNA_median_Zscores.txt') %>%
    mutate(probe_id = paste0('probe_', 1:nrow(.))) %>%
    column_to_rownames('probe_id')

  ## annotation and missingness of the gene expression information

  bcan$annotation <- bcan$raw_expression %>%
    transmute(entrez_id = as.character(Entrez_Gene_Id)) %>%
    mutate(probe_id = rownames(bcan$raw_expression),
           n_missing = unname(rowMissing(bcan$raw_expression[-1])),
           percent_missing = n_missing/ncol(bcan$raw_expression) * 100,
           gene_symbol = mapIds(org.Hs.eg.db,
                                keys = entrez_id,
                                keytype = 'ENTREZID',
                                column = 'SYMBOL')) %>%
    select(probe_id, gene_symbol, entrez_id, n_missing, percent_missing) %>%
    filter(complete.cases(.)) %>%
    as_tibble

  ## formatting the expression data frame

  bcan$raw_expression <- integrate_expression(bcan$raw_expression[-1],
                                              bcan$annotation)

  ## removal of samples with almost all missing samples

  bcan$sample_missing <- bcan$raw_expression %>%
    column_to_rownames('sample_id') %>%
    rowMissing %>%
    compress(names_to = 'sample_id',
             values_to = 'n_missing') %>%
    mutate(percent_missing = n_missing/ncol(bcan$raw_expression) * 100)

  bcan$messy_sample_id <- bcan$sample_missing %>%
    filter(percent_missing >= 80) %>%
    .$sample_id

  bcan$raw_expression <- bcan$raw_expression %>%
    filter(!sample_id %in% bcan$messy_sample_id)

# Imputation of the expression data set: variables with < 50% missingness ---------

  insert_msg('Imputation of the expression data set')

  ## the imputed data set will be used for
  ## predictions of the cluster assignment by an Elastic Net classifier
  ##
  ## done for proteins with less than 50% missing values

  bcan$impute_vars <- bcan$annotation %>%
    filter(percent_missing < 50) %>%
    .$gene_symbol %>%
    unname

  bcan$expression <- bcan$raw_expression %>%
    column_to_rownames('sample_id') %>%
    select(all_of(bcan$impute_vars)) %>%
    t %>%
    impute.knn(k = 9)

  bcan$expression <- bcan$expression$data %>%
    t %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

# Mutation data -------

  insert_msg('Mutation data')

  bcan$mutation_detail <- read_tsv('./data/BCAN/data_mutations.txt')

  bcan$mutation <- bcan$mutation %>%
    filter(!is.na(Protein_position),
           !stri_detect(Hugo_Symbol, regex = '^\\d{1}-Sep')) %>%
    extract_mutations(sample_id = 'Tumor_Sample_Barcode',
                      gene_id = 'Hugo_Symbol')

  bcan$mutation_detail <- bcan$mutation_detail %>%
    filter(Hugo_Symbol %in% paste0('FGFR', 1:4)) %>%
    mutate(sample_id = Tumor_Sample_Barcode,
           Protein_position = stri_split_fixed(Protein_position,
                                               pattern = '/',
                                               simplify = TRUE)[, 1],
           Protein_position = as.numeric(Protein_position))

# Restricting of to samples with the complete expression information -------

  insert_msg('Samples with complete expression information')

  bcan[c("raw_expression",
         "expression")] <-
    bcan[c("raw_expression",
           "expression")] %>%
    map(~inner_join(bcan$clinic[, c("sample_id", "patient_id", "tissue_type")],
                    .x,
                    by = 'sample_id'))

  bcan[c("clinic",
         "raw_expression",
         "expression",
         "mutation")] <-
    bcan[c("clinic",
           "raw_expression",
           "expression",
           "mutation")] %>%
    map(filter, sample_id %in% bcan$expression$sample_id)

# Caching --------

  insert_msg('Caching')

  bcan <- bcan[c("clinic",
                 "raw_expression",
                 "annotation",
                 "expression",
                 "mutation",
                 "mutation_detail")]

  save(bcan, file = './data/bcan.RData')

# END -----

  insert_tail()
