# Import of clinical and genetic data for the MSK cohort

  insert_head()

# container ------

  msk <- list()

# clinical data -------

  insert_msg('Clinical data')

  ## patient

  msk$clinic$patient <-
    read_tsv('./data/MSK 2022/data_clinical_patient.txt',
             skip = 4)

  msk$clinic$patient <- msk$clinic$patient %>%
    transmute(patient_id = PATIENT_ID,
              race = factor(RACE),
              sex = factor(tolower(SEX), c('female', 'male')),
              death = ifelse(OS_STATUS == 'DECEASED', 1, 0),
              os_months = as.numeric(OS_MONTHS),
              os_days = os_months * 30.437,
              metastasis_location = car::recode(MET_LOCATION,
                                                "'Distant' = 'distant';
                                                'Ln' = 'LN'"),
              metastasis_location = factor(metastasis_location,
                                           c('LN', 'distant')),
              smoking_history = car::recode(SMOKER,
                                            "'Active' = 'yes';
                                            'Former' = 'yes';
                                            'Never' = 'no'"),
              smoking_history = factor(smoking_history, c('no', 'yes')))

  ## sample: restricting to primary cancers, in case of duplicated samples
  ## the first is taken for the analysis

  msk$clinic$sample <-
    read_tsv('./data/MSK 2022/data_clinical_sample.txt',
             skip = 4) %>%
    filter(SAMPLE_TYPE == 'Primary')

  msk$clinic$sample <- msk$clinic$sample %>%
    transmute(sample_id = SAMPLE_ID,
              patient_id = PATIENT_ID,
              age = as.numeric(AGE_AT_DX),
              tissue_type = factor('tumor', c('normal', 'tumor')),
              pt_stage = paste0('T', SPECIMEN_STAGE),
              pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')),
              invasiveness = ifelse(pt_stage == 'T1',
                                    'non-muscle invasive', 'muscle invasive'),
              invasiveness = factor(invasiveness,
                                    c('non-muscle invasive', 'muscle invasive')),
              tissue = ifelse(PRIMARY_SITE %in% c('bladder', 'Bladder'),
                              'bladder', 'non-bladder'),
              tissue = factor(tissue, c('bladder', 'non-bladder')),
              tm_purity = as.numeric(TUMOR_PURITY),
              msi_score = as.numeric(MSI_SCORE),
              intravesical_treatment = ifelse(INTRAVESICAL_TREATMENT == 0,
                                              'no', 'yes'),
              intravesical_treatment = factor(intravesical_treatment,
                                              c('no', 'yes')),
              systemic_treatment = ifelse(SYSTEMIC_TREATMENT == 0,
                                          'no', 'yes'),
              systemic_treatment = factor(systemic_treatment,
                                          c('no', 'yes')),
              tmb_per_mb = as.numeric(TMB_NONSYNONYMOUS))

  msk$clinic$sample <- msk$clinic$sample %>%
    arrange(patient_id, sample_id) %>%
    blast(patient_id) %>%
    map_dfr(~.x[1, ])

  ## merging the samples and patients together,
  ## removal of pediatric cases/nonsense age entries

  msk$clinic <- msk$clinic %>%
    reduce(inner_join, by = 'patient_id') %>%
    filter(age >= 18)

# Mutation data -------

  insert_msg('Mutation data')

  ## a detailed data frame with mutations to be used later for
  ## an analysis of FGFR1 - 4 mutations

  msk$mutation_detail <-
    read_tsv('./data/MSK 2022/data_mutations.txt') %>%
    filter(Variant_Classification != 'Silent')

  ## simplified mutation data: 0 codes for WT and 1 for
  ## any protein-affecting mutation

  msk$mutation <- msk$mutation_detail %>%
    extract_mutations

  ## restricting of the detailed mutation frame with FGFR1 - 4

  msk$mutation_detail <- msk$mutation_detail %>%
    filter(Hugo_Symbol %in% paste0('FGFR', 1:4)) %>%
    mutate(sample_id = Tumor_Sample_Barcode)

# Copy number alterations --------

  insert_msg('CNA')

  msk$cna <- read_tsv('./data/MSK 2022/data_cna.txt')

  msk$cna <- msk$cna %>%
    filter(!duplicated(Hugo_Symbol)) %>%
    column_to_rownames('Hugo_Symbol') %>%
    select(any_of(msk$clinic$sample_id)) %>%
    t

  ## separate data frames for homozyguous deletions and amplifications
  ## 0/1 coding. NA: assumed no copy number alteration

  msk$deletion <- ifelse(msk$cna < -1, 1, 0)
  msk$amplification <- ifelse(msk$cna > 1, 1, 0)

  msk[c("deletion", "amplification")] <-
    msk[c("deletion", "amplification")] %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

  for(i in c('deletion', 'amplification')) {

    msk[[i]][, -1] <- msk[[i]][, -1] %>%
      map_dfc(~ifelse(is.na(.x), 0, .x))

  }

# Restricting the data sets to the patients with complete information ------

  insert_msg('Complete cases')

  msk$complete_ids <-
    msk[c("clinic", "mutation", "deletion", "amplification")] %>%
    map(~.x$sample_id) %>%
    reduce(intersect)

  msk[c("clinic",
         "mutation", "mutation_detail",
         "deletion", "amplification")] <-
    msk[c("clinic",
           "mutation", "mutation_detail",
           "deletion", "amplification")] %>%
    map(filter, sample_id %in% msk$complete_ids)

# Caching the results -------

  insert_msg('Caching the results')

  msk <- msk[c("clinic",
               "mutation_detail", "mutation",
               "deletion", "amplification")]

  save(msk, file = './data/msk.RData')

# END ------

  insert_tail()

