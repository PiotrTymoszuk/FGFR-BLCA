# Import of the GENIE BLA and urothelial cohort data.
#
# In case of multiple samples per patient, the first one is included in the
# analysis.
# Approximate survival is calculated based on age at sequencing and time to the
# last contact or death

  insert_head()

# container ------

  genie <- list()

# raw data fetched from API --------

  insert_msg('Raw data fetched from API')

  load('./data/GENIE/genie_raw.RData')

# Clinical information for the patients --------

  insert_msg('Clinical information: patient')

  genie$clinic$patient <- genie_raw$patient %>%
    select(patientId, clinicalAttributeId, value) %>%
    pivot_wider(id_cols = 'patientId',
                names_from = 'clinicalAttributeId',
                values_from = 'value')

  genie$clinic$patient <- genie$clinic$patient %>%
    transmute(patient_id = patientId,
              center = factor(CENTER),
              death = as.numeric(as.logical(toupper(DEAD))),
              sex = factor(tolower(SEX), c('female', 'male')),
              sample_count = as.numeric(SAMPLE_COUNT),
              ## age at the last contact and death in days
              last_contact_days = as.numeric(INT_CONTACT),
              death_days = as.numeric(INT_DOD))

# Clinical information for the samples -------

  insert_msg('Clinical information for the samples')

  ## only primary tumors!

  genie$clinic$sample <- genie_raw$sample %>%
    select(sampleId, patientId, clinicalAttributeId, value) %>%
    pivot_wider(id_cols = c('sampleId', 'patientId'),
                names_from = 'clinicalAttributeId',
                values_from = 'value') %>%
    filter(SAMPLE_TYPE_DETAILED == 'Primary tumor',
           ONCOTREE_CODE %in% c('BLAD',
                                'BLCA',
                                'BLSC',
                                'UA', 'UAD', 'UCA', 'UCU',
                                'UPA', 'URCA', 'USCC',
                                'UTUC'))

  ## TMB: in this case this is the fraction genome altered

  genie$clinic$sample <- genie$clinic$sample %>%
    transmute(sample_id = sampleId,
              patient_id = patientId,
              age = stri_extract(AGE_AT_SEQ_REPORT, regex = '\\d+'),
              age = as.numeric(age),
              age_days = age * 365.25,
              oncotree = ONCOTREE_CODE,
              tissue = ifelse(stri_detect(ONCOTREE_CODE, regex = '^BL'),
                              'bladder', 'non-bladder'),
              tissue = factor(tissue, c('bladder', 'non-bladder')),
              mutation_count = as.numeric(MUTATION_COUNT),
              tmb_per_mb = as.numeric(FRACTION_GENOME_ALTERED))

# A common clinical information frame, multiplicates and survival ------

  insert_msg('A common clinical information frame, multiplicates, survival')

  ## multiplicates: the first sample per patient is included in  teh analysis

  genie$clinic <- genie$clinic %>%
    reduce(inner_join, by = 'patient_id') %>%
    arrange(patient_id, sample_id) %>%
    blast(patient_id) %>%
    map_dfr(~.x[1, ])

  genie$clinic <- genie$clinic %>%
    mutate(os_days = ifelse(!is.na(death_days),
                            death_days - age_days,
                            last_contact_days - age_days),
           os_days = ifelse(os_days < 0, NA, os_days),
           os_months = os_days/30.437)

# Mutation data -------

  insert_msg('Mutation data')

  genie$mutation_detail <- genie_raw$mutations

  ## a simplified data frame with 0/1 coded non-silent mutations

  genie$mutation <- genie$mutation_detail %>%
    extract_mutations(id_variable = 'sampleId',
                      gene_id = 'hugoGeneSymbol')

  ## detailed mutations for FGFR1 - 4

  genie$mutation_detail <- genie$mutation_detail %>%
    filter(hugoGeneSymbol %in% paste0('FGFR', 1:4)) %>%
    mutate(sample_id = sampleId)

# Copy number alterations -------

  insert_msg('Copy number alterations')

  genie[c('deletion', 'amplification')] <-
    list(filter(genie_raw$cna, alteration < 0),
         filter(genie_raw$cna, alteration > 0)) %>%
    map(transmute,
        ID = sampleId,
        gene_symbol = hugoGeneSymbol,
        alteration = alteration) %>%
    map(pivot_wider,
        id_cols = 'ID',
        names_from = 'gene_symbol',
        values_from = 'alteration') %>%
    map(column_to_rownames, 'ID') %>%
    map(as.matrix)

  ## padding with NA for samples absent from the CNA frames
  ## it's assumed that those specimens do not have gene amplifications
  ## or deletions

  genie[c('pad_losses', 'pad_gains')] <-
    genie[c('deletion', 'amplification')] %>%
    map(rownames) %>%
    map(~filter(genie$mutation, !sample_id %in% .x)) %>%
    map(~.x$sample_id) %>%
    map2(., genie[c('deletion', 'amplification')],
         ~matrix(NA,
                 nrow = length(.x),
                 ncol = ncol(.y),
                 dimnames = list(.x, colnames(.y))))

  genie[c("deletion", "amplification")] <-
    map2(genie[c("deletion", "amplification")],
         genie[c("pad_losses", "pad_gains")],
         rbind) %>%
    map(~ifelse(is.na(.x), 0, 1))

  ## data frame conversion, gene variables as factors

  genie[c("deletion", "amplification")] <-
    genie[c("deletion", "amplification")] %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Complete cases -------

  insert_msg('Complete cases')

  genie$complete_ids <-
    genie[c("clinic", "mutation", "deletion", "amplification")] %>%
    map(~.x$sample_id) %>%
    reduce(intersect)

  genie[c("clinic",
          "mutation", "mutation_detail",
          "deletion", "amplification")] <-
    genie[c("clinic",
            "mutation", "mutation_detail",
            "deletion", "amplification")] %>%
    map(filter, sample_id %in% genie$complete_ids)

# Caching the results ------

  insert_msg('Caching the results')

  genie <-
    genie[c("clinic",
            "mutation", "mutation_detail",
            "deletion", "amplification")]

  save(genie, file = './data/genie.RData')

# END -----

  rm(genie_raw)

  insert_tail()
