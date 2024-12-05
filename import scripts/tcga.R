# Import of the TCGA BLCA data set

  insert_head()

# container ------

  tcga <- list()

# Expression data and annotation -----

  insert_msg('Expression and annotation')

  tcga$expression <- read_tsv('./data/TCGA PanCancer/data_mrna_seq_v2_rsem.txt')

  ## annotation, manual correction for the 3730 gene (Entrez ID)

  tcga$annotation <- tcga$expression %>%
    transmute(probe_id = paste0('probe_', 1:nrow(.)),
              gene_symbol = Hugo_Symbol,
              entrez_id = Entrez_Gene_Id) %>%
    filter(complete.cases(.))

  tcga$annotation <- tcga$annotation %>%
    mutate(gene_symbol = ifelse(entrez_id == '3730',
                                'ANOS1', gene_symbol))

  ## integration of expression by arithmetic mean

  tcga$expression <- tcga$expression %>%
    select(-Hugo_Symbol, -Entrez_Gene_Id) %>%
    mutate(probe_id = paste0('probe_', 1:nrow(.))) %>%
    column_to_rownames('probe_id') %>%
    integrate_expression(annotation = tcga$annotation)

  ## patient's ID and tissue type: I'm working anyway only
  ## with the tumor samples of the bladder

  tcga$expression <- tcga$expression %>%
    column_to_rownames('sample_id') %>%
    as.matrix

  tcga$expression <- log2(tcga$expression + 1)

  tcga$expression <- tcga$expression %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

  tcga$expression <- tcga$expression %>%
    mutate(tissue_type = factor('tumor', c('tumor', 'normal')),
           patient_id = stri_extract(sample_id, regex = '^TCGA-\\w{2}-\\w{4}')) %>%
    relocate(tissue_type) %>%
    relocate(patient_id) %>%
    relocate(sample_id)

# Clinical and survival data -------

  insert_msg('Clinical and survival data')

  ## patient's data

  tcga$clinic$patient <-
    read_tsv('./data/TCGA PanCancer/data_clinical_patient.txt', skip = 4)

  tcga$clinic$patient <- tcga$clinic$patient %>%
    transmute(patient_id = PATIENT_ID,
              age = as.numeric(AGE),
              sex = factor(tolower(SEX),
                           c('female', 'male')),
              race = factor(RACE),
              p_stage = stri_replace(AJCC_PATHOLOGIC_TUMOR_STAGE,
                                     regex = '^STAGE\\s{1}',
                                     replacement = ''),
              p_stage = factor(p_stage, c('I', 'II', 'III', 'IV')),
              neoadjuvant = factor(tolower(HISTORY_NEOADJUVANT_TRTYN),
                                   c('no', 'yes')),
              location = car::recode(ICD_10,
                                     "'C67.0' = 'trigone/bladder';
                                     'C67.1' = 'dome/bladder';
                                     'C67.2' = 'lateral wall/bladder';
                                     'C67.3' = 'anterior wall/bladder';
                                     'C67.4' = 'posterior wall/bladder';
                                     'C67.5' = 'neck/bladder';
                                     'C67.6' = 'ureteric orifice/bladder';
                                     'C67.9' = 'unspecified/bladder'"),
              tissue = factor('bladder', c('non-bladder', 'bladder')),
              pm_stage = stri_extract(PATH_M_STAGE, regex = 'M\\d{1}'),
              pm_stage = factor(pm_stage),
              pn_stage = stri_extract(PATH_N_STAGE, regex = 'N\\d{1}'),
              pn_stage = factor(pn_stage),
              pt_stage = stri_extract(PATH_T_STAGE, regex = 'T\\d{1}'),
              pt_stage = factor(pt_stage),
              invasiveness = ifelse(pt_stage %in% c('T0', 'T1'),
                                    'non-muscle invasive',
                                    ifelse(is.na(pt_stage),
                                           NA, 'muscle invasive')),
              invasiveness = factor(invasiveness,
                                    c('non-muscle invasive', 'muscle invasive')),
              radiation = factor(tolower(RADIATION_THERAPY),
                                 c('no', 'yes')),
              death = stri_extract(OS_STATUS, regex = '^\\d{1}'),
              death = as.numeric(death),
              os_days = as.numeric(OS_MONTHS) * 30.436875,
              tumor_death = stri_extract(DSS_STATUS, regex = '^\\d{1}'),
              tumor_death = as.numeric(tumor_death),
              tss_days = as.numeric(DSS_MONTHS) * 30.436875,
              relapse = stri_extract(PFS_STATUS, regex = '^\\d{1}'),
              relapse = as.numeric(relapse),
              rfs_days = as.numeric(PFS_MONTHS) * 30.436875)

  ## sample data

  tcga$clinic$sample <-
    read_tsv('./data/TCGA PanCancer/data_clinical_sample.txt', skip = 4)

  tcga$clinic$sample <- tcga$clinic$sample %>%
    transmute(sample_id = SAMPLE_ID,
              patient_id = PATIENT_ID,
              hist_grade = stri_extract(GRADE, regex = 'Low|High'),
              hist_grade = factor(tolower(hist_grade), c('low', 'high')),
              aneuploidy_score = as.numeric(ANEUPLOIDY_SCORE),
              msi_score = as.numeric(MSI_SCORE_MANTIS),
              tmb_per_mb = as.numeric(TMB_NONSYNONYMOUS))

  ## a common frame with the clinical information

  tcga$clinic <- tcga$clinic %>%
    reduce(inner_join, by = 'patient_id')

# Mutation data ------

  insert_msg('Mutation data')

  ## a detailed data frame with mutations to be used later for
  ## an analysis of FGFR1 - 4 mutations

  tcga$mutation_detail <-
    read_tsv('./data/TCGA PanCancer/data_mutations.txt') %>%
    filter(Variant_Classification != 'Silent')

  ## simplified mutation data: 0 codes for WT and 1 for
  ## any protein-affecting mutation

  tcga$mutation <- tcga$mutation_detail %>%
    extract_mutations

  ## restricting of the detailed mutation frame with FGFR1 - 4

  tcga$mutation_detail <- tcga$mutation_detail %>%
    filter(Hugo_Symbol %in% paste0('FGFR', 1:4)) %>%
    mutate(sample_id = Tumor_Sample_Barcode)

# Copy number alterations --------

  insert_msg('CNA')

  tcga$cna <- read_tsv('./data/TCGA PanCancer/data_cna.txt')

  tcga$cna <- tcga$cna %>%
    filter(!duplicated(Hugo_Symbol)) %>%
    select(-Entrez_Gene_Id) %>%
    column_to_rownames('Hugo_Symbol') %>%
    t

  ## separate data frames for homozyguous deletions and amplifications
  ## 0/1 coding

  tcga$deletion <- ifelse(tcga$cna < -1, 1, 0)
  tcga$amplification <- ifelse(tcga$cna > 1, 1, 0)

  tcga[c("deletion", "amplification")] <-
    tcga[c("deletion", "amplification")] %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Restricting the data sets to the patients with complete information ------

  insert_msg('Complete cases')

  tcga$complete_ids <-
    tcga[c("clinic", "expression", "mutation", "deletion", "amplification")] %>%
    map(~.x$sample_id) %>%
    reduce(intersect)

  tcga[c("clinic", "expression",
         "mutation", "mutation_detail",
         "deletion", "amplification")] <-
    tcga[c("clinic", "expression",
           "mutation", "mutation_detail",
           "deletion", "amplification")] %>%
    map(filter, sample_id %in% tcga$complete_ids)

# Caching the results -------

  insert_msg('Caching the results')

  tcga <- tcga[c("clinic",
                 "expression", "annotation",
                 "mutation_detail", "mutation",
                 "deletion", "amplification")]

  save(tcga, file = './data/tcga.RData')

# END ------

  insert_tail()
