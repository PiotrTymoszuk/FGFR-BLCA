# Container -------

  imvigor <- list()

# Reading the data sets -------

  insert_msg('Reading the datasets')

  data(fmone)
  data(cds)

# Clinical metadata for the expression dataset -------

  insert_msg('Clinical metadata for the expression dataset')

  ## clearing and restricting only to the bladder and ureter samples

  imvigor$exp_clinic <- cds@phenoData@data %>%
    rownames_to_column('sample_id') %>%
    as_tibble %>%
    transmute(sample_id = sample_id,
              response = ifelse(`Best Confirmed Overall Response` == 'NE',
                                NA,
                                as.character(`Best Confirmed Overall Response`)),
              response = factor(response, c('PD', 'SD', 'PR', 'CR')),
              bi_response = factor(binaryResponse, c('SD/PD', 'CR/PR')),
              ic_enrollment = `Enrollment IC`,
              ic_level = `IC Level`,
              tc_level = `TC Level`,
              immune_pheno = `Immune phenotype`,
              tmb_per_mb = `FMOne mutation burden per MB`,
              sex = car::recode(as.character(Sex),
                                "'M' = 'male'; 'F' = 'female'"),
              sex = factor(sex, c('female', 'male')),
              race = factor(tolower(Race)),
              invasiveness = factor('muscle invasive',
                                    c('non-muscle invasive', 'muscle invasive')),
              bcg_therapy = car::recode(as.character(`Intravesical BCG administered`),
                                        "'N' = 'no'; 'Y' = 'yes'"),
              bcg_therapy = factor(bcg_therapy, c('no', 'yes')),
              intravesical_treatment = bcg_therapy,
              ecog_baseline = `Baseline ECOG Score`,
              smoking = factor(tolower(`Tobacco Use History`),
                               c('never', 'previous', 'current')),
              smoking_history = ifelse(smoking %in% c('previous', 'current'),
                                       'yes', 'no'),
              smoking_history = factor(smoking_history, c('no', 'yes')),
              meta_site = car::recode(as.character(`Met Disease Status`),
                                      "'Liver' = 'liver';
                                      'LN Only' = 'LN';
                                      'Visceral' = 'visceral'"),
              meta_site = factor(meta_site, c('LN', 'visceral', 'liver')),
              sample_age = `Sample age`,
              tissue = Tissue,
              platinum = car::recode(as.character(`Received platinum`),
                                     "'Y' = 'yes'; 'N' = 'no'"),
              platinum = factor(platinum, c('no', 'yes')),
              neo_platinum = car::recode(as.character(`Sample collected pre-platinum`),
                                         "'Y' = 'yes'; 'N' = 'no'"),
              neo_platinum = factor(neo_platinum, c('no', 'yes')),
              neoadjuvant = neo_platinum,
              systemic_treatment = factor('yes', c('no', 'yes')), ## chemo-resistant tumors in the cohort
              systemic_chemotherapy = ifelse(platinum == 'yes' |
                                               neo_platinum == 'yes',
                                             'yes', 'no'),
              systemic_chemotherapy = factor(systemic_chemotherapy,
                                             c('no', 'yes')),
              immunotherapy = factor('yes', c('no', 'yes')), ## as per protocol
              neoag_per_mb = `Neoantigen burden per MB`,
              size_factor = sizeFactor,
              patient_id = ANONPT_ID,
              os_months = os,
              os_days = os_months * 30.437, ## average month length
              death = censOS,
              lund_subtype = Lund,
              lund2_subtype = Lund2,
              tcga_subtype = `TCGA Subtype`)

  ## excluding non-bladder and non-ureter samples

  imvigor$exp_clinic <- imvigor$exp_clinic %>%
    filter(tissue %in% c('bladder', 'ureter')) %>%
    mutate(tissue = ifelse(tissue == 'bladder',
                           tissue, 'non-bladder'),
           tissue = factor(tissue, c('bladder', 'non-bladder')))

# Expression data and its annotation -------

  insert_msg('Expression data and annotation')

  ## raw count table

  imvigor$expression <- cds@assayData$counts %>%
    as.data.frame %>%
    rownames_to_column('entrez_id') %>%
    as_tibble

  ## expression annotation
  ## manual substitution of
  ## the gene '102800317' and the gene '3730'
  ## according to the recent NCBI Gene

  imvigor$exp_annotation <-
    cds@featureData@data[c('entrez_id',
                           'symbol')] %>%
    transmute(entrez_id = entrez_id,
              gene_symbol = symbol) %>%
    as_tibble

  imvigor$exp_annotation[imvigor$exp_annotation$entrez_id == '102800317', 'gene_symbol'] <-
    'TPTEP2-CSNK1E'

  imvigor$exp_annotation[imvigor$exp_annotation$entrez_id == '3730', 'gene_symbol'] <-
    'ANOS1'

  ## re-naming the expression data set,
  ## removing features without gene symbol assignment
  ## (checked: mostly discontinued EntrezID features),
  ## restricting the expression data only to the participants
  ## with bladder or ureter samples
  ## log2 transformation

  imvigor$expression <- imvigor$expression %>%
    mutate(gene_symbol = exchange(entrez_id,
                                  dict = imvigor$exp_annotation,
                                  key = 'entrez_id',
                                  value = 'gene_symbol'),
           gene_symbol = unname(gene_symbol),
           .before = entrez_id) %>%
    filter(!is.na(gene_symbol)) %>%
    filter(gene_symbol != '') %>%
    select( - entrez_id) %>%
    column_to_rownames('gene_symbol') %>%
    t

  imvigor$expression <- imvigor$expression[imvigor$exp_clinic$sample_id, ]

  imvigor$expression <- log2(imvigor$expression + 1)

  imvigor$expression <- imvigor$expression %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

  ## restricting the annotation data frame to the genes present
  ## in the clean expression data set

  imvigor$exp_annotation <- imvigor$exp_annotation %>%
    filter(gene_symbol %in% names(imvigor$expression))

  ## appending with the patient's ID

  imvigor$expression <-
    inner_join(imvigor$exp_clinic[c('sample_id', 'patient_id')],
               imvigor$expression,
               by = 'sample_id')

# Mutation metadata ------

  insert_msg('Mutation metadata')

  ## restricted to individuals with the expression metadata available

  imvigor$mut_clinic <- fmone@phenoData@data %>%
    rownames_to_column('sample_id') %>%
    as_tibble %>%
    transmute(sample_id = sample_id,
              patient_id = as.character(ANONPT_ID)) %>%
    filter(patient_id %in% imvigor$exp_clinic$patient_id) %>%
    left_join(imvigor$exp_clinic[names(imvigor$exp_clinic) != 'sample_id'],
              by = 'patient_id')

# Mutation data --------

  insert_msg('Mutation data')

  ## mutation data recoding. If at least one mutation is found
  ## its coded as 1, WT allele is represented by 0

  imvigor$mutation <- fmone@assayData %>%
    as.list %>%
    map(~ifelse(.x != '', 1, 0)) %>%
    map(t) %>%
    reduce(`+`)

  imvigor$mutation <- ifelse(imvigor$mutation != 0, 1, 0)

  ## restricting to the individuals with clinical metadata available
  ## renaming and annotation data

  imvigor$mutation <- imvigor$mutation[imvigor$mut_clinic$sample_id, ]

  imvigor$mut_annotation <-
    tibble(variable = paste0(colnames(imvigor$mutation), '_mut'),
           label = colnames(imvigor$mutation))

  imvigor$mutation <- imvigor$mutation %>%
    as.data.frame %>%
    rownames_to_column('sample_id') %>%
    as_tibble

  ## appending with the patient ID

  imvigor$mutation <-
    inner_join(imvigor$mut_clinic[c('sample_id', 'patient_id')],
               imvigor$mutation,
               by = 'sample_id')

# Samples with the complete mutation and expression information ------

  insert_msg('Samples with complete mutation and expression data')

  ## and removal of single duplicates

  imvigor$complete_ids <- imvigor[c("expression", "mutation")] %>%
    map(~.x$patient_id) %>%
    reduce(intersect)

  imvigor[c("exp_clinic", "expression", "mut_clinic", "mutation")] <-
    imvigor[c("exp_clinic", "expression", "mut_clinic", "mutation")] %>%
    map(filter,
        patient_id %in% imvigor$complete_ids,
        !duplicated(patient_id)) %>%
    map(mutate,
        sample_id = patient_id)

# Simplifying the clinical information, MIBC bladder samples -------

  insert_msg('Simplifying the clinical information, selection of bladder MIBC')

  imvigor$clinic <- imvigor$mut_clinic

  imvigor$analysis_ids <- imvigor$clinic %>%
    filter(tissue == 'bladder',
           invasiveness == 'muscle invasive') %>%
    .$sample_id

  imvigor[c("clinic", "expression", "mutation")] <-
    imvigor[c("clinic", "expression", "mutation")] %>%
    map(filter, sample_id %in% imvigor$analysis_ids)

# Caching ------

  insert_msg('Caching')

  imvigor <- imvigor[c("clinic",
                       "expression", "exp_annotation",
                       "mutation", "mut_annotation",
                       "analysis_ids")]

  save(imvigor, file = './data/imvigor.RData')

# END ------

  insert_tail()
