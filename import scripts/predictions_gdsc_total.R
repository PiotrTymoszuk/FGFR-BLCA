# Prediction of individual IC50 values for bladder cancer samples with
# RIDGE regression and the total GDSC drug screening experiment training data.
#
# The strategy:
#
# 1) selection of compounds and cell lines for the training data set.
# The cell lines are expected to be of epithelial origin, bone, nervous system,
# and blood cancers are excluded. In case of compounds shared by both the GDSC1
# and GDSC2 experiments, the GDSC2 data are preferred, as recommended by the
# project authors.
#
# 2) selection of explanatory factors: 50% of genes with the highest expression
# in the training GDSC data set (i.e. median split) are selected and restricted
# to genes detected by all studies. The rationale: inclusion of genes that
# are abundantly expressed by cancer cells in the culture without any
# non-malignant components of the tumor microenvironment.
# Next, numeric matrices with expression of the selected modeling genes are
# subjected to batch adjustment with ComBat.
#
# 3) Training of RIDGE models. Lambdas are selected with the minimal
# mean squared errors in cross-validation. Successful models are defined
# by cross-validated R^2 >= 0.13 and correlation of the predicted and observed
# outcome with r >= 0.5 (Cohen 1988 and Cohen 2013).
#
# 4) Subsequently, log IC50 [µM] are predicted for the cancer samples based
# on expression of the selected genes (ComBat-adjusted).

  insert_head()

# container -------

  gdsct <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# test data: gene expression ------

  insert_msg('Test data: gene expression')

  gdsct$cancer_genes <- globals$cohort_expr %>%
    eval %>%
    map(~.x[c('exp_annotation', 'annotation')]) %>%
    map(compact) %>%
    map(function(x) if(length(x) == 0) NULL else x) %>%
    compact %>%
    map(~.x[[1]]) %>%
    map(filter, !duplicated(gene_symbol)) %>%
    map(~.x$gene_symbol) %>%
    map(unname)

  ## conversion to numeric matrices with genes in rows
  ## and observations in columns, data are already log2-transformed

  gdsct$x <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(column_to_rownames, 'sample_id')

  gdsct$x <-
    map2(gdsct$x,
         gdsct$cancer_genes[names(gdsct$x)],
         ~select(.x, any_of(.y))) %>%
    map(t)

# Training data: reading the processed data -------

  insert_msg('Reading the processed training data')

  list(cache_path = c('./data/GDSC/gdsc_exp.RData',
                      './data/GDSC/gdsc_drugs.RData'),
       script_path = c('./import scripts/gdsc_expression.R',
                       './import scripts/gdsc_drugs.R'),
       message = c('Cached expression for the GDSC cell lines',
                   'Cached drug sensitivity for the GDSC cell lines')) %>%
    pwalk(access_cache)

# Training data: IC50 ----------

  insert_msg('Training data: IC50')

  ## Goal: matrices with cell IDs in row names
  ## and drugs in column names.
  ##
  ## in case of drugs shared by both GDSC1 and GDSC2 experiments,
  ## GDSC2 is preferred, as recommended by the study authors

  gdsct$train_drugs <- gdsc_drugs$lexicons

  gdsct$train_drugs$gdsc1 <- gdsct$train_drugs$gdsc1 %>%
    filter(!variable %in% gdsc_drugs$common_drugs)

  gdsct$train_drugs <- gdsct$train_drugs %>%
    map(~.x$variable)

  ## Only solid, epithelial origin cell lines are included
  ## in the training data sets

  gdsct$train_tissues <- gdsc_drugs$tissue_counts %>%
    filter(complete.cases(.),
           !tissue %in% c('blood', 'bone', 'nervous_system')) %>%
    .$tissue

  gdsct$train_subtypes <- gdsc_drugs$tissue_subtype_counts %>%
    filter(complete.cases(.),
           tissue_subtype %in% c('adrenal_gland',
                                 'biliary_tract',
                                 'bladder',
                                 'breast',
                                 'cervix',
                                 'digestive_system_other',
                                 'head_and_neck',
                                 'kidney',
                                 'large_intestine',
                                 'liver',
                                 'lung_NSCLC_adenocarcinoma',
                                 'lung_NSCLC_carcinoid',
                                 'lung_NSCLC_large_cell',
                                 'lung_NSCLC_not_specified',
                                 'lung_NSCLC_squamous_cell_carcinoma',
                                 'lung_other',
                                 'lung_small_cell_carcinoma',
                                 'melanoma',
                                 'oesophagus',
                                 'pancreas',
                                 'prostate',
                                 'skin_other',
                                 'soft_tissue_other',
                                 'stomach',
                                 'thyroid',
                                 'urogenital_system_other',
                                 'uterus')) %>%
    .$tissue_subtype

  gdsct$train_samples <- gdsc_drugs$design %>%
    map(filter,
        !duplicated(sample_id),
        tissue %in% gdsct$train_tissues,
        tissue_subtype %in% gdsct$train_subtypes) %>%
    map(~.x$sample_id)

  ## training matrices of log IC50

  gdsct$train_ic50 <- gdsc_drugs$log_ic50 %>%
    map(column_to_rownames, 'sample_id') %>%
    map(as.matrix)

  gdsct$train_ic50 <-
    list(x = gdsct$train_ic50,
         y = gdsct$train_samples,
         z = gdsct$train_drugs) %>%
    pmap(function(x, y, z) x[y, z])

  ## normality check: the assumption is, in spite of the log transformation,
  ## violated for some compounds (10 - 24 drugs with W < 0.9)

  gdsct$normality_check <- gdsct$train_ic50 %>%
    map(as.data.frame) %>%
    map(map_dfr, shapiro_test)

  gdsct$normality_check <-
    map2(gdsct$normality_check,
         gdsct$train_ic50,
         ~mutate(.x, variable = colnames(.y)))

  ## conversion to pM to avoid negative values upon log transformation

  gdsct$train_ic50 <- gdsct$train_ic50 %>%
    map(~log(exp(.x) * 1e6))

# Training data: gene expression --------

  insert_msg('Training data: gene expression')

  ## expression data are already log2 transformed.
  ## matrix format: genes in rows and samples in columns

  gdsct$x$train <- gdsc_exp$expression %>%
    column_to_rownames('sample_id') %>%
    t

# Selection of modeling genes -------

  insert_msg('Selection of explanatory factors')

  ## identification of genes:
  ##
  ## by analysis of histogram of median expression per gene in the
  ## training data set with expression values for cell lines,
  ## roughly 50% of investigated genes are not expressed at all
  ## and hence filtered out. The upper half is selected for modeling

  gdsct$train_median <- gdsct$x$train %>%
    rowMedians %>%
    median

  gdsct$train_genes <- gdsct$x$train %>%
    row_stats %>%
    filter(median >= gdsct$train_median) %>%
    .$variable

  gdsct$common_genes <- gdsct$x %>%
    map(rownames) %>%
    reduce(intersect)

  gdsct$train_genes <-
    gdsct$train_genes[gdsct$train_genes %in% gdsct$common_genes]

  ## batch effect removal for the common genes selected above

  gdsct$x <- gdsct$x %>%
    map(~.x[gdsct$train_genes, ])

  gdsct$data_object <-
    multi_process(train = gdsct$x$train,
                  test = gdsct$x[names(gdsct$x) != 'train'])

# Gene expression distribution stats -------

  insert_msg('Gene expression distribution stats')

  gdsct$gene_stats <- components(gdsct$data_object, type = 'stats')

  gdsct$gene_stats$before_adjustment <-
    gdsct$gene_stats$before_adjustment %>%
    filter(modeling_variable == 'yes') %>%
    select(-modeling_variable)

  gdsct$gene_stats <- gdsct$gene_stats %>%
    set_names(c('raw', 'adjusted')) %>%
    compress(names_to = 'combat')

  gdsct$gene_stats <- gdsct$gene_stats %>%
    mutate(dataset = ifelse(data_set == 'train',
                            'training', 'test'),
           data_set = factor(data_set,
                             c('train', globals$analysis_cohorts)),
           combat = factor(combat, c('raw', 'adjusted')))

  plan('sequential')

# QC plots for the pre-processing step --------

  insert_msg('Pre-processing plots')

  gdsct$gene_plot <- gdsct$gene_stats %>%
    select(dataset, combat, data_set, median) %>%
    ggplot(aes(x = data_set,
               y = median,
               fill = dataset)) +
    geom_boxplot(outlier.alpha = 0.25) +
    facet_grid(combat ~ .,
               scales = 'free') +
    scale_fill_manual(values = c(training = 'indianred3',
                                 test = 'steelblue'),
                      name = 'Data set') +
    scale_x_discrete(labels = c(globals$cohort_labs,
                                'train' = 'GDSC1')) +
    guides(x = guide_axis(angle = 45)) +
    globals$common_theme +
    theme(axis.title.x = element_blank()) +
    labs(title = 'Median expression, modeling variables',
         subtitle = paste('n = ', nobs(gdsct$data_object)$n[19], 'genes'),
         y = expression('log'[2] * ' median expression'))

# Cleaning the data space prior to modeling ---------

  insert_msg('Cleaning the data space')

  ## cleaning and deactivating the parallel backend
  ## to space memory

  gdsct <-
    gdsct[c("x", "data_object",
            "train_drugs", "train_samples", "train_ic50",
            "normality_check", "gene_stats", "gene_plot")]

  plan('sequential')

# Training and evaluation of the RIDGE models ---------

  insert_msg('Model selection, training and evaluation')

  for(i in names(gdsct$train_ic50)) {

    plan('multisession')

    gdsct$train_objects[[i]] <-
      train.default(x = gdsct$data_object$train,
                    y = gdsct$train_ic50[[i]],
                    standardize = TRUE,
                    alpha = 0,
                    family = 'gaussian',
                    type.measure = 'mse')

    plan('sequential')

  }

# Model evaluation stats --------

  insert_msg('Model evaluation stats')

  ## quality criteria for the models:
  ## cross-validated R^2 of at least 0.13.
  ## correlation of the predicted and observed response with r >= 0.5

  gdsct$model_stats <- gdsct$train_objects %>%
    map(summary) %>%
    map(mutate,
        qc_passed = ifelse(rsq_oof >= 0.13 & pearson >= 0.5,
                           'yes', 'no'))

  ## plotting of the cross-validated R^2 and Pearson's correlation

  gdsct$model_plots <- gdsct$train_objects %>%
    map(plot) %>%
    map(~.x$rsq_spearman) %>%
    map(~.x +
          geom_vline(xintercept = 0.13,
                     linetype = 'dashed') +
          geom_hline(yintercept = 0.5,
                     linetype = 'dashed') +
          globals$common_theme)

  ## additional styling

  gdsct$model_plots <-
    list(x = gdsct$model_plots,
         y = c('Model performance, GDSC1',
               'Model performance, GDSC2'),
         z = gdsct$model_stats) %>%
    pmap(function(x, y, z) x  +
           labs(title = y,
                subtitle = paste0('total: n = ', nrow(z),
                                  ', QC passed: n = ', table(z$qc_passed)['yes'])))

# Compound lexicons ---------

  insert_msg('Compound lexicon for the QC-passing models')

  gdsct$lexicons <- gdsct$model_stats %>%
    map(filter, qc_passed == 'yes') %>%
    map(transmute,
        variable = response) %>%
    map2(., gdsc_drugs$lexicons,
         left_join, by = 'variable')

# Predictions -------

  insert_msg('Predictions of IC50')

  for(i in names(gdsct$train_objects)) {

    gdsct$predictions[[i]] <-
      predict(gdsct$train_objects[[i]],
              newdata = gdsct$data_object,
              type = 'response')

    ## formatting of the predictions:
    ## they are converted to log(IC50 [µM])
    ## stored as a tibble

    gdsct$predictions[[i]] <- gdsct$predictions[[i]] %>%
      map(~exp(.x) * 1e-6) %>%
      map(log)

    gdsct$predictions[[i]] <- gdsct$predictions[[i]] %>%
      map(as.data.frame) %>%
      map(rownames_to_column, 'sample_id') %>%
      map(as_tibble)

  }

# Caching the results --------

  insert_msg('Caching the results')

  gdsct <-
    gdsct[c("normality_check",
            "gene_stats", "gene_plot",
            "train_objects", "model_stats", "model_plots",
            "lexicons", "predictions")]

  save(gdsct, file = './data/gdsct.RData')

# END ------

  rm(gdsc_drugs, gdsc_exp)

  rm(i, common_samples)

  insert_tail()
