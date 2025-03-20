# Prediction of individual AUC values for bladder cancer samples with
# RIDGE regression and GDSC2 drug screening experiment training data.
#
# The strategy:
#
# 1) selection of compounds and cell lines for the training data set.
# The cell lines are expected to be of epithelial origin, bone, nervous system,
# and blood cancers are excluded.
#
# 2) selection of explanatory factors: 60% of genes with the highest expression
# in the training GDSC data set (i.e. 0.4 quantile split) are selected
# and restricted to genes detected by all studies.
# The rationale: inclusion of genes that
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
# 4) Subsequently, AUC are predicted for the cancer samples based
# on expression of the selected genes (ComBat-adjusted).

  insert_head()

# container -------

  ctrp2 <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# test data: gene expression ------

  insert_msg('Test data: gene expression')

  ctrp2$cancer_genes <- globals$cohort_expr %>%
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
  ## and observations in columns

  ctrp2$x <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(column_to_rownames, 'sample_id')

  ctrp2$x <-
    map2(ctrp2$x,
         ctrp2$cancer_genes[names(ctrp2$x)],
         ~select(.x, any_of(.y))) %>%
    map(t)

# Reading the processed training data ------

  insert_msg('Reading the processed training data')

  list(cache_path = c('./data/CTRP2/ctrp2_exp.RData',
                      './data/CTRP2/ctrp2_drugs.RData'),
       script_path = c('./import scripts/ctrp2_expression.R',
                       './import scripts/ctrp2_drugs.R'),
       message = c('Cached expression for the CTRP2 cell lines',
                   'Cached drug sensitivity for the CTRP2 cell lines')) %>%
    pwalk(access_cache)

# Training data: AUC --------

  insert_msg('Training data, IC50')

  ## In the training, cell lines of epithelial origin are included

  ctrp2$train_tissues <- ctrp2_drugs$tissue_counts %>%
    filter(!tissue %in% c('Ampulla of Vater',
                          'Bone',
                          'CNS/Brain',
                          'Lymphoid',
                          'Myeloid',
                          'Ovary/Fallopian Tube',
                          'Peripheral Nervous System')) %>%
    .$tissue

  ctrp2$train_samples <- ctrp2_drugs$samples %>%
    filter(tissue %in% ctrp2$train_tissues) %>%
    .$sample_id

  ## AUC: neither log nor sqrt conversion help to improve normality
  ## approx. for 100 compounds, the normality assumption is violated
  ## Hence, AUC is used in modeling without any transformation

  ctrp2$train_auc <- ctrp2_drugs$auc %>%
    filter(sample_id %in% ctrp2$train_samples) %>%
    column_to_rownames('sample_id') %>%
    as.matrix

  ctrp2$normality_check <- ctrp2$train_auc %>%
    as.data.frame %>%
    map_dfr(shapiro_test) %>%
    mutate(variable = colnames(ctrp2$train_auc))

# Training data: expression --------

  insert_msg('Training data: expression')

  ## conversion to a matrix format with genes as rows
  ## and samples as columns

  ctrp2$x$train <- ctrp2_exp$expression %>%
    filter(sample_id %in% ctrp2$train_samples) %>%
    column_to_rownames('sample_id') %>%
    t

# Selection of modeling genes -------

  insert_msg('Selection of explanatory factors')

  ## identification of genes:
  ##
  ## by analysis of histogram of median expression per gene in the
  ## training data set with expression values for cell lines,
  ## roughly 40% of investigated genes are not expressed at all
  ## and hence filtered out. The upper half is selected for modeling

  ctrp2$train_cutoff <- ctrp2$x$train %>%
    rowMedians %>%
    quantile(0.4)

  ctrp2$train_genes <- ctrp2$x$train %>%
    row_stats %>%
    filter(median >= ctrp2$train_cutoff) %>%
    .$variable

  ctrp2$common_genes <- ctrp2$x %>%
    map(rownames) %>%
    reduce(intersect)

  ctrp2$train_genes <-
    ctrp2$train_genes[ctrp2$train_genes %in% ctrp2$common_genes]

  ## batch effect removal for the common genes selected above

  ctrp2$x <- ctrp2$x %>%
    map(~.x[ctrp2$train_genes, ])

  ctrp2$data_object <-
    multi_process(train = ctrp2$x$train,
                  test = ctrp2$x[names(ctrp2$x) != 'train'])

# Gene expression distribution stats -------

  insert_msg('Gene expression distribution stats')

  ctrp2$gene_stats <- components(ctrp2$data_object, type = 'stats')

  ctrp2$gene_stats$before_adjustment <-
    ctrp2$gene_stats$before_adjustment %>%
    filter(modeling_variable == 'yes') %>%
    select(-modeling_variable)

  ctrp2$gene_stats <- ctrp2$gene_stats %>%
    set_names(c('raw', 'adjusted')) %>%
    compress(names_to = 'combat')

  ctrp2$gene_stats <- ctrp2$gene_stats %>%
    mutate(dataset = ifelse(data_set == 'train',
                            'training', 'test'),
           data_set = factor(data_set,
                             c('train', globals$analysis_cohorts)),
           combat = factor(combat, c('raw', 'adjusted')))

  plan('sequential')

# QC plots for the pre-processing step --------

  insert_msg('Pre-processing plots')

  ctrp2$gene_plot <- ctrp2$gene_stats %>%
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
                                'train' = 'CTRP2')) +
    guides(x = guide_axis(angle = 45)) +
    globals$common_theme +
    theme(axis.title.x = element_blank()) +
    labs(title = 'Median expression, modeling variables',
         subtitle = paste('n = ', nobs(ctrp2$data_object)$n[19], 'genes'),
         y = 'median expression')

# Clearing of the data space ------

  insert_msg('Cleaning the data space')

  ## cleaning and deactivating the parallel backend
  ## to space memory

  ctrp2 <-
    ctrp2[c("x", "data_object",
            "train_drugs", "train_samples", "train_auc",
            "normality_check", "gene_stats", "gene_plot")]

  plan('sequential')

# Training and evaluation of the RIDGE models ---------

  insert_msg('Model selection, training and evaluation')

  plan('multisession')

  ctrp2$train_object <-
    train.default(x = ctrp2$data_object$train,
                  y = ctrp2$train_auc,
                  standardize = TRUE,
                  alpha = 0,
                  family = 'gaussian',
                  type.measure = 'mse')

  plan('sequential')

# Model evaluation stats --------

  insert_msg('Model evaluation stats')

  ## quality criteria for the models:
  ## cross-validated R^2 of at least 0.13.
  ## correlation of the predicted and observed response with r >= 0.5

  ctrp2$model_stats <- ctrp2$train_object %>%
    summary %>%
    mutate(qc_passed = ifelse(rsq_oof >= 0.13 & pearson >= 0.5,
                              'yes', 'no'))

  ## plotting of the cross-validated R^2 and Pearson's correlation

  ctrp2$model_plot <- ctrp2$train_object %>%
    plot

  ctrp2$model_plot <- ctrp2$model_plot$rsq_spearman +
    geom_vline(xintercept = 0.13,
               linetype = 'dashed') +
    geom_hline(yintercept = 0.5,
               linetype = 'dashed') +
    globals$common_theme +
    labs(title = 'Model performance, CV and training',
         subtitle = paste0('total: n = ', nrow(ctrp2$model_stats),
                           ', QC passed: n = ', table(ctrp2$model_stats$qc_passed)['yes']))

# Compound lexicon ---------

  insert_msg('Compound lexicon for the QC-passing models')

  ctrp2$lexicon <- ctrp2$model_stats %>%
    filter(qc_passed == 'yes') %>%
    transmute(variable = response) %>%
    left_join(., ctrp2_drugs$lexicon,
              by = 'variable')

# Predictions -------

  insert_msg('Predictions of AUC')

  ctrp2$predictions <-
    predict(ctrp2$train_object,
            newdata = ctrp2$data_object,
            type = 'response')

  ## formatting of the predictions:
  ## stored as a tibble

  ctrp2$predictions <- ctrp2$predictions %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Caching the results --------

  insert_msg('Caching the results')

  ctrp2 <-
    ctrp2[c("normality_check",
            "gene_stats", "gene_plot",
            "train_object", "model_stats", "model_plot",
            "lexicon", "predictions")]

  save(ctrp2, file = './data/ctrp2.RData')

# END ------

  rm(ctrp2_drugs, ctrp2_exp)

  insert_tail()
