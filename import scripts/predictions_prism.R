# Prediction of individual resistance scores for pancreatic cancer samples with
# RIDGE regression and GDSC2 drug screening experiment training data.
#
# The strategy:
#
# 1) selection of compounds and cell lines for the training data set.
# The cell lines are expected to be of epithelial origin, bone, nervous system,
# and blood cancers are excluded. To improve normality and suppress the two-peak
# characteristic of the predicted response, the AUC are min-max normalized
# and, subsequently, transformed with asin(sqrt(x)) function.
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
# 4) Subsequently, resistance scores are predicted for the cancer samples based
# on expression of the selected genes (ComBat-adjusted).

# container -------

  prism <- list()

# parallel backend -------

  insert_msg('Parallel backend')

  plan('multisession')

# test data: gene expression ------

  insert_msg('Test data: gene expression')

  prism$cancer_genes <- globals$cohort_expr %>%
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

  prism$x <- globals$cohort_expr %>%
    eval %>%
    map(~.x$expression) %>%
    compact %>%
    map(column_to_rownames, 'sample_id')

  prism$x <-
    map2(prism$x,
         prism$cancer_genes[names(prism$x)],
         ~select(.x, any_of(.y))) %>%
    map(t)

# Reading the processed training data ------

  insert_msg('Reading the processed training data')

  list(cache_path = c('./data/PRISM/prism_exp.RData',
                      './data/PRISM/prism_drugs.RData'),
       script_path = c('./import scripts/prism_expression.R',
                       './import scripts/prism_drugs.R'),
       message = c('Cached expression for the PRISM cell lines',
                   'Cached drug sensitivity for the PRISM cell lines')) %>%
    pwalk(access_cache)

# Training data: AUC --------

  insert_msg('Training data, IC50')

  ## In the training, cell lines of epithelial origin are included

  prism$train_tissues <- prism_drugs$tissue_counts %>%
    filter(!tissue %in% c('Ampulla of Vater',
                          'Bone',
                          'CNS/Brain',
                          'Lymphoid',
                          'Myeloid',
                          'Ovary/Fallopian Tube',
                          'Peripheral Nervous System')) %>%
    .$tissue

  prism$train_samples <- prism_drugs$samples %>%
    filter(tissue %in% prism$train_tissues) %>%
    .$sample_id

  ## AUC: normality is not really improved by log or sqrt transformations
  ## resorting to the asin-transformation, which at least reduces
  ## the problem with two distribution peaks
  ## Normality is substantially violated (W < 0.9)
  ## for 842 substances out of ca. 1480 substances

  prism$train_auc <- prism_drugs$auc %>%
    filter(sample_id %in% prism$train_samples) %>%
    column_to_rownames('sample_id') %>%
    min_max %>%
    as.matrix

  prism$train_auc <- asin(sqrt(prism$train_auc))

  prism$normality_check <- prism$train_auc %>%
    as.data.frame %>%
    map_dfr(shapiro_test) %>%
    mutate(variable = colnames(prism$train_auc))

# Training data: expression --------

  insert_msg('Training data: expression')

  ## conversion to a matrix format with genes as rows
  ## and samples as columns

  prism$x$train <- prism_exp$expression %>%
    filter(sample_id %in% prism$train_samples) %>%
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

  prism$train_cutoff <- prism$x$train %>%
    rowMedians %>%
    quantile(0.4)

  prism$train_genes <- prism$x$train %>%
    row_stats %>%
    filter(median >= prism$train_cutoff) %>%
    .$variable

  prism$common_genes <- prism$x %>%
    map(rownames) %>%
    reduce(intersect)

  prism$train_genes <-
    prism$train_genes[prism$train_genes %in% prism$common_genes]

  ## batch effect removal for the common genes selected above

  prism$x <- prism$x %>%
    map(~.x[prism$train_genes, ])

  prism$data_object <-
    multi_process(train = prism$x$train,
                  test = prism$x[names(prism$x) != 'train'])

# Gene expression distribution stats -------

  insert_msg('Gene expression distribution stats')

  prism$gene_stats <- components(prism$data_object, type = 'stats')

  prism$gene_stats$before_adjustment <-
    prism$gene_stats$before_adjustment %>%
    filter(modeling_variable == 'yes') %>%
    select(-modeling_variable)

  prism$gene_stats <- prism$gene_stats %>%
    set_names(c('raw', 'adjusted')) %>%
    compress(names_to = 'combat')

  prism$gene_stats <- prism$gene_stats %>%
    mutate(dataset = ifelse(data_set == 'train',
                            'training', 'test'),
           data_set = factor(data_set,
                             c('train', globals$analysis_cohorts)),
           combat = factor(combat, c('raw', 'adjusted')))

  plan('sequential')

# QC plots for the pre-processing step --------

  insert_msg('Pre-processing plots')

  prism$gene_plot <- prism$gene_stats %>%
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
                                'train' = 'PRISM')) +
    guides(x = guide_axis(angle = 45)) +
    globals$common_theme +
    theme(axis.title.x = element_blank()) +
    labs(title = 'Median expression, modeling variables',
         subtitle = paste('n = ', nobs(prism$data_object)$n[19], 'genes'),
         y = 'median expression')

# Clearing of the data space ------

  insert_msg('Cleaning the data space')

  ## cleaning and deactivating the parallel backend
  ## to space memory

  prism <-
    prism[c("x", "data_object",
            "train_drugs", "train_samples", "train_auc",
            "normality_check", "gene_stats", "gene_plot")]

  plan('sequential')

# Tuning, training, evaluation ---------

  insert_msg('Model selection, training and evaluation')

  plan('multisession')

  prism$train_object <-
    train.default(x = prism$data_object$train,
                  y = prism$train_auc,
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

  prism$model_stats <- prism$train_object %>%
    summary %>%
    mutate(qc_passed = ifelse(rsq_oof >= 0.13 & pearson >= 0.5,
                              'yes', 'no'))

  ## plotting of the cross-validated R^2 and Pearson's correlation

  prism$model_plot <- prism$train_object %>%
    plot

  prism$model_plot <- prism$model_plot$rsq_spearman +
    geom_vline(xintercept = 0.13,
               linetype = 'dashed') +
    geom_hline(yintercept = 0.5,
               linetype = 'dashed') +
    globals$common_theme +
    labs(title = 'Model performance, CV and training',
         subtitle = paste0('total: n = ',
                           nrow(prism$model_stats),
                           ', QC passed: n = ',
                           table(prism$model_stats$qc_passed)['yes']))

# Compound lexicon ---------

  insert_msg('Compound lexicon for the QC-passing models')

  prism$lexicon <- prism$model_stats %>%
    filter(qc_passed == 'yes') %>%
    transmute(variable = response) %>%
    left_join(.,
              prism_drugs$lexicon,
              by = 'variable')

# Predictions -------

  insert_msg('Predictions of AUC')

  prism$predictions <-
    predict(prism$train_object,
            newdata = prism$data_object,
            type = 'response')

  ## formatting of the predictions:
  ## back to
  ## stored as a tibble

  prism$predictions <- prism$predictions %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'sample_id') %>%
    map(as_tibble)

# Model coefficients for erdafitinib, used in later analyses -------

  insert_msg('Selected coefficients')

  prism$coefs <- prism$train_object %>%
    coef

  prism$coefs <- prism$coefs[, 'erdafitinib'] %>%
    compress(names_to = 'gene_symbol',
             values_to = 'beta') %>%
    filter(beta != 0)

# Caching the results --------

  insert_msg('Caching the results')

  prism <-
    prism[c("normality_check",
            "gene_stats", "gene_plot",
            "train_object", "model_stats", "model_plot",
            "lexicon", "predictions", "coefs")]

  save(prism, file = './data/prism.RData')

# END ------

  rm(prism_drugs, prism_exp)

  insert_tail()
