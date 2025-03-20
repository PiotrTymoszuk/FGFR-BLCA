# Data for ML modeling of molecular subtypes with expression of FGF, FGFR,
# and FGFBP genes. The expression levels are ComBat-processed.
#
# The consensus classes to be modeled are LumP, LumU, LumNS, Stroma-rich,
# and Ba/Sq. Not assigned samples, NE-like and LumNS subset specimens are
# excluded; the later ones because of their low frequency (n = 2 and n = 14).

  insert_head()

# container -------

  mibc_globals <- list()

# analysis data ------

  insert_msg('Analysis globals')

  ## variables of interest and their lexicon

  mibc_globals$lexicon <-
    tibble(variable = consensusMIBC::centroids$hgnc_symbol) %>%
    left_join(ml_globals$lexicon, by = 'variable') %>%
    mutate(is.na(protein_type), 'other')

  ## expression data

  mibc_globals$data <- combat$expression[c("tcga", "imvigor", "bcan")] %>%
    map(select, sample_id, any_of(mibc_globals$lexicon$variable))

  mibc_globals$variables <- mibc_globals$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  mibc_globals$lexicon <- mibc_globals$lexicon %>%
    filter(variable %in% mibc_globals$variable)

  ## consensus class information

  mibc_globals$assignment <- subtypes$assignment %>%
    map(~.x[c('sample_id', 'consensusClass')]) %>%
    map(~filter(.x,
                complete.cases(.x),
                !consensusClass %in% c('not assigned', 'NE-like', 'LumNS'))) %>%
    map(mutate,
        consensusClass = droplevels(consensusClass))

  mibc_globals$data <-
    map2(mibc_globals$data,
         mibc_globals$assignment,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

  mibc_globals$class_levels <- levels(mibc_globals$data[[1]]$consensusClass)

# Definition of the training and test subset of the TCGA cohort -------

  insert_msg('Definition of the training and test subset of the TCGA cohort')

  set.seed(1234)

  mibc_globals$tcga_train_indexes <-
    createDataPartition(mibc_globals$data$tcga$consensusClass,
                        p = 2/3)[[1]]

  mibc_globals$data$tcga_train <-
    mibc_globals$data$tcga[mibc_globals$tcga_train_indexes, ]

  mibc_globals$data$tcga_test <-
    mibc_globals$data$tcga[-mibc_globals$tcga_train_indexes, ]

# Response vectors and model matrices ------

  insert_msg('Survival objects and x matrices')

  mibc_globals$y <- mibc_globals$data %>%
    map(~.x$consensusClass)

  mibc_globals$x <- mibc_globals$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(mibc_globals$variables)) %>%
    map(as.matrix)

# Augmentation of the training data by SMOTE: expansion of the minorities -------

  insert_msg('Augmentation of the minority groups by SMOTE')

  set.seed(12345)

  mibc_globals$train_data <- mibc_globals$data$tcga_train %>%
    column_to_rownames('sample_id') %>%
    augment_multi(minorities = c('LumU', 'Stroma-rich'),
                  perc.over = c(2.7, 2),
                  perc.under = c(1, 1)) %>%
    remove_missing

  mibc_globals$train_x <- mibc_globals$train_data %>%
    select(all_of(mibc_globals$variables)) %>%
    as.matrix

  mibc_globals$train_y <- mibc_globals$train_data$consensusClass

# CV fold used for tuning of Elastic Net models ---------

  insert_msg('CV folds used for tuning of Elastic Net models')

  mibc_globals$n_rep <- 200

  set.seed(1234)

  mibc_globals$cv_folds <- 1:mibc_globals$n_rep %>%
    map(function(x) createFolds(mibc_globals$train_data$consensusClass,
                                list = FALSE,
                                returnTrain = TRUE)) %>%
    set_names(paste0('rep_', 1:mibc_globals$n_rep))

# Tune grid for Random Forests ------

  insert_msg('Tune grid for random forests')

  mibc_globals$ranger_grid <-
    expand.grid(mtry = seq(2:20),
                min.node.size = c(1, 2, 3),
                splitrule = c('gini', 'extratrees'),
                stringsAsFactors = FALSE)

# CV folds used for construction of the ensembles --------

  insert_msg('CV fold for construction of the ensembles')

  mibc_globals$ens_cv_id <-
    createFolds(mibc_globals$train_data$consensusClass,
                k = 5,
                list = TRUE,
                returnTrain = TRUE)

  mibc_globals$ens_cv_x$train <- mibc_globals$ens_cv_id %>%
    map(~mibc_globals$train_data[.x, ]) %>%
    map(select, -consensusClass) %>%
    map(as.matrix)

  mibc_globals$ens_cv_x$test <- mibc_globals$ens_cv_id %>%
    map(~mibc_globals$train_data[-.x, ]) %>%
    map(select, -consensusClass) %>%
    map(as.matrix)

  mibc_globals$ens_cv_y$train <- mibc_globals$ens_cv_id %>%
    map(~mibc_globals$train_data[.x, 'consensusClass'])

  mibc_globals$ens_cv_y$test <- mibc_globals$ens_cv_id %>%
    map(~mibc_globals$train_data[-.x, 'consensusClass'])

# END -----

  mibc_globals$clinic <- NULL

  mibc_globals <- compact(mibc_globals)

  insert_tail()
