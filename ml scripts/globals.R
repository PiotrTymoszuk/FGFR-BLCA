# Data for ML modeling of molecular subtypes with expression of FGF, FGFR,
# and FGFBP genes. The expression levels are ComBat-processed.
#
# The consensus classes to be modeled are LumP, LumU, LumNS, Stroma-rich,
# and Ba/Sq. Not assigned samples, NE-like and LumNS subset specimens are
# excluded; the later ones because of their low frequency (n = 2 and n = 14).

  insert_head()

# container -------

  ml_globals <- list()

# analysis data ------

  insert_msg('Analysis globals')

  ## variables of interest and their lexicon

  ml_globals$lexicon <-
    globals[c("ligands", "receptors", "binding_proteins")] %>%
    map(~tibble(variable = .x)) %>%
    compress(names_to = 'protein_type')

  ml_globals$data <- combat$expression[c("tcga", "imvigor", "bcan")] %>%
    map(select, sample_id, any_of(ml_globals$lexicon$variable))

  ml_globals$variables <- ml_globals$data %>%
    map(select, -sample_id) %>%
    map(names) %>%
    reduce(intersect)

  ml_globals$lexicon <- ml_globals$lexicon %>%
    filter(variable %in% ml_globals$variable)

  ## consensus class information

  ml_globals$assignment <- subtypes$assignment %>%
    map(~.x[c('sample_id', 'consensusClass')]) %>%
    map(~filter(.x,
                complete.cases(.x),
                !consensusClass %in% c('not assigned', 'NE-like', 'LumNS'))) %>%
    map(mutate,
        consensusClass = droplevels(consensusClass))

  ml_globals$data <-
    map2(ml_globals$data,
         ml_globals$assignment,
         inner_join, by = 'sample_id') %>%
    map(~filter(.x, complete.cases(.x)))

  ml_globals$class_levels <- levels(ml_globals$data[[1]]$consensusClass)

# Response vectors and model matrices ------

  insert_msg('Survival objects and x matrices')

  ml_globals$y <- ml_globals$data %>%
    map(~.x$consensusClass)

  ml_globals$x <- ml_globals$data %>%
    map(column_to_rownames, 'sample_id') %>%
    map(select, all_of(ml_globals$variables)) %>%
    map(as.matrix)

# Augmentation of the training data by SMOTE: expansion of the minorities -------

  insert_msg('Augmentation of the minority groups by SMOTE')

  set.seed(12345)

  ml_globals$train_data <- ml_globals$data$tcga %>%
    column_to_rownames('sample_id') %>%
    augment_multi(minorities = c('LumU', 'Stroma-rich'),
                  perc.over = c(2.7, 2),
                  perc.under = c(1, 1)) %>%
    remove_missing

  ml_globals$train_x <- ml_globals$train_data %>%
    select(all_of(ml_globals$variables)) %>%
    as.matrix

  ml_globals$train_y <- ml_globals$train_data$consensusClass

# CV fold used for tuning of Elastic Net models ---------

  insert_msg('CV folds used for tuning of Elastic Net models')

  ml_globals$n_rep <- 200

  set.seed(1234)

  ml_globals$cv_folds <- 1:ml_globals$n_rep %>%
    map(function(x) createFolds(ml_globals$train_data$consensusClass,
                                list = FALSE,
                                returnTrain = TRUE)) %>%
    set_names(paste0('rep_', 1:ml_globals$n_rep))

# Tune grid for Random Forests ------

  insert_msg('Tune grid for random forests')

  ml_globals$ranger_grid <-
    expand.grid(mtry = seq(2:20),
                min.node.size = c(1, 2, 3),
                splitrule = c('gini', 'extratrees'),
                stringsAsFactors = FALSE)

# CV folds used for construction of the ensembles --------

  insert_msg('CV fold for construction of the ensembles')

  ml_globals$ens_cv_id <-
    createFolds(ml_globals$train_data$consensusClass,
                k = 5,
                list = TRUE,
                returnTrain = TRUE)

  ml_globals$ens_cv_x$train <- ml_globals$ens_cv_id %>%
    map(~ml_globals$train_data[.x, ]) %>%
    map(select, -consensusClass) %>%
    map(as.matrix)

  ml_globals$ens_cv_x$test <- ml_globals$ens_cv_id %>%
    map(~ml_globals$train_data[-.x, ]) %>%
    map(select, -consensusClass) %>%
    map(as.matrix)

  ml_globals$ens_cv_y$train <- ml_globals$ens_cv_id %>%
    map(~ml_globals$train_data[.x, 'consensusClass'])

  ml_globals$ens_cv_y$test <- ml_globals$ens_cv_id %>%
    map(~ml_globals$train_data[-.x, 'consensusClass'])

# END -----

  ml_globals$clinic <- NULL

  ml_globals <- compact(ml_globals)

  insert_tail()
