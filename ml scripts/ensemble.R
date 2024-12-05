# Ensemble model that operates with predictions of the consensus classes
# made by Elastic Net and Random Forest algorithms

  insert_head()

# container --------

  ml_ens <- list()

# modeling data --------

  insert_msg('Modeling data')

  ## class levels

  ml_ens$class_levels <- ml_globals$class_levels

  ## out-of-fold predictions for the SMOTE data

  ml_ens$data <-
    list(elnet = ml_train$oof_predictions$elnet %>%
           select(-obs, -pred),
         ranger = ml_train$oof_predictions$ranger) %>%
    reduce(inner_join, by = c('sample_id', 'fold_id'))

  ## modeling matrix and response vector

  ml_ens$x <- ml_ens$data %>%
    column_to_rownames('sample_id') %>%
    select(ends_with('elnet'), ends_with('ranger')) %>%
    as.matrix

  ml_ens$y <- ml_ens$data$obs

  ## CV folds for tuning of the lambda

  ml_ens$n_rep <- 200

  set.seed(1234)

  ml_ens$cv_folds <- 1:ml_ens$n_rep %>%
    map(function(x) createFolds(ml_ens$data$obs,
                                list = FALSE,
                                returnTrain = TRUE)) %>%
    set_names(paste0('rep_', 1:ml_ens$n_rep))

# Elastic Net ensembles: tuning and training in the TCGA cohort ---------

  insert_msg('Tuning and training')

  plan('multisession')

  ml_ens$tune_obj <-
    tune_glmnet(ml_ens$x,
                ml_ens$y,
                indexes = ml_ens$cv_folds,
                family = 'multinomial',
                alpha = 0.5)

  plan('sequential')

# Caching -------

  insert_msg('Caching')

  ml_ens <- ml_ens["tune_obj"]

  save(ml_ens, file = './cache/ml_ens.RData')

# END ---------

  insert_tail()
