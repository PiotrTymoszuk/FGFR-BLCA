# Ensemble model that operates with predictions of the consensus classes
# made by Elastic Net and Random Forest algorithms.

  insert_head()

# container --------

  mibc_ens <- list()

# modeling data --------

  insert_msg('Modeling data')

  ## class levels

  mibc_ens$class_levels <- mibc_globals$class_levels

  ## out-of-fold predictions for the SMOTE data

  mibc_ens$data <-
    list(elnet = mibc_train$oof_predictions$elnet %>%
           select(-obs, -pred),
         ranger = mibc_train$oof_predictions$ranger) %>%
    reduce(inner_join, by = c('sample_id', 'fold_id'))

  ## modeling matrix and response vector

  mibc_ens$x <- mibc_ens$data %>%
    column_to_rownames('sample_id') %>%
    select(ends_with('elnet'), ends_with('ranger')) %>%
    as.matrix

  mibc_ens$y <- mibc_ens$data$obs

  ## CV folds for tuning of the lambda

  mibc_ens$n_rep <- 200

  set.seed(1234)

  mibc_ens$cv_folds <- 1:mibc_ens$n_rep %>%
    map(function(x) createFolds(mibc_ens$data$obs,
                                list = FALSE,
                                returnTrain = TRUE)) %>%
    set_names(paste0('rep_', 1:mibc_ens$n_rep))

# Elastic Net ensembles: tuning and training in the TCGA cohort ---------

  insert_msg('Tuning and training')

  plan('multisession')

  mibc_ens$tune_obj <-
    tune_glmnet(mibc_ens$x,
                mibc_ens$y,
                indexes = mibc_ens$cv_folds,
                family = 'multinomial',
                alpha = 0.5)

  plan('sequential')

# Caching -------

  insert_msg('Caching')

  mibc_ens <- mibc_ens["tune_obj"]

  save(mibc_ens, file = './cache/mibc_ens.RData')

# END ---------

  insert_tail()
