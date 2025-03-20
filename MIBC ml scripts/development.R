# Development of and Elastic Net and Random Forest models of consensus
# molecular classes in the TCGA cohort.
# The models use the full set of transcriptional predictors proposed by
# Kamoun et al.

  insert_head()

# container --------

  mibc_train <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  mibc_train$x <- mibc_globals$train_x
  mibc_train$y <- mibc_globals$train_y

  mibc_train$cv_folds <- mibc_globals$cv_folds
  mibc_train$ranger_grid <- mibc_globals$ranger_grid
  mibc_train$gbm_grid <- mibc_globals$gbm_grid

# tuning of the models --------

  insert_msg('Tuning')

  plan('multisession')

  mibc_train$tune_obj$elnet <- tune_glmnet(mibc_train$x,
                                         mibc_train$y,
                                         indexes = mibc_train$cv_folds,
                                         family = 'multinomial',
                                         alpha = 0.5)

  mibc_train$tune_obj$ranger <- tune_ranger(mibc_train$x,
                                          mibc_train$y,
                                          tune_grid = mibc_train$ranger_grid,
                                          tune_stat = 'prediction.error',
                                          num.trees = 500)

  plan('sequential')

# OOF predictions used by the Ensemble model ---------

  insert_msg('OOF predictions')

  mibc_train$oof_predictions$elnet <-
    oof_glmnet(mibc_globals$ens_cv_x$train,
               mibc_globals$ens_cv_x$test,
               mibc_globals$ens_cv_y$train,
               mibc_globals$ens_cv_y$test,
               family = 'multinomial',
               alpha = 0.5,
               lambda = mibc_train$tune_obj$elnet$best_tune$lambda)

  mibc_train$oof_predictions$ranger <-
    call2(.fn = 'oof_ranger',
          mibc_globals$ens_cv_x$train,
          mibc_globals$ens_cv_x$test,
          mibc_globals$ens_cv_y$train,
          mibc_globals$ens_cv_y$test,
          num.trees = 500,
          !!!mibc_train$tune_obj$ranger$best_tune) %>%
    eval

# Caching -------

  insert_msg('Caching')

  mibc_train <- mibc_train[c("tune_obj", "oof_predictions")]

  save(mibc_train, file = './cache/mibc_train.RData')

# END ---------

  insert_tail()
