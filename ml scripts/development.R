# Development of and Elastic Net and Random Forest models of consensus
# molecular classes in the TCGA cohort.

  insert_head()

# container --------

  ml_train <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ml_train$x <- ml_globals$train_x
  ml_train$y <- ml_globals$train_y

  ml_train$cv_folds <- ml_globals$cv_folds
  ml_train$ranger_grid <- ml_globals$ranger_grid

# tuning of the lambda --------

  insert_msg('Tuning')

  plan('multisession')

  ml_train$tune_obj$elnet <- tune_glmnet(ml_train$x,
                                         ml_train$y,
                                         indexes = ml_train$cv_folds,
                                         family = 'multinomial',
                                         alpha = 0.5)

  ml_train$tune_obj$ranger <- tune_ranger(ml_train$x,
                                          ml_train$y,
                                          tune_grid = ml_train$ranger_grid,
                                          tune_stat = 'prediction.error',
                                          num.trees = 500)

  plan('sequential')

# OOF predictions used by the Ensemble model ---------

  insert_msg('OOF predictions')

  ml_train$oof_predictions$elnet <-
    oof_glmnet(ml_globals$ens_cv_x$train,
               ml_globals$ens_cv_x$test,
               ml_globals$ens_cv_y$train,
               ml_globals$ens_cv_y$test,
               family = 'multinomial',
               alpha = 0.5,
               lambda = ml_train$tune_obj$elnet$best_tune$lambda)

  ml_train$oof_predictions$ranger <-
    call2(.fn = 'oof_ranger',
          ml_globals$ens_cv_x$train,
          ml_globals$ens_cv_x$test,
          ml_globals$ens_cv_y$train,
          ml_globals$ens_cv_y$test,
          num.trees = 500,
          !!!ml_train$tune_obj$ranger$best_tune) %>%
    eval

# Caching -------

  insert_msg('Caching')

  ml_train <- ml_train[c("tune_obj", "oof_predictions")]

  save(ml_train, file = './cache/ml_train.RData')

# END ---------

  insert_tail()
