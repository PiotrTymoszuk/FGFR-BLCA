# Development of and Elastic Net model of overall survival in the TCBA cohort.

  insert_head()

# container --------

  surv_train <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  surv_train$x <- surv_globals$x$tcga
  surv_train$y <- surv_globals$y$tcga

  surv_train$cv_folds <- surv_globals$cv_folds

# tuning of the lambda --------

  insert_msg('Tuning')

  plan('multisession')

  surv_train$tune_obj <- tune_glmnet(surv_train$x,
                                     surv_train$y,
                                     indexes = surv_train$cv_folds,
                                     family = 'cox',
                                     alpha = 0.5)

  plan('sequential')

# Caching -------

  insert_msg('Caching')

  surv_train <- surv_train["tune_obj"]

  save(surv_train, file = './cache/surv_train.RData')

# END -------

  insert_tail()
