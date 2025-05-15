# Investigating if and how the consensus molecular subtypes, at least the major
# ones (LumP, LumU, Stroma-rich, and basal/squamous) can be emulated by
# expression of FGF, FGFR, and FGFBP genes.
#
# To this end, we're developing a multinomial Elastic Net and a Random Forest
# model of the consensus classes in the training portion TCGA BLCA cohort
# (2/3 of observations).
# ComBat log2-transformed expression of the
# genes of interest serve as explanatory factors.
# Because of an imbalance of consensus class sizes, we're amplifying the
# minority classes (LumU and Stroma-rich) by SMOTE.
# The models are tuned (optimal lambda and RF hyper-parameters) by
# minimizing deviance in repeated cross-validation and by minimizing
# classification errors of OOB predictions for, respectively, Elastic Net and
# Random Forest models.
# Performance of the model is evaluated in the TCGA (training and test portion),
# IMvigor, and BCAN cohorts, where the consensus cluster assignment predicted
# by consensusMIBC package is assumed as the ground truth.

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(glmnet)
  library(ranger)

  library(caret)
  library(performanceEstimation)
  library(exda)

  library(kernelshap)
  library(shapviz)

  library(furrr)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map
  set_rownames <- trafo::set_rownames
  rename <- dplyr::rename

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis globals --------

  insert_msg('Analysis globals')

  c('./ml scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts -------

  insert_msg('Analysis scripts')

  ## model tuning and training, SHAP importance

  list(cache_path = c('./cache/ml_train.RData',
                      './cache/ml_ens.RData',
                      './cache/ml_shap.RData'),
       script_path = c('./ml scripts/development.R',
                       './ml scripts/ensemble.R',
                       './ml scripts/shap.R'),
       message = c('Cached model tuning and training results',
                   'Cached ensemble tuning and training results',
                   'Cached SHAP variable importance')) %>%
    pwalk(access_cache)

  ## predictions, evaluation, and variable importance

  c('./ml scripts/predictions.R',
    './ml scripts/evaluation.R',
    './ml scripts/importance.R',
    './ml scripts/shap_plots.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END --------

  insert_tail()
