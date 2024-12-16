# Investigating if and how the consensus molecular subtypes, at least the major
# ones (LumP, LumU, Stroma-rich, and basal/squamous) can be emulated by
# expression of FGF, FGFR, and FGFBP genes.
#
# To this end, we're developing a multinomial Elastic Net model of the consensus
# classes in the TCGA BLCA cohort; ComBat log2-transformed expression of the
# genes of interest serve as explanatory factors. Because of an imbalance of
# consensus class sizes, we're amplifying the minority classes
# (LumU and Stroma-rich) by SMOTE. The model is tuned (optimal
# lambda) by minimizing deviance in repeated cross-validation. Performance of
# the model is evaluated in the TCGA, IMvigor, and BCAN cohorts, where
# the consensus cluster assignment predicted by consensusMIBC package is
# assumed as the ground truth.

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(glmnet)
  library(caret)
  library(ranger)
  library(performanceEstimation)

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

  ## model tuning and training

  list(cache_path = c('./cache/ml_train.RData',
                      './cache/ml_ens.RData'),
       script_path = c('./ml scripts/development.R',
                       './ml scripts/ensemble.R'),
       message = c('Cached model tuning and training results',
                   'Cached ensemble tuning and training results')) %>%
    pwalk(access_cache)


  ## predictions, evaluation, and variable importance

  c('./ml scripts/predictions.R',
    './ml scripts/evaluation.R',
    './ml scripts/importance.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END --------

  insert_tail()
