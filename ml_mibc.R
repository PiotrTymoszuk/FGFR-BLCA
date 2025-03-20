# Machine learning modeling of MIBC consensus classes with genes proposed by
# Kamoun et al. as explanatory factors
#
# To this end, we're developing a multinomial Elastic Net and a Random Forest
# model of the consensus classes in the training portion TCGA BLCA cohort
# (2/3 of observations).
# ComBat log2-transformed expression of the genes of interest serve as
# explanatory factors.
# Because of an imbalance of consensus class sizes, we're amplifying the
# minority classes (LumU and Stroma-rich) by SMOTE.
# The models are tuned (optimal lambda and RF hyper-parameters) by
# minimizing deviance in repeated cross-validation and by minimizing
# classification errors of OOB predictions for, respectively, Elastic Net and
# Random Forest models.
# Performance of the model is evaluated in the TCGA (training and test portion),
# IMvigor, and BCAN cohorts, where the consensus cluster assignment predicted
# by consensusMIBC package is assumed as the ground truth.
#
# We're essentially replicating the machine learning pipeline developed for the
# modeling with FGFR/FGF/FGFBP-related genes and compare these models with
# the models developed with the full set of genes proposed by Kamoun et al.

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

  c('./MIBC ml scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts -------

  insert_msg('Analysis scripts')

  ## model tuning and training, SHAP importance

  list(cache_path = c('./cache/mibc_train.RData',
                      './cache/mibc_ens.RData',
                      './cache/mibc_shap.RData'),
       script_path = c('./MIBC ml scripts/development.R',
                       './MIBC ml scripts/ensemble.R',
                       './MIBC ml scripts/shap.R'),
       message = c('Cached model tuning and training results',
                   'Cached ensemble tuning and training results',
                   'Cached SHAP variable importance')) %>%
    pwalk(access_cache)


  ## predictions, evaluation, and variable importance
  ##
  ## comparison of performance of the models employing the FGFR-related genes
  ## and the models employing the full set of genes proposed by Kamoun et al.

  c('./MIBC ml scripts/predictions.R',
    './MIBC ml scripts/evaluation.R',
    './MIBC ml scripts/importance.R',
    './MIBC ml scripts/shap_plots.R',
    './MIBC ml scripts/comparison.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
