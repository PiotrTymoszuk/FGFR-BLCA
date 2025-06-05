# Overall survival analysis for FGF-, FGFR-, and FGFBP-coding genes expression
# in the TCGA, IMvigor, and BCAN cohort. Two approaches:
#
# 1) Multi-parameter modeling via Elastic Net Cox regression with
# ComBat-processed log2 expression of the genes of interest as explanatory
# factors. The model is tuned by repeated CV and trained in the TCGA BLCA cohort.
# Subsequently, its linear predictor scores in the TCGA, IMvigor, and BCAN
# cohorts serve as explanatory variables for 'canonical' Cox models,
# whose performance is evaluated by C-index, IBC, and R-square, as well as
# graphical performance and calibration checks.
#
# 2) Univariable Kaplan-Meier analysis of overall survival in the TCGA, IMvigor,
# and BCAN cohort. Statistical significance is determined by Peto-Peto test.
# Effect shared by at least two cohorts are of interest.

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(survival)
  library(survminer)
  library(coxExtensions)

  library(glmnet)
  library(caret)

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

# analysis globals -------

  insert_msg('Analysis globals')

  c('./survival scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts ---------

  insert_msg('Analysis scripts')

  ## multi-parameter modeling

  access_cache(cache_path = './cache/surv_train.RData',
               script_path = './survival scripts/development.R',
               message = 'Cached results of survival model development')

  c('./survival scripts/evaluation.R',
    './survival scripts/importance.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## univariable KM analysis

  c('./survival scripts/univariable.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END --------

  insert_tail()
