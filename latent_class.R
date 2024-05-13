# Analysis of co-occurrence and clustering of bladder cancers
# in respect to the most common genetic alterations
# (mutations, deletions, amplification) present in at least 5% of samples.
# The analysis is done solely for the GENIE and TCGA cohorts.
#
# Clustering is by latent class modeling with poLCA/polcaExtra with the minimal
# BIC as the model selection criterion.

# tools -------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(polcaExtra)
  library(clustTools)

  library(ranger)
  library(caret)

  library(exda)

  library(survival)
  library(survminer)

  library(microViz)

  library(ggrepel)
  library(ggtext)

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

# analysis globals: ready to use data frames with genetic features -------

  insert_msg('Analysis globals')

  c('./latent class scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scrips: co-occurrence and consensus molecular subtypes -------

  insert_msg('Co-occurrence and consensus molecular subtypes')

  ## general co-occurrence of the most common genetic features

  access_cache(cache_path = './cache/lca_occur.RData',
               script_path = './latent class scripts/co_occurrence.R',
               message = 'Loading cached co-occurrence results')

  ## presence/absence of genetic alterations in the consensus classes

  c('./latent class scripts/subtypes.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Definition and evaluation of the genetic classes -----

  insert_msg('Definition and evaluation of the latant classes')

  ## development of the clusters by hierarchical clustering
  ## of spectralization results in the training GENIE cohort

  access_cache(cache_path = './cache/lca_dev.RData',
               script_path = './latent class scripts/development.R',
               message = 'Loading cached spectral clustering results')

  ## prediction of the clusters in the TCAG cohort,
  ## characteristic of the clusters

  c('./latent class scripts/prediction.R',
    './latent class scripts/evaluation.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Clinical, pathological and prognostic characteristic -------

  insert_msg('Clinical and pathological characteristic of the clusters')

  ## survival analysis: the 11q13amp and 11q13amp + 9p21del groups are merged

  c('./latent class scripts/clinic.R',
    './latent class scripts/tmb.R',
    './latent class scripts/overall_survival.R',
    './latent class scripts/disease_relapse_free_survival.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## analysis of consensus molecular classes and differential gene expression,
  ## possible only for the TCGA collective

  c('./latent class scripts/genetic_molecular_subtypes.R',
    './latent class scripts/dge.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
