# Analysis of co-occurrence and clustering of bladder cancers
# in respect to the most common genetic alterations
# (mutations, deletions, amplification) present in at least 5% of samples.
# The analysis is done solely for the GENIE and TCGA cohorts.
#
# Clustering is by latent class modeling with poLCA/polcaExtra with the minimal
# BIC as the model selection criterion.

## Work on it, propably the five cluster solution is better!!!

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
  library(fastTest)
  library(perich)

  library(survival)
  library(survminer)
  library(coxExtensions)

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
  ## and network analysis
  ##
  ## presence/absence of genetic alterations in the consensus classes

  list(cache_path = c('./cache/lca_occur.RData',
                      './cache/lca_sub.RData'),
       script_path = c('./latent class scripts/co_occurrence.R',
                       './latent class scripts/subtypes.R'),
       message = c('Loading cached co-occurrence results',
                   'Cached mutation enrichment for consensus classes')) %>%
    pwalk(access_cache)

  c('./latent class scripts/co_occurence_networks.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Definition and evaluation of the genetic classes -----

  insert_msg('Definition and evaluation of the latant classes')

  ## development of the clusters by latent class modeling
  ## predictions
  ## evaluation by comparing the class-defining variables
  ## between the classes

  list(cache_path = c('./cache/lca_dev.RData'),
       script_path = c('./latent class scripts/development.R'),
       message = c('Loading cached latent class analysis')) %>%
    pwalk(access_cache)

  c('./latent class scripts/prediction.R',
    './latent class scripts/evaluation.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Clinical, pathological, prognostic, and transcriptomic characteristic -------

  insert_msg('Clinical and pathological characteristic of the clusters')

  ## clinical characteristic, TMB, and survival analysis

  c('./latent class scripts/clinic.R',
    './latent class scripts/tmb.R',
    './latent class scripts/overall_survival.R',
    './latent class scripts/disease_relapse_free_survival.R',
    './latent class scripts/overall_survival_cox.R',
    './latent class scripts/overall_survival_lca_subsets.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## analysis of consensus molecular classes and differential gene expression,
  ## possible only for the TCGA collective

  c('./latent class scripts/genetic_molecular_subtypes.R',
    './latent class scripts/dge.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Genetic background of the clusters -------

  insert_msg('Genetic background of the clusters')

  ## analyses with weighted permutation testing

  list(cache_path = c('./cache/lca_mut.RData',
                      './cache/lca_del.RData',
                      './cache/lca_amp.RData'),
       script_path = c('./latent class scripts/mutations.R',
                       './latent class scripts/deletions.R',
                       './latent class scripts/amplifications.R'),
       message = paste('Loading cached',
                       c('mutation enrichment analysis',
                         'deletion enrichment analysis',
                         'amplification enrichment analysis'))) %>%
    pwalk(access_cache)

  ## visualizations

  c('./latent class scripts/genetic_plots.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
