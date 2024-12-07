# Exploratory data analysis:
#
# 1) Clinical and pathological characteristic of the cohorts.
#
# 2) Distribution of the consensus molecular subtypes.
#
# 3) Differences in survival between the cohorts.
#
# 4) General frequency of mutations and copy number alterations. Details for
# mutations and CNA of the FGFR genes. Expression levels of the genes and
# proteins of interest
#
# 5) Characteristic of the mutation types and location for the FGFR genes.
#
# 6) Co-expression of FGF and FGFR gene mRNA in the cancer tissue: correlation
# and network analysis.
#
# 7) Differential expression of FGFR and FGF genes in cancer samples with and
# without mutations in the FGFR1 - 4.
#
# 8) Co-occurrence of alterations in FGF and FGFR genes: multi-dimensional
# scaling and a network analysis.

# tools -------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(exda)
  library(survival)
  library(survminer)
  library(fastTest)

  library(igraph)
  library(graphExtra)

  library(jsonlite)
  library(cowplot)

  library(microViz)

  library(ggrepel)
  library(ggtext)
  library(figur)

  library(furrr)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Scripts for the exploratory data analysis --------

  insert_msg('Analysis scripts')

  ## comparison of the cohorts

  c('./exploration scripts/cohorts.R',
    './exploration scripts/subtypes.R',
    './exploration scripts/survival.R') %>%
    source_all(message = TRUE, crash = TRUE) %>%
    print

  ## mutations and copy number alterations

  c('./exploration scripts/genetics.R',
    './exploration scripts/details.R') %>%
    source_all(message = TRUE, crash = TRUE) %>%
    print

  ## overall expression levels of the features of interest: RNA and IHC
  ## co-expression and differential gene expression
  ## co-occurrence of alterations in the FGF and FGFR genes

  c('./exploration scripts/expression.R',
    './exploration scripts/ihc_representative.R',
    './exploration scripts/coexpression.R',
    './exploration scripts/dge.R') %>%
    source_all(message = TRUE)

  access_cache(cache_path = './cache/expl_occur.RData',
               script_path = './exploration scripts/co_occurrence.R',
               message = 'Loading cached co-occurrence analysis results')

  ## network analyses

  c('./exploration scripts/co_expression_networks.R',
    './exploration scripts/alteration_networks.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
