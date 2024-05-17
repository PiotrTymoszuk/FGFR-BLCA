# Exploratory data analysis:
#
# 1) Clinical and pathological characteristic of the cohorts.
#
# 2) Distribution of the consensus molecular subtypes.
#
# 3) Differences in survival between the cohorts.
#
# 4) General frequency of mutations and copy number alterations. Details for
# mutations and CNA of the FGFR genes.
#
# 5) Characteristic of the mutation types and location for the FGFR genes.
#
# 6) Co-expression of FGF and FGFR gene mRNA in the cancer tissue.
#
# 7) Differential expression of FGFR and FGF genes in cancer samples with and
# without mutations in the FGFR1 - 4.
#
# 8) Co-occurrence of alterations in FGF and FGFR genes.

# tools -------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(exda)
  library(survival)
  library(survminer)

  library(jsonlite)

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

  ## co-expression and differential gene expression
  ## co-occurrence of alterations in the FGF and FGFR genes

  c('./exploration scripts/coexpression.R',
    './exploration scripts/dge.R') %>%
    source_all(message = TRUE)

  access_cache(cache_path = './cache/expl_occur.RData',
               script_path = './exploration scripts/co_occurrence.R',
               message = 'Loading cached co-occurrence analysis results')

# END -----

  insert_tail()
