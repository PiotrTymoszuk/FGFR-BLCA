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
# 6) Expression of mRNA and proteins of FGFR-, FGF-, and FGFBP-coding genes.
#
# 7) Co-expression of FGF and FGFR gene mRNA in the cancer tissue: correlation
# and network analysis.
#
# 8) Differential expression of FGFR and FGF genes in cancer samples with and
# without mutations in the FGFR1 - 4.
#
# 9) Co-occurrence of alterations in FGF and FGFR genes: multi-dimensional
# scaling and a network analysis.
#
# 10) Characteristic and distribution tests for urothelial cancer cell lines
# in the GDSC experiment, choice of the modeling features,
# correlation of the expression variables.
#
# 11) Distribution tests for ssGSEA scores of FGFR-related gene signatues.

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

# Analysis globals -------

  insert_msg('Analysis globals')

  c('./exploration scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Scripts for the exploratory data analysis --------

  insert_msg('Analysis scripts')

  ## comparison of the cohorts

  c('./exploration scripts/cohorts.R',
    './exploration scripts/subtypes.R',
    './exploration scripts/survival.R') %>%
    source_all(message = TRUE, crash = TRUE) %>%
    print

  ## mutations and copy number alterations, identification
  ## of mutational hot spots in the FGFR genes

  c('./exploration scripts/genetics.R',
    './exploration scripts/details.R',
    './exploration scripts/mutation_hotspots.R') %>%
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

  ## network analyses for the bulk cancers:
  ## FGFR/FGF-related transcripts and genetic alterations

  c('./exploration scripts/coexpression_networks.R',
    './exploration scripts/alteration_networks.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## predictions of drug resistance: distribution tests
  ## PCA, correlation of predictions for drugs that target FGFRs,
  ## key predictor genes for models of anti-FGFR drugs

  c('./exploration scripts/drugs.R') %>%
    source_all(message = TRUE, crash = TRUE)

  rm(gdsct, ctrp2, prism)

  ## cell lines: characteristic, distribution tests,
  ## choice of the modeling features, co-expression networks,
  ## and analyses of distribution/background of CRISPR and RNAi gene effects
  ## and drug resistance

  c('./exploration scripts/cells_clinic.R',
    './exploration scripts/cells.R',
    './exploration scripts/network_cells.R',
    './exploration scripts/crispr.R',
    './exploration scripts/resistance.R',
    './exploration scripts/rnai.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## distribution tests for the ssGSEA scores of FGFR-related gene signatures

  c('./exploration scripts/signatures.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
