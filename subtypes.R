# Differences in frequency of genetic alterations and expression of the FGF
# and FGFR genes between the consensus molecular classes.

# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(exda)
  library(microViz)
  library(fastTest)
  library(perich)

  library(furrr)

  library(ggtext)
  library(ggrepel)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts -----

  insert_msg('Analysis scripts')

  ## genetics of the consensus subtypes

  access_cache(cache_path = './cache/sub_genet.RData',
               script_path = './subtype scripts/genetic.R',
               message = 'Cached results for genetic alterations')

  ## clinical characteristic of the consensus classes
  ## differential expression of FGF/FGFR/FGFBP-related genes
  ## co-expression networks for the genes of interest
  ## differences in predicted drug resistance between the consensus classes,
  ## enrichment analysis of mechanisms of action and molecular targets
  ## for the common significant drugs.
  ## differences in ssGSEA scores of FGFR-reated gene signatures between
  ## consensus classes of MIBC

  c('./subtype scripts/clinic.R',
    './subtype scripts/dge.R',
    './subtype scripts/coexpression_networks.R',
    './subtype scripts/drugs.R',
    './subtype scripts/moa_enrichment.R',
    './subtype scripts/signatures.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()


