# Analysis scripts for patients stratified by
# presence or absence of the FGFR3 mutations (at least one per gene).
# For stratification details, please refer
# to the 'globals.R' script in the 'FGFR analysis scripts' folder.
#
# 1) Demographic, clinical, and pathological characteristics of the patients
#
# 2) Total mutation burden, an therapy response in the mutation strata
#
# 3) Relapse-free, disease-specific and overall survival in the mutation strata
#
# 4) Effects of amplification of the FGF3/4/19 genes on expression of genes
# coding for FGFR, FGF, and FGFBP proteins


# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(exda)
  library(microViz)
  library(fastTest)

  library(survival)
  library(survminer)

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

# analysis globals ------

  insert_msg('Analysis globals')

  c('./FGFR3 analysis scripts/globals.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts -------

  insert_msg('Analysis scripts')

  ## clinical characteristic, therapy response, and survival

  c('./FGFR3 analysis scripts/clinic.R',
    './FGFR3 analysis scripts/icb_response.R',
    './FGFR3 analysis scripts/overall_survival.R',
    './FGFR3 analysis scripts/disease_relapse_free_survival.R') %>%
    source_all(message = TRUE, crash = TRUE)

  ## TMB, effects of FGF3/4/19 amplification on gene expression

  c('./FGFR3 analysis scripts/tmb.R',
    './FGFR3 analysis scripts/amplification_gene_expression.R') %>%
    source_all(message = TRUE)

# END -----

  insert_tail()
