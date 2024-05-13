# Import of the following data sets of mutations and expression:
#
# 1) GENIE, the subset with bladder and urothelial cancers
#
# 2) TCGA BLCA, the Pan Cancer version deposited at cBioportal
#
# 3) Imvigor obtained via package IMvigor210CoreBiologies
#
# Assignment of the TCGA and IMvigor data to the molecular consensus clusters
# is done with consensusMIBC package based on ComBat-adjusted whole genome
# expression

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(IMvigor210CoreBiologies)
  library(cbioportalR)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  library(htGLMNET)
  library(consensusMIBC)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# import scripts -------

  insert_msg('Import scripts')

  list(cache_path = c('./data/genie.RData',
                      './data/tcga.RData',
                      './data/imvigor.RData'),
       script_path = c('./import scripts/genie.R',
                       './import scripts/tcga.R',
                       './import scripts/imvigor.R'),
       message = paste('Loading cached',
                       c('GENIE BLCA study',
                         'TCGA BLCA study',
                         'IMvigor study'))) %>%
    pwalk(access_cache)

# ComBat expression adjustment and assignment to the molecular subtypes -------

  insert_msg('ComBat and molecular subtypes')

  access_cache(cache_path = './data/combat.RData',
               script_path = './import scripts/combat.R',
               message = 'Loading cached ComBat adjutment results')

  c('./import scripts/subtypes.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
