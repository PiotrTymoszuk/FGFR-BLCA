# Import of the following data sets of mutations and expression:
#
# 1) GENIE, the subset with bladder and urothelial cancers
#
# 2) MSK 2022 bladder cancer cohort from cBioportal
#
# 3) TCGA BLCA, the Pan Cancer version deposited at cBioportal
#
# 4) Imvigor obtained via package IMvigor210CoreBiologies
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
  library(xml2)

  library(AnnotationDbi)
  library(org.Hs.eg.db)

  library(htGLMNET)
  library(consensusMIBC)
  library(impute)

  library(furrr)

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
                      './data/msk.RData',
                      './data/tcga.RData',
                      './data/imvigor.RData',
                      './data/bcan.RData'),
       script_path = c('./import scripts/genie.R',
                       './import scripts/msk.R',
                       './import scripts/tcga.R',
                       './import scripts/imvigor.R',
                       './import scripts/bcan.R'),
       message = paste('Loading cached',
                       c('GENIE BLCA study',
                         'MSK 2022 study',
                         'TCGA BLCA study',
                         'IMvigor study',
                         'BCAN study'))) %>%
    pwalk(access_cache)

# ComBat expression adjustment and assignment to the molecular subtypes -------

  insert_msg('ComBat and molecular subtypes')

  access_cache(cache_path = './data/combat.RData',
               script_path = './import scripts/combat.R',
               message = 'Loading cached ComBat adjustment results')

  c('./import scripts/subtypes.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Import of IHC scores for Human Protein Atlas -------

  insert_msg('Import of HPA data')

  access_cache(cache_path = './data/hpa.RData',
               script_path = './import scripts/hpa.R',
               message = paste('Cached Human Protein Atlas data'))

# Import of GDSC drug screening information -------

  insert_msg('Import of GDSC data')

  access_cache(cache_path = './data/depmap.RData',
               script_path = './import scripts/cells.R',
               message = 'Cached GDSC data')

  depmap$raw_expression <- NULL
  depmap <- compact(depmap)

# END ------

  insert_tail()
