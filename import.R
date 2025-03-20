# Import of the following bulk cancer data sets of mutations and expression:
#
# 1) GENIE, the subset with bladder and urothelial cancers
#
# 2) MSK 2022 bladder cancer cohort from cBioportal
#
# 3) TCGA BLCA, the Pan Cancer version deposited at cBioportal
#
# 4) Imvigor obtained via package IMvigor210CoreBiologies
#
# 5) BCAN metastatic bladder cancer cohort from cBioportal
#
# Only MIBC samples of the bladder are selected for further analyses.
#
# Assignment of the TCGA and IMvigor data to the molecular consensus clusters
# is done with consensusMIBC package based on ComBat-adjusted whole genome
# expression.
#
# We're importing IHC staining intensity scores for the Human Protein Atlas,
#
# For urothelial cancer cell lines of the Broad Institute collection (DepMap),
# whole-transcriptome gene expression data (RNA seq), mutation profiles,
# drug resistance data (GDSCT), as well as results of CRISPR and RNAi screening
# are imported.

# tools ------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(rstatix)
  library(figur)

  library(IMvigor210CoreBiologies)
  library(cbioportalR)
  library(xml2)
  library(readxl)
  library(gseaTools)

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

  ## only MIBC (selection not always possible) of the bladder
  # are selected for further analyses

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

# Prediction of drug sensitivity with the bulk cancer transcriptomes --------

  insert_msg('Prediction of drug response with the bulk cancer transcriptomes')

  list(cache_path = c('./data/gdsct.RData',
                      './data/ctrp2.RData',
                      './data/prism.RData'),
       script_path = c('./import scripts/predictions_gdsc_total.R',
                       './import scripts/predictions_ctrp2.R',
                       './import scripts/predictions_prism.R'),
       message = c('Loading drug sensitivity predictions, GDSCT',
                   'Loading drug sensitivity predictions, CTRP2',
                   'Loading drug sensitivity predictions, PRISM')) %>%
    pwalk(access_cache)

  ## building a single lexicon

  drugs <- list()

  drugs$lexicon <- list(GDSC1 = gdsct$lexicons$gdsc1,
                        GDSC2 = gdsct$lexicons$gdsc2,
                        CTRP2 = ctrp2$lexicon,
                        PRISM = prism$lexicon)

  drugs$lexicon[c("GDSC1", "GDSC2")] <-
    drugs$lexicon[c("GDSC1", "GDSC2")] %>%
    map(mutate,
        label = drug_name,
        moa = pathway_name)

  drugs$lexicon$PRISM <- drugs$lexicon$PRISM %>%
    mutate(drug_id = variable)

  drugs$lexicon <- drugs$lexicon %>%
    map(mutate, drug_id = as.character(drug_id)) %>%
    map(select,
        variable, drug_id, drug_name, label, moa, targets) %>%
    compress(names_to = 'experiment')

  ## appending the lexicon with the unit information
  ## plot titles, and labels for the heat maps with color-coding
  ## of the experiment

  drugs$lexicon <- drugs$lexicon %>%
    mutate(label = stri_capitalize_first(label),
           unit = ifelse(experiment != 'PRISM',
                         'log IC50', 'AUC'),
           experiment_color = car::recode(experiment,
                                          "'GDSC1' = '#581845';
                                          'GDSC2' = 'firebrick';
                                          'CTRP2' = 'darkolivegreen';
                                          'PRISM' = 'steelblue'"),
           plot_title = ifelse(drug_id != variable,
                               paste(label, drug_id, experiment, sep = ', '),
                               paste(label, experiment, sep = ', ')),
           plot_title = stri_capitalize_first(plot_title),
           hm_label = ifelse(drug_id != variable,
                             paste(label, drug_id, sep = ', '),
                             label),
           hm_label = stri_capitalize_first(hm_label),
           hm_label = paste0('<span style = "color:', experiment_color, '">',
                             hm_label, '</span>'))

  ## unified classification of mechanisms of action for the anti-cancer drugs.
  ## dictionaries: lists of drugs sharing a mechanism of action or a molecular
  ## target

  drugs$lexicon  <- drugs$lexicon %>%
    classify_moa

  drugs$moa_dictionary <- drugs$lexicon %>%
    blast(moa_short) %>%
    map(~.x$variable)

  drugs$target_dictinoary <- drugs$lexicon %>%
    make_target_dict

  ## a list of data frames with predictions
  ## by the successfully validated models

  drugs$predictions <-
    list(gdsc1 = gdsct$predictions$gdsc1,
         gdsc2 = gdsct$predictions$gdsc2,
         ctrp2 = ctrp2$predictions,
         prism = prism$predictions) %>%
    transpose %>%
    map(reduce, left_join, by = 'sample_id') %>%
    map(select, sample_id, all_of(drugs$lexicon$variable))

  ## identification of drugs that target FGFR according to
  ## the experiment's annotation

  drugs$fgfr_drugs <- drugs$lexicon %>%
    reglook('FGFR(1|2|3|4)') %>%
    .$variable

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

# Calculation of ssGSEA scores of FGFR-related gene signatures --------

  insert_msg('Caclulation of ssGSEA scvores of FGFR-related gene signatures')

  access_cache(cache_path = './data/sig.RData',
               script_path = './import scripts/signatures.R',
               message = 'Cached FGFR-related gene signature scores')

# END ------

  #rm(gdsct, ctrp2, prism)

  insert_tail()
