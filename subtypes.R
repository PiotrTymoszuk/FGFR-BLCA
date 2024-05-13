# Differences in frequency of genetic alterations and expression of the FGF
# and FGFR genes between the consensus molecular classes.

# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)

  library(exda)
  library(microViz)

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

  c('./subtype scripts/genetic.R',
    './subtype scripts/dge.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()


