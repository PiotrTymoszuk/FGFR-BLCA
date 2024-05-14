# Figures, tables and the analysis report

# tools -------

  library(tidyverse)
  library(trafo)
  library(rlang)
  library(stringi)

  library(cowplot)
  library(figur)
  library(flextable)
  library(writexl)

  library(rmarkdown)
  library(knitr)
  library(bookdown)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Utils, figures and tables -----

  insert_msg('Utils, figures and tables')

  c('./report scripts/links.R',
    './report scripts/html.R') %>%
    source_all(message = TRUE, crash = TRUE) %>%
    print

  c('./report scripts/figures.R',
    './report scripts/tables.R') %>%
    source_all(message = TRUE, crash = TRUE) %>%
    print

# Analysis report -------

  insert_msg('Analysis report')

  c('./report scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
