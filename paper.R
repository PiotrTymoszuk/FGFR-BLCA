# Elements of the short correspondence paper

  insert_head()

# tools --------

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

# figure, supplementary figures and supplementary tables -------

  insert_msg('Figure, Supplementary Figures, Supplementary Tables')

  c('./paper scripts/figure.R') %>%
    source_all(message = TRUE, crash = TRUE)

# rendering the legends and supplementary material ------

  insert_msg('Rendering manuscript parts and Supplementary Material')

  c('./paper scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
