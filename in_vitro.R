# Analyses and visualizations of results of in vitro treatments of urothelial
# cancer cell lines with WT and mutated FGFR3 with Erdafitinib

# tools -------

  library(tidyverse)
  library(trafo)
  library(stringi)
  library(rlang)

  library(ggrepel)
  library(ggtext)
  library(figur)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  map <- purrr::map

  insert_head()

  c('./tools/globals.R',
    './tools/functions.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis and visualization scripts ------

  insert_msg("Analysis and visualization scripts")

  c("./in vitro scripts/plots.R") %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
