# characteristic of the urothelial cancer cell line collective in the GDSC
# drug screening experiment.

  insert_head()

# container --------

  expl_cells <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  expl_cells$lexicon <- depmap$lexicon

  expl_cells$data <-
    depmap[c("gdsc", "prism", "crispr", "expression", "mutation")] %>%
    reduce(full_join, by = 'sample_id')

  expl_cells$lexicon <- expl_cells$lexicon

  expl_cells$variables <- expl_cells$lexicon %>%
    blast(type) %>%
    map(~.x$variable)

# descriptive stats for the drug resistance -------

  insert_msg('Descriptive stats')

  expl_cells$stats <-
    explore(expl_cells$data,
            split_factor = NULL,
            variables = reduce(expl_cells$variables, union),
            what = 'table',
            pub_styled = TRUE) %>%
    format_res_tbl(dict = expl_cells$lexicon,
                   value = 'table_lab')

# Variability of expression and mutation variables -------

  insert_msg('Variability of expression variables')

  expl_cells$data[expl_cells$variables$mutation] <-
    expl_cells$data[expl_cells$variables$mutation] %>%
    map_dfc(~as.numeric(.x) - 1)

  expl_cells$distribution <- expl_cells$data %>%
    select(all_of(expl_cells$variables$expression),
           all_of(expl_cells$variables$mutation)) %>%
    distr_stats

  ## selection of modeling variables: non-zero variance features

  expl_cells$mod_variables <- expl_cells$distribution %>%
    filter(freqRatio > 0) %>%
    .$variable

# Normality of distribution of the modeling responses --------

  insert_msg('Normality of distribution of the responses')

  ## there are, in part, serious violations of normality,
  ## despite log2 transformation

  ## Shapiro-Wilk test

  expl_cells$normality <-
    expl_cells$data[c(expl_cells$variables$resistance,
                      expl_cells$variables$crispr)] %>%
    f_shapiro_test(as_data_frame = TRUE) %>%
    as_tibble

  ## QQ plots

  expl_cells$qq_plots <-
    expl_cells$data[c(expl_cells$variables$resistance,
                      expl_cells$variables$crispr)] %>%
    map(eda) %>%
    map(plot, type = 'qq')

# END -------

  expl_cells$data <- NULL
  expl_cells$lexicon <- NULL

  expl_cells <- compact(expl_cells)

  insert_tail()
