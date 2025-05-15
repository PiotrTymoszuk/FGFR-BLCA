# Descriptive demographic and clinical characteristic of the DepMap urothelial
# cancer cell lines

  insert_head()

# container ------

  expl_cecli <- list()

# analysis globals --------

  insert_msg('analysis globals')

  ## variable lexicon and demo/clinical data

  expl_cecli$lexicon <- depmap$clinic_lexicon %>%
    mutate(table_lab = ifelse(!is.na(unit),
                              paste(label, unit, sep = ', '),
                              label))

  expl_cecli$data <- depmap$cell_lines %>%
    map_dfc(function(x) if(is.character(x)) factor(x) else x)

# Descriptive stats -------

  insert_msg('descriptive stats')

  expl_cecli$stats <- expl_cecli$data %>%
    explore(variables = expl_cecli$lexicon$variable,
            what = 'table',
            pub_styled = TRUE) %>%
    format_res_tbl(dict = expl_cecli$lexicon,
                   value = 'table_lab')

# Result table ------

  insert_msg('Result table')

  ## appending the stats table with the total N number

  expl_cecli$result_tbl <-
    rbind(tibble(variable = 'Cell lines, N',
                 statistic = nrow(expl_cecli$data)),
          expl_cecli$stats) %>%
    format_res_tbl(dict = NULL) %>%
    set_names(c('Variable', 'Statistic'))

# END -------

  expl_cecli$data <- NULL
  expl_cecli$lexicon <- NULL

  expl_cecli <- compact(expl_cecli)

  insert_tail()
