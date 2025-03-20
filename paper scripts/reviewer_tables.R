# Tables for the Reviewers

  insert_head()

# container -------

  rev_tabs <- list()

# characteristic of the cell lines ------

  insert_msg('Characteristic of the cell lines')

  rev_tabs$cell_lines <- expl_cecli$result_tbl %>%
    map_dfc(stri_replace_all, fixed = 'NaN%', replacement = 'not available') %>%
    as_mdtable(label = 'cell_lines',
               ref_name = 'cell_lines',
               caption = paste('Demographic, clinical, and pathological',
                               'characteristic of donors of the DepMap',
                               'urothelial cancer cell lines.',
                               'Numeric variables are presented as medians',
                               'with interquartile ranges and ranges.',
                               'Qualitative variables are shown as percentages',
                               'of the categories.'))

# Saving the tables on the disc --------

  insert_msg('Saving the tables on the disc')

  rev_tabs %>%
    save_excel(path = './report/reviewer_tables.xlsx',
               prefix = 'Reviewer Table ')

# END ------

  insert_tail()
