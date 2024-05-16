# Renders the manuscript and supplements

  insert_head()

# Bibliography -----

  insert_msg('Reading the bibliography')

  blca_bib <- read_bib('./report/markdown/blca.bib') %>%
    as_mdbib

# Rendering the paper ------

  insert_msg('Rendering the analysis report')

  render('./report/markdown/report.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
           output_dir = './report')

# END ----

  insert_tail()
