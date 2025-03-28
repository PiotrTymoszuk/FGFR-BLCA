# Renders the manuscript and supplements

  insert_head()

# Bibliography -----

  insert_msg('Reading the bibliography')

  blca_bib <- read_bib('./report/markdown/blca.bib') %>%
    as_mdbib

# Rendering: figure 1 legend, rebuttal letter and supplements ------

  insert_msg('Rendering manuscript parts')

  render('./report/markdown/figure1_legend.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
           output_dir = './report')

  render('./report/markdown/rebuttal_letter.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './report')

  render('./report/markdown/supplementary_material.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './report')

# END ----

  insert_tail()
