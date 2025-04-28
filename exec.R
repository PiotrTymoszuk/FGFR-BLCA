# Launches the entire analysis pipeline

  library(soucer)

  print(source_all(c('import.R', ## import of the clinical, genomic, and transcriptomic data
                     'exploration.R', ## exploratory analysis
                     'FGFR3_analysis.R', ## characteristic of cancers with and without FGFR3 mutations
                     'subtypes.R', ## clinical, genetic, and transcriptomic characteristic of MIBC consensus classes
                     'ml.R', ## machine learning: prediction of MIBC consensus classes with FGFR/FGF/FGFBP genes
                     'ml_mibc.R', ## machine learning: prediction of MIBC consensus classes with genes by Kamoun et al. 2020
                     'paper.R'), ## figures, supplementary material
                   message = TRUE, crash = TRUE))
