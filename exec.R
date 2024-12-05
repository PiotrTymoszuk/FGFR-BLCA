# Launches the entire analysis pipeline

  library(soucer)

  print(source_all(c('import.R',
                     'exploration.R',
                     'FGFR3_analysis.R',
                     'subtypes.R',
                     'ml.R',
                     'survival.R',
                     'latent_class.R',
                     'report.R'),
                   message = TRUE, crash = TRUE))
