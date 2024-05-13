# analysis globals:
#
# The cancer samples are stratified in the GENIE, TCGA, and Imvigor cohorts
# by presence of at least one mutation


  insert_head()

# container -------

  fgfr_globals <- list()

# stratification ------

  insert_msg('Stratification scheme')

  fgfr_globals$genetic_groups <- globals$cohort_expr %>%
    eval %>%
    map(~.x$mutation[c('sample_id', 'FGFR3')]) %>%
    map(set_names, c('sample_id', 'FGFR3_mutation')) %>%
    map(mutate,
        FGFR3_mutation = car::recode(FGFR3_mutation, "0 = 'WT'; 1 = 'mutated'"),
        FGFR3_mutation = factor(FGFR3_mutation, c('WT', 'mutated')))

# END ------

  insert_tail()
