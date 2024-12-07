# Tables for the analysis report

  insert_head()

# container ------

  rep_tabs <- list()

# Genes of interest --------

  insert_msg('Genes of interest')

  ## FGF15 is absent from the human genoe

  rep_tabs$genes <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    map(~tibble(`Gene symbol` = .x)) %>%
    map(filter, `Gene symbol` != 'FGF15') %>%
    compress(names_to = 'Protein type') %>%
    mutate(`Protein type` = car::recode(`Protein type`,
                                        "'binding_proteins' = 'binding proteins'"),
           `Entrez ID` = mapIds(org.Hs.eg.db,
                                keys = `Gene symbol`,
                                keytype = 'SYMBOL',
                                column = 'ENTREZID')) %>%
    relocate(`Protein type`) %>%
    as_mdtable(label = 'genes',
               ref_name = 'genes',
               caption = paste('FGFR-, FGF-, and FGFBP-coding genes',
                               'investigated in the current report.'))

# characteristic of the cohorts -------

  insert_msg('Cohort characteristic')

  rep_tabs$cohorts <- expl_cohorts$result_tbl %>%
    set_names(c('Variable',
                globals$cohort_labs[c("genie", "msk",
                                      "tcga", "imvigor",
                                      "bcan", "hpa")],
                'Significance',
                'Effect size')) %>%
    mdtable(label = 'cohort_characteristic',
            ref_name = 'cohorts',
            caption = paste('Characteristic of the investigated cohorts.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))

# Frequency of genetic alterations -------

  insert_msg('Frequency of genetic features')

  rep_tabs$genetics <- expl_genet$result_tbl %>%
    compress(names_to = 'Alteration type') %>%
    relocate(`Alteration type`) %>%
    mdtable(label = 'genetic_features',
            ref_name = 'genetics',
            caption = paste('Frequency of the most common somatic mutations',
                            'and copy number variants in the',
                            'investigated cohorts.',
                            'Frequencies of alterations present in at least 5%',
                            'of cancer samples are displayed.',
                            'Statistical significance for differences between',
                            'the cohorts was assessed by chi-squared test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Frequency of genetic alterations of the FGF- and FGFR-coding genes -------

  insert_msg('Frequency of alterations in FGF and FGFR genes')

  rep_tabs$fgf_genetics <- expl_genet$fgf_result_tbl %>%
    compress(names_to = 'Alteration type') %>%
    relocate(`Alteration type`) %>%
    mdtable(label = 'genetic_fgf_fgfr_features',
            ref_name = 'fgf_genetics',
            caption = paste('Frequency of somatic mutations and copy number',
                            'variants of FGF-, FGFR-, and FGFBP-coding genes',
                            'in the investigated cohorts.',
                            'Statistical significance for differences between',
                            'the cohorts was assessed by chi-squared test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Classification of somatic mutations of the FGFR genes ------

  insert_msg('Classification of mutations')

  rep_tabs$fgfr_class <- expl_fgfr$result_tbl

  for(i in names(rep_tabs$fgfr_class)) {

    names(rep_tabs$fgfr_class[[i]])[3] <-
      'Mutation classification'

  }

  rep_tabs$fgfr_class <- rep_tabs$fgfr_class %>%
    set_names(c('Mutation type', 'Variant type', 'Protein domain')) %>%
    compress(names_to = 'Classification scheme') %>%
    relocate(`Classification scheme`)

  rep_tabs$fgfr_class <- rep_tabs$fgfr_class %>%
    mdtable(label = 'fgfr_mutation_classification',
            ref_name = 'fgfr_class',
            caption = paste('Frequency of mutations in the FGFR1/2/3/4 genes',
                            'in the GENIE BLCA, MSK, TCGA BLCA, IMvigor,',
                            'and BCAN cohorts',
                            'classified by mutation type, variant type, and',
                            'protein domain.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Expression of genes and proteins -------

  insert_msg('Expression, genes and proteins')

  rep_tabs$expression <- expl_expr$stats %>%
    set_names(c('Variable', globals$cohort_labs[names(expl_expr$stats)[-1]])) %>%
    as_mdtable(label = 'expression',
               ref_name = 'expression',
               caption = paste('Expression of FGFR-, FGF-, and FGFBP-coding',
                               'genes at mRNA level was quantified in the',
                               'TCGA BLCA, IMvigor, and BCAN cohort by',
                               'RNA sequencing. The mRNA amounts were',
                               'provided as log2-transformed gene counts',
                               '(TCGA BLCA and IMvigor) or Z-scores of log2',
                               'gene counts (BCAN).',
                               'Protein expression of the genes of interest was',
                               'measured in the Human Protein Atlas (HPA) by',
                               'immunohistochemistry (IHC).',
                               'The protein amounts are metered by the IHC score',
                               'defined as a product of staining',
                               'intensity and quantity of positive cells.',
                               'Median mRNA and protein amounts with',
                               'interquartile ranges, ranges, and numbers of',
                               'cancer samples are listed.',
                               'The table is available as a supplementary Excel',
                               'file.'))

# Differential gene expression, FGFR3 WT and FGFR3 mutated -------

  insert_msg('Differential gene expresssion, FGFR3 WT and mutated')

  rep_tabs$fgfr3_dge <- expl_dge$result_tbl$FGFR3_mut %>%
    mdtable(label = 'fgfr3_differential_gene_expression',
            ref_name = 'fgfr3_dge',
            caption = paste('Differential expression of FGF-, FGFR-, and',
                            'FGFBP-coding genes in urothelial cancers with',
                            'and without FGFR3 mutations.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance was determined by',
                            'two-tailed T test with',
                            "Cohen's d effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Clinical characteristic of samples with and without FGFR3 mutations ------

  insert_msg('Clinical characteristic, FGFR3 WT and FGFR3 mutated samples')

  rep_tabs$fgfr3_clinic <- fgfr_clinic$result_tbl

  names(rep_tabs$fgfr3_clinic)[c(3, 4)] <-
    c('FGFR3 wild-type', 'FGFR3 mutated')

  rep_tabs$fgfr3_clinic <- rep_tabs$fgfr3_clinic %>%
    mdtable(label = 'fgfr3_clinical_features',
            ref_name = 'fgfr3_clinic',
            caption = paste('Demographic, clinical, and pathological',
                            'characteristic of cancers with WT and mutated',
                            'FGFR3.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))


# Genetic alterations in the consensus molecular subtypes ----------

  insert_msg('Gnetics of the consensus classes')

  rep_tabs$cons_genetics <- sub_genet$result_tbl %>%
    mdtable(label = 'consensus_subsets_genetics',
            ref_name = 'cons_genetics',
            caption = paste('Frequency of somatic mutations and copy number',
                            'variants of FGFR-, FGF-, and FGFBP-coding genes',
                            'in the consensus molecular classes of',
                            'urothelial cancers.',
                            'The frequencies are presented as percentages of the',
                            'consensus class.',
                            'Statistical significance of differences between',
                            'the consensus classes was determined by',
                            'weighted permutation test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Differential gene expression in the consensus subsets -------

  insert_msg('Differential gene expression in the consensus subtypes')

  rep_tabs$cons_dge <- sub_dge$result_tbl %>%
    mdtable(label = 'consensus_subsets_expression',
            ref_name = 'cons_dge',
            caption = paste('Expression of FGF-, FGFR-, and FGFBP-coding genes',
                            'in the consensus molecular classes of urothelial',
                            'cancers.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance of differences between',
                            'the molecular classes was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Modeling of the consensus subsets: model evaluation -----

  insert_msg('Evaluation of models of the consensus subsets')

  rep_tabs$ml_eval <- ml_eval$stats %>%
    filter(algorithm != 'ensemble') %>%
    select(algorithm, cohort, Accuracy, Kappa, brier_score,
           AUC, Mean_Sensitivity, Mean_Specificity) %>%
    mutate(cohort = globals$cohort_labs[cohort],
           algorithm = car::recode(algorithm,
                                   "'elnet' = 'Elastic Net';
                                   'ranger' = 'Random Forest';
                                   'ensemble' = 'Ensemble'")) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
    set_names(c('Algorithm', 'Cohort', 'Accuracy', 'Kappa', 'Brier score',
                'AUC', 'Mean sensitivity', 'Mean specificity')) %>%
    as_mdtable(label = 'ml_eval',
               ref_name = 'ml_eval',
               caption = paste('Performance of Elastic Net and Random Forest',
                               'models of LumP, LUMU, stroma-rich,',
                               'and BA/Sq consensus',
                               'molecular classes of urothelial cancers that',
                               'uses expression of FGFR-, FGF-, and FGFBP-coding',
                               'genes as the sole explanatory factors.',
                               'The models were trained in the TCGA BLCA cohort,',
                               'the gene expression values were ComBat-processed',
                               'to minimize inter-cohort variability.',
                               'Minority consensus classes were amplified by the',
                               'SMOTE algorithm to improve their identification',
                               'by the model.',
                               'The models were evaluated in the TCBA BLCA,',
                               'IMvigor, and BCAN cohorts; the consensus class',
                               'assignment predicted by a nearest centroid',
                               'algorithm from consensusMIBC package served as',
                               'the ground truth.',
                               'Metrics of model accuracy, calibration, and',
                               'receiver-operating characteristic are listed.'))

# Non-zero beta coefficients of the Elastic Net model of consensus classes -------

  insert_msg('Non-zero betas')

  rep_tabs$ml_imp_elnet <- ml_imp$coefs %>%
    map(select, variable, beta) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('Variable', names(ml_imp$coefs))) %>%
    map_dfc(function(x) if(is.numeric(x)) as.character(signif(x, 3)) else x) %>%
    as_mdtable(label = 'ml_imp_elnet',
               ref_name = 'ml_imp_elnet',
               caption = paste('Coefficients of the Elastic Net model of',
                               'consensus molecular classes of urothelial',
                               'cancers.',
                               'The Elastic Net model was developed in the',
                               'TCGA BLCA cohort with expression values of genes',
                               'coding for FGFR, FGF, and FGFBP proteins as',
                               'the sole explanatory factors.',
                               'Values of beta of the coefficients are listed.'))

# Permutation variable importance of the Random Forest model --------

  insert_msg('Permutation importance for the Random Forest model')

  rep_tabs$ml_imp_ranger <- ml_imp$ranger_importance %>%
    arrange(-importance) %>%
    mutate(importance = as.character(signif(importance, 2))) %>%
    set_names(c('Varaible', 'Permutation importance')) %>%
    as_mdtable(label = 'ml_imp_ranger',
               ref_name = 'ml_imp_ranger',
               caption = paste('Permutation variable importance of the',
                               'Random Forest model of consensus molecular',
                               'classes of urothelial cancers.',
                               'The Random Forest model was developed in the',
                               'TCGA BLCA cohort with expression values of genes',
                               'coding for FGFR, FGF, and FGFBP proteins as',
                               'the sole explanatory factors.',
                               'Permutation importance was computed as',
                               'differences of errors of the genuine model',
                               'and models, in which the explanatory variables',
                               'were re-shuffled at random.'))

# Evaluation of the Elastic Net Cox model of overall survival --------

  insert_msg('Evaluation of the Elastic Net model of overall survival')

  rep_tabs$surv_eval <- surv_eval$stats %>%
    mutate(cohort = globals$cohort_labs[cohort],
           n_events = as.character(n_events),
           n_complete = as.character(n_complete)) %>%
    select(cohort, n_complete, n_events, c_index, raw_rsq, ibs_model) %>%
    map_dfc(function(x) if(is.numeric(x)) as.character(signif(x, 2)) else x) %>%
    set_names(c('Cohort',
                'Patients, N',
                'Events, N',
                'Concordance index',
                'R-square',
                'Integrated Brier score')) %>%
    as_mdtable(label = 'surv_eval',
               ref_name = 'surv_eval',
               caption = paste('Performance of an Elastic Net Cox model',
                               'using expression of FGFR-, FGF-, and FGFBP-coding',
                               'genes as explanatory factors at prediction of',
                               'overall survival.',
                               'The Elastic Net model was trained in the TCGA BLCA',
                               'cohort with ComBat-processed log2-transformed',
                               'expression values as firs- and second-order',
                               'terms.',
                               'Metrics of concordance, explanatory performance,',
                               'and overall calibration are summarized.'))

# Non-zero coefficients of the survival model -------

  insert_msg('Non-zero coefficients of the survival model')

  rep_tabs$surv_imp <- surv_imp$coefs %>%
    arrange(-beta) %>%
    transmute(`Gene symbol` = stri_replace(variable,
                                           regex = '_sec$',
                                           replacement = ''),
              `Term order` = ifelse(stri_detect(variable, regex = '_sec$'),
                                    'second', 'first'),
              `Coefficient beta` = as.character(signif(beta, 3)),
              `Hazard ratio` = as.character(signif(hr, 3))) %>%
    as_mdtable(label = 'surv_imp',
               ref_name = 'surv_imp',
               caption = paste('Non-zero coefficients of the Elastic Net Cox',
                               'model of overall survival with FGFR-, FGF-,',
                               'and FGFBP-coding genes as explanatory factors.'))

# KM analysis of survival as a function of expression --------

  insert_msg('KM analysis for FGFR, FGF, and FGFBP gene expression')

  rep_tabs$surv_km <- surv_km$result_tbl %>%
    filter(`Gene symbol` %in% c('TNFAIP6',
                                'GPC1',
                                'FIBP',
                                'FGF11',
                                'FGFBP1')) %>%
    as_mdtable(label = 'surv_km',
               ref_name = 'surv_km',
               caption = paste('Univariable, model-free analysis of overall',
                               'survival in urothelial cancer patients',
                               'stratified by expression of FGFR-,',
                               'FGF-, and FGFR-coding genes.',
                               'Patients were stratified as low and high',
                               'expressors of a gene by a cutoff value that',
                               'corresponds to the maximum log rank statistics.',
                               'Statistical significance of differences in',
                               'overall survival between low and high',
                               'expressors was determined by Peto-Peto test',
                               'corrected for multiple comparisons with the',
                               'false discovery rate (FDR) method.',
                               'The stratification cutoffs, numbers of high',
                               'and low expressors, numbers of events in the',
                               'strata, median survival times with 95%',
                               'confidence intervals, as well as raw and',
                               'FDR-corrected p values are summarized.',
                               'The results are shown for the genes with',
                               'significant effects in at least two out of',
                               'the TCGA, IMvigor, and BCAN cohorts.'))

# Genetics of the genetic subsets --------

  insert_msg('Genetics of the genetic subsets')

  rep_tabs$lca_genetics <- lca_eval$result_tbl %>%
    mdtable(label = 'genetic_subsets_genetics',
            ref_name = 'lca_genetics',
            caption = paste('Frequency of the class-defining somatic mutations',
                            'and copy number alterations in the genetic subsets',
                            'of urothelial cancers.',
                            'The frequencies are presented as percentages of',
                            'samples with genetic alterations in the genetic',
                            'subsets.',
                            'Differences between the genetic subsets were',
                            'assessed by chi-square test with',
                            "Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Differential gene expression in the genetic subsets -------

  insert_msg('Genetic subsets: differential gene expression')

  rep_tabs$lca_dge <- lca_dge$result_tbl %>%
    mdtable(label = 'genetic_subsets_gene_expression',
            ref_name = 'lca_dge',
            caption = paste('Expression of FGF-, FGFR-, and FGFBP-coding genes',
                            'in the genetic subsets of urothelial cancers of',
                            'the TCGA BLCA cohort.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'Statistical significance of differences between',
                            'the genetic subsets was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Clinical characteristic of the genetic subsets ------

  insert_msg('Clinical characteristic of the genetic subsets')

  rep_tabs$lca_clinic <- lca_clinic$result_tbl %>%
    mdtable(label = 'genetic_subsets_clinics',
            ref_name = 'lca_clinic',
            caption = paste('Demographic, clinical, and pathological',
                            'characteristic of the genetic subsets of',
                            'urothelial cancers.',
                            'Numeric features are presented as medians with',
                            'interquartile ranges and ranges.',
                            'Qualitative variables are shown as percentages',
                            'and observation counts of the categories.'))

# Saving the tables on the disc ------

  insert_msg('Saving the tables on the disc')

  rep_tabs %>%
    save_excel(path = './report/tables.xlsx')

# END -----

  insert_tail()
