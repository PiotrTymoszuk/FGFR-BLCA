# Supplementary Tables for the brief report manuscript

  insert_head()

# container ------

  suppl_tabs <- list()

# Characteristic of the bulk cancer and IHC cohorts -------

  insert_msg('Charcterisric of bulk cancer and IHC cohorts')

  suppl_tabs$cohorts <- expl_cohorts$result_tbl %>%
    as_mdtable(label = 'cohorts',
               ref_name = 'cohorts',
               caption = paste('Characteristic of the investigated cohorts.',
                               'Numeric features are presented as medians with',
                               'interquartile ranges and ranges.',
                               'Qualitative variables are shown as percentages',
                               'and observation counts of the categories.'))

# Clinical characteristic of the consensus molecular classes of MIBC ------

  insert_msg('Clinical charactaristic of the consensus molecular classes')

  suppl_tabs$sub_clinic <- sub_clinic$result_tbl %>%
    as_mdtable(label = 'sub_clinic',
               ref_name = 'sub_clinic',
               caption = paste('Characteristic of the consensus molecular',
                               'classes of MIBC in the TCGA BLCA, IMvigor,',
                               'and BCAN cohorts.',
                               'Numeric features are presented as medians with',
                               'interquartile ranges and ranges.',
                               'Qualitative variables are shown as percentages',
                               'and observation counts of the categories.'))

# Genes of interest --------

  insert_msg('Genes of interest')

  ## FGF15 is absent from the human genome

  suppl_tabs$genes <-
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

# Frequency of genetic alterations in the FGFR/FGF/FGFBP genes -------

  insert_msg('Frequency of genetic alterations in FGFR/FGF/FGFBP genes')

  suppl_tabs$genetic_alterations <- expl_genet$fgf_result_tbl %>%
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

# Alterations of the FGFR/FGF/FGFBP genes in the consensus classes of MIBC ------

  insert_msg('Genetics of FGFR/FGF/FGFBP in the consensus MIBC classes')

  suppl_tabs$sub_genetics <- sub_genet$result_tbl %>%
    mdtable(label = 'sub_genetics',
            ref_name = 'sub_genetics',
            caption = paste('Frequency of somatic mutations and copy number',
                            'variants of FGFR-, FGF-, and FGFBP-coding genes',
                            'in the consensus molecular classes of MIBC',
                            'in the TCGA BLCA, IMvigor, and BCAN cohorts.',
                            'The frequencies are presented as percentages of',
                            'the consensus class.',
                            'Statistical significance of differences between',
                            'the consensus classes was determined by',
                            'weighted permutation test',
                            "with Cramer's V effect size statistic.",
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'The table is available as a supplementary',
                            'Excel file.'))

# Expression of genes and proteins -------

  insert_msg('Expression, genes and proteins')

  suppl_tabs$expression <- expl_expr$stats %>%
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

# Differential expression of FGFR/FGF/FGFBP genes in the consensus classes ------

  insert_msg('Differential expression of FGFR/FGF/FGFBP genes in the classes')

  suppl_tabs$sub_dge <- sub_dge$result_tbl %>%
    mdtable(label = 'sub_dge',
            ref_name = 'sub_dge',
            caption = paste('Expression of FGF-, FGFR-, and FGFBP-coding genes',
                            'in the consensus molecular classes of MIBC',
                            'in the TCGA BLCA, IMvigor, and BCAN cohorts.',
                            'Statistical significance of differences between',
                            'the molecular classes was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'log2-transformed expression levels are presented',
                            'as medians with interquartile ranges and ranges.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Signatures of FGFR-related genes in the consensus classes ------

  insert_msg('Signatures of FGFR-related genes in the classes')

  suppl_tabs$sub_sig <- sub_sig$result_tbl %>%
    mdtable(label = 'sub_sig',
            ref_name = 'sub_sig',
            caption = paste('Single sample gene set enrichment analysis',
                            '(ssGSEA) scores of FGFR-related gene signatures',
                            'in the consensus molecular classes of MIBC',
                            'in the TCGA BLCA, IMvigor, and BCAN cohorts.',
                            'Statistical significance of differences between',
                            'the molecular classes was determined by one-way',
                            'ANOVA with eta-square effect size statistic.',
                            'P values were corrected for multiple testing with',
                            'the false discovery rate method.',
                            'ssGSEA scores are presented',
                            'as medians with interquartile ranges and ranges.',
                            'The table is available as a supplementary Excel',
                            'file.'))

# Forecasts of resistance in the consensus molecular classes, anti-FGFR drugs -------

  insert_msg('Forecasts of anti-FGFR drug resistance, consensus classes')

  suppl_tabs$sub_drugs <- sub_drugs$result_tbl %>%
    filter(Variable %in% exchange(drugs$fgfr_drugs,
                                  drugs$lexicon,
                                  value = 'plot_title')) %>%
    mdtable(label = 'sub_drugs',
            ref_name = 'sub_drugs',
            caption = paste('Predicted resistance metrics to pan-FGFR',
                            'inhibitors in the consensus molecular classes of',
                            'MIBC in the TCGA BLCA, IMvigor, and BCAN cohorts.',
                            'The predictions of were made with RIDGE regularized',
                            'linear models trained with transcriptome data of',
                            'epithelial cancer cell lines treated with FGFR',
                            'inhibitors in GDSC1, GDSC2, CTRP2, and PRISM',
                            'drug screening experiments.',
                            'Resistance metrics are presented',
                            'as medians with interquartile ranges and ranges.'))

# Tuning of the machine learning models --------

  insert_msg('tuning of the machine learning models')

  suppl_tabs$tuning <-
    list(mibc_genes = mibc_train,
         fgfr_genes = ml_train) %>%
    map(~.x$tune_obj) %>%
    map(map,
        ~.x$best_tune) %>%
    map(map, map_dfc, function(x) if(is.numeric(x)) signif(x, 3) else x) %>%
    map(map,
        ~.x[!names(.x) %in% c('cvm',
                              'mse_oob',
                              'class_error_oob',
                              'model_id')]) %>%
    map(~map2(map(.x, names), .x, paste, sep = ' = ')) %>%
    map(map_chr, paste, collapse = '\n') %>%
    map(compress,
        names_to = 'algorithm',
        values_to = 'hyper_params') %>%
    compress(names_to = 'model') %>%
    transmute(`Explanatory factors/genes` = mibc_comp$model_labels[model],
              Algorithm = globals$algo_labs[algorithm],
              `Selection criterion` = ifelse(algorithm == 'elnet',
                                             'model deviance\nrepeated cross-validation',
                                             'classification error\nout-of-bag predictions'),
              `Hyper-parameters` = hyper_params)

  suppl_tabs$tuning <- suppl_tabs$tuning %>%
    as_mdtable(label = 'tuning',
               ref_name = 'tuning',
               caption = paste('Selection of the optimal sets of',
                               'hyper-parameters for Elastic Net and Random',
                               'Forest models of the MIBC consensus molecular',
                               'classes by cross-validation and out-of-bag',
                               'prediction tuning in the training subset of',
                               'the TCGA BLCA cohort.'))

# Performance of machine learning models -------

  insert_msg('Performance of machine learning models')

  suppl_tabs$ml <- mibc_comp$stats %>%
    mutate(cohort = factor(cohort,
                           c('tcga_train', 'tcga_test', 'imvigor', 'bcan')),
           algorithm = factor(algorithm, c('elnet', 'ranger')),
           model = factor(model, c('mibc_genes', 'fgfr_genes'))) %>%
    filter(!is.na(cohort), !is.na(algorithm), !is.na(model)) %>%
    arrange(model, algorithm, cohort) %>%
    transmute(`Explanatory factors/genes` = mibc_comp$model_labels[as.character(model)],
              Algorithm = globals$algo_labs[as.character(algorithm)],
              `Cohort, subset` = globals$cohort_labs[as.character(cohort)],
              Accuracy = Accuracy,
              Kappa = Kappa,
              `Brier score` = brier_score,
              AUC = AUC,
              `Mean sensitivity` = Mean_Sensitivity,
              `Mean specificity` = Mean_Specificity) %>%
    map_dfc(function(x) if(is.numeric(x)) as.character(signif(x, 2)) else x)

  suppl_tabs$ml <- suppl_tabs$ml %>%
    as_mdtable(label = 'ml',
               ref_name = 'ml',
               caption = paste('Performance of Elastic Net and Random Forest',
                               'models of LumP, LUMU, stroma-rich,',
                               'and BA/Sq consensus',
                               'molecular classes of MIBC.',
                               'The models were trained in a random subset of',
                               'two-thirds observations of the TCGA BLCA cohort',
                               'with ComBat-processed log2 expression values of',
                               'the consensus class-refining genes by Kamoun',
                               'et al. and the FGFR/FGF/FGFBP genes as',
                               'explanatory variables.',
                               'Minority consensus classes were amplified by the',
                               'SMOTE algorithm to improve their identification',
                               'by the model.',
                               'The models were evaluated in the training',
                               'and test subsets of the TCGA BLCA,',
                               'IMvigor, and BCAN cohorts; the consensus class',
                               'assignment predicted by a nearest centroid',
                               'classifier from consensusMIBC package served as',
                               'the ground truth.',
                               'Metrics of model accuracy, calibration, and',
                               'receiver-operating characteristic are listed.'))

# SHAP variable importance -------

  insert_msg('SHAP variable importance')

  suppl_tabs$shap <-
    list(mibc_genes = mibc_splots$stats %>%
           map(transmute,
               variable = variable,
               class = class,
               importance = abs(importance)),
         fgfr_genes = ml_splots$stats) %>%
    map(compress, names_to = 'algorithm') %>%
    compress(names_to = 'model') %>%
    left_join(ml_globals$lexicon) %>%
    transmute(`Explanatory factors/genes` = mibc_comp$model_labels[model],
              Algorithm = globals$algo_labs[algorithm],
              `Consensus class` = class,
              `Gene symbol` = variable,
              `FGFR/FGF/FGFBP protein type` = ifelse(is.na(protein_type),
                                                     'non-FGFR/FGF/FGFBP',
                                                     protein_type),
              `Importance, mean absolute SHAP` = signif(importance, 2)) %>%
    as_mdtable(label = 'shap',
               ref_name = 'shap',
               caption = paste('Importance of explanatory variables for the',
                               'Elastic Net and Random Forest models of',
                               'MIBC consensus molecular classes estimated by',
                               'the SHAP algorithm.',
                               'Contribution of the genuine MIBC consensus',
                               'class-defining genes and the FGFR-, FGF-, and',
                               'FGFBP-coding genes to prediction of consensus',
                               'molecular classes was assessed by mean absolute',
                               'values of SHAP metrics (Shapley additive',
                               'explanations).',
                               'The table is available as a supplementary',
                               'Excel file.'))

# Saving the tables on the disc --------

  insert_msg('Saving the tables on the disc')

  suppl_tabs %>%
    save_excel(path = './report/supplementary_tables.xlsx',
               prefix = 'Supplementary Table S')

# END -------

  insert_tail()
