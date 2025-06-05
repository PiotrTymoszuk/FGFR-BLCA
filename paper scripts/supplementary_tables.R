# Supplementary Tables for the brief report manuscript

  insert_head()

# containers ------

  suppl_tabs <- list()
  suppl_globals <- list()

# some globals ----

  insert_msg('Table globals')

  ## levels of classes and cohorts
  ## for arr<anging the supplementary tables

  suppl_globals$class_levels <- c('LumP', 'LumU', 'Stroma-rich', 'Ba/Sq')

  suppl_globals$cohort_levels <- globals$analysis_cohorts

  suppl_globals$cohort_lab_levels <-
    unname(globals$cohort_labs[suppl_globals$cohort_levels])

# Characteristic of the bulk cancer and IHC cohorts -------

  insert_msg('Charcterisric of bulk cancer and IHC cohorts')

  ## we don't show hypothesis testing, according to the EAU guidelines
  ## it is quite obvious that the cohorts are pretty different
  ## due to the study designs

  suppl_tabs$cohorts <- expl_cohorts$result_tbl %>%
    select(-`FDR p value`, -`Effect size`) %>%
    as_mdtable(label = 'cohorts',
               ref_name = 'cohorts',
               caption = paste("Characteristic of the investigated cohorts.",
                               "Numeric variables are presented as medians with",
                               "interquartile ranges.",
                               "Qualitative variables are shown as counts of",
                               "the categories and percentages of the categories",
                               "within the complete observation set."))

# Clinical characteristic of the consensus molecular classes of MIBC ------

  insert_msg('Clinical charactaristic of the consensus molecular classes')

  suppl_tabs$sub_clinic <- sub_clinic$result_tbl %>%
    as_mdtable(label = 'sub_clinic',
               ref_name = 'sub_clinic',
               caption = paste("Characteristic of the major consensus molecular",
                               "classes of MIBC (luminal papillary [LumP],",
                               "luminal unstable [LumU], stroma-rich,",
                               "basal/squamous-like [Ba/Sq]) in the",
                               "TCGA BLCA, IMvigor, and BCAN cohorts.",
                               "Numeric variables are presented as medians with",
                               "interquartile ranges.",
                               "Qualitative variables are shown as counts of",
                               "the categories and percentages of the categories",
                               "within the complete observation set.",
                               "Statistical significance of differences between",
                               "the consensus molecular classes were assessed by",
                               "Kruskal-Wallis and \u03C7\u00B2 test with,",
                               "respectively, \u03B7\u00B2 and Cramer's V",
                               "effect size statistics.",
                               "P values were corrected for multiple testing",
                               "with the false discovery rate method (FDR)."))

# Genes of interest --------

  insert_msg('Genes of interest')

  ## FGF15 is absent from the human genome
  ## compacted representation: gene symbol (Entrez ID)

  suppl_tabs$genes <-
    globals[c("receptors", "ligands", "binding_proteins")] %>%
    map(~mapIds(org.Hs.eg.db,
                keys = .x,
                keytype = 'SYMBOL',
                column = 'ENTREZID')) %>%
    map(function(x) map2_chr(names(x), x,
                             ~paste0("_", .x, "_ (", .y, ")"))) %>%
    map_chr(paste, collapse = ", ") %>%
    compress(names_to = "Protein type",
             values_to = "Gene symbol (Entrez ID)") %>%
    mutate(`Protein type` = stri_replace_all(`Protein type`,
                                             fixed = "_",
                                             replacement = " "))

  suppl_tabs$genes <- suppl_tabs$genes %>%
    as_mdtable(label = 'genes',
               ref_name = 'genes',
               caption = paste('FGFR-, FGF-, and FGFBP-coding genes',
                               'investigated in the current report.'))

# Frequency of genetic alterations in the FGFR/FGF/FGFBP genes -------

  insert_msg('Frequency of genetic alterations in FGFR/FGF/FGFBP genes')

  ## to meet the request by a reviewer: showing only the alterations
  ## presenting at least 1% of the samples

  suppl_tabs$genetic_alterations <-
    map2(expl_genet$fgf_result_tbl,
         expl_genet$fgf_top_alterations,
         ~filter(.x, `Gene symbol` %in% .y)) %>%
    compress(names_to = 'alteration') %>%
    mutate(`Gene symbol, alteration` =
             paste0("_", `Gene symbol`, "_ ", alteration),
           `Gene group` = factor(`Gene group`,
                                 c("receptor", "ligand", "binding protein")),
           alteration = factor(alteration,
                               c("mutation", "deletion", "amplification"))) %>%
    arrange(alteration, `Gene group`, `Gene symbol`) %>%
    select(-`Gene symbol`, -alteration) %>%
    relocate(`Gene group`, `Gene symbol, alteration`)

  suppl_tabs$genetic_alterations <- suppl_tabs$genetic_alterations %>%
    mdtable(label = 'genetic_fgf_fgfr_features',
            ref_name = 'fgf_genetics',
            caption = paste("Frequency of somatic mutations and copy number",
                            "variants of FGFR-, FGF-, and FGFBP-coding genes.",
                            "Counts and percentages of samples with alterations",
                            "for the mutations, amplifications, and deletions",
                            "found in at least 1%",
                            "of specimens in the BCAN, GENIE BLCA, IMvigor,",
                            "MSK IMPACT, and TCGA BLCA are presented."))

# Alterations of the FGFR/FGF/FGFBP genes in the consensus classes of MIBC ------

  insert_msg('Genetics of FGFR/FGF/FGFBP in the consensus MIBC classes')

  ## showing the same genetic alterations as in the previous table
  ## (mutations of FGFR3, amplification of the receptors,
  ## amplification of FGF3/4/19)

  suppl_tabs$sub_genetics <- sub_genet$result_tbl %>%
    mutate(Alteration = factor(Alteration,
                               c('mutation', 'deletion', 'amplification'))) %>%
    blast(Alteration)

  suppl_tabs$sub_genetics <-
    map2_dfr(suppl_tabs$sub_genetics,
             expl_genet$fgf_top_alterations[names(suppl_tabs$sub_genetics)],
             ~filter(.x, `Gene symbol` %in% .y)) %>%
    mutate(`Gene symbol, alteration` = paste0("_", `Gene symbol`,
                                              "_ ", Alteration)) %>%
    relocate(`Gene symbol, alteration`, Cohort)

  suppl_tabs$sub_genetics <- suppl_tabs$sub_genetics %>%
    mutate(Alteration = factor(Alteration, unique(Alteration)),
           Cohort = factor(Cohort, suppl_globals$cohort_lab_levels),
           consensusClass = stri_split_fixed(`MIBC consensus class (n samples)`,
                                             pattern = " ",
                                             simplify = TRUE)[, 1],
           consensusClass = factor(consensusClass,
                                   suppl_globals$class_levels)) %>%
    arrange(Alteration, `Gene symbol`, Cohort, consensusClass) %>%
    select(-consensusClass, -`Gene symbol`, -Alteration)

  suppl_tabs$sub_genetics <- suppl_tabs$sub_genetics %>%
    mdtable(label = 'sub_genetics',
            ref_name = 'sub_genetics',
            caption = paste("Frequency of somatic mutations and copy number",
                            "variants of FGFR1/3 and FGF3/4/19 genes",
                            "in the luminal papillary (LumP), luminal unstable (LumU),",
                            "stroma-rich, and basal/squamous like (Ba/Sq)",
                            "consensus molecular classes of MIBC",
                            "in the TCGA BLCA, IMvigor, and BCAN cohorts.",
                            "The frequencies are presented as counts and",
                            "percentages of samples with alterations within",
                            "the consensus class.",
                            "Statistical significance of differences between",
                            "the consensus classes was determined by",
                            "weighted permutation test",
                            "with Cramer's V effect size statistic.",
                            "P values were corrected for multiple testing with",
                            "the false discovery rate (FDR) method."))

# Expression of genes and proteins -------

  insert_msg('Expression, genes and proteins')

  ## showing the entire table as requested by the reviewer

  suppl_tabs$expression <- expl_expr$result_tbl %>%
    mutate(!!names(.)[ncol(.)] :=
             stri_replace(.data[[names(.)[ncol(.)]]],
                          fixed = "\n",
                          replacement = ", ")) %>%
    as_mdtable(label = 'expression',
               ref_name = 'expression',
               caption = paste("Expression of FGFR-, FGF-, and FGFBP-coding",
                               "genes and proteins in the MIBC tissue and",
                               "urothelial cancer cell lines.",
                               "Expression at mRNA level was quantified in the",
                               "TCGA BLCA, IMvigor and BCAN cohorts,",
                               "and in DepMap cell lines by RNA sequencing.",
                               "The mRNA amounts were",
                               "provided as log2-transformed gene counts",
                               "(TCGA BLCA, IMvigor, DepMap) or Z-scores of",
                               "log2 gene counts (BCAN).",
                               "Protein expression was measured in the Human",
                               "Protein Atlas (HPA) by immunohistochemistry",
                               "(IHC).",
                               "The protein amounts was metered by the IHC score",
                               "defined as a product of staining",
                               "intensity and quantity of positive cells.",
                               "Median mRNA and protein amounts with",
                               "interquartile ranges, and numbers of",
                               "cancer samples are listed.",
                               "Numbers of all available samples in the",
                               "cohorts are provided in the table caption.",
                               "For the IHC data in the HPA cohort,",
                               "numbers of specimens with evaluable staining",
                               "results are indicated for each gene."))

# Differential expression of FGFR/FGF/FGFBP genes in the consensus classes ------

  insert_msg('Differential expression of FGFR/FGF/FGFBP genes in the classes')

  ## showing the to genes found to be differentially regulated between the
  ## consensus classes in at least two cohorts

  suppl_tabs$sub_dge <- sub_dge$result_tbl %>%
    filter(Variable %in% c('Samples, N', sub_dge$top_common_significant)) %>%
    mdtable(label = 'sub_dge',
            ref_name = 'sub_dge',
            caption = paste("Expression of FGFR-, FGF-, and FGFBP-coding genes",
                            "in the consensus molecular classes of MIBC",
                            "(luminal papillary [LumP], luminal unstable [LumU],",
                            "stroma-rich, basal/squamous-like [Ba/Sq])",
                            "in the TCGA BLCA, IMvigor, and BCAN cohorts.",
                            "Statistical significance of differences between",
                            "the consensus molecular classes was determined by",
                            "one-way ANOVA with \u03B7\u00B2 effect size statistic.",
                            "P values were corrected for multiple testing with",
                            "the false discovery rate (FDR) method.",
                            "log2-transformed expression levels are presented",
                            "as medians with interquartile ranges for",
                            "the top genes whose expression differentiated",
                            "between the molecular consensus classes",
                            "(top 10 largest \u03B7\u00B2 values shared by",
                            "all cohorts)."))

# Signatures of FGFR-related genes in the consensus classes ------

  insert_msg('Signatures of FGFR-related genes in the classes')

  ## showing the top differentially regulated signatures
  ## to meet the reviewer's request

  suppl_tabs$sub_sig <- sub_sig$result_tbl %>%
    filter(Variable %in%
             c('Samples, N',
               sub_sig$top_common_significant %>%
                 reduce(union) %>%
                 exchange(sig$lexicon))) %>%
    select(-`Member genes`) %>%
    map_dfc(stri_replace, fixed = "(", replacement = "\n(") %>%
    mdtable(label = 'sub_sig',
            ref_name = 'sub_sig',
            caption = paste("Single sample gene set enrichment analysis",
                            "(ssGSEA) scores of FGFR-related gene signatures",
                            "in the consensus molecular classes of MIBC",
                            "(luminal papillary [LumP], luminal unstable [LumU],",
                            "stroma-rich, basal/squamous-like [Ba/Sq])",
                            "in the TCGA BLCA, IMvigor, and BCAN cohorts.",
                            "Statistical significance of differences between",
                            "the consensus molecular classes was determined by",
                            "one-way ANOVA with eta-square effect size statistic.",
                            "P values were corrected for multiple testing with",
                            "the false discovery rate (FDR) method.",
                            "ssGSEA scores are presented as medians with",
                            "interquartile ranges for",
                            "the top signatures whose scores differentiating",
                            "between the molecular consensus classes",
                            "(top 10 largest \u03B7\u00B2 values shared by",
                            "all cohorts)."))

# Forecasts of resistance in the consensus molecular classes, anti-FGFR drugs -------

  insert_msg('Forecasts of anti-FGFR drug resistance, consensus classes')

  suppl_tabs$sub_drugs <- sub_drugs$result_tbl %>%
    filter(Variable %in% c("Samples, N",
                           exchange(drugs$fgfr_drugs,
                                  drugs$lexicon,
                                  value = 'plot_title'))) %>%
    mdtable(label = 'sub_drugs',
            ref_name = 'sub_drugs',
            caption = paste("Predicted metrics of resistance to pan-FGFR",
                            "inhibitors in the consensus molecular classes of MIBC",
                            "(luminal papillary [LumP], luminal unstable [LumU],",
                            "stroma-rich, basal/squamous-like [Ba/Sq])",
                            "in the TCGA BLCA, IMvigor, and BCAN cohorts.",
                            "The predictions of IC50 (concentration resulting",
                            "in 50% growth inhibition) and AUC (area under the",
                            "dose-response curve) were made by RIDGE regularized",
                            "linear models trained with transcriptome data of",
                            "epithelial cancer cell lines treated with FGFR",
                            "inhibitors in GDSC1, GDSC2, CTRP2, and PRISM",
                            "drug screening experiments.",
                            "Statistical significance of differences between",
                            "the consensus molecular classes was determined by",
                            "one-way ANOVA with \u03B7\u00B2 effect size statistic.",
                            "P values were corrected for multiple testing with",
                            "the false discovery rate (FDR) method.",
                            "Resistance metrics are presented",
                            "as medians with interquartile ranges."))

# Tuning of the machine learning models --------

  insert_msg('tuning of the machine learning models')

  suppl_tabs$tuning <-
    list(mibc_genes = mibc_train,
         fgfr_genes = ml_train) %>%
    map(~.x$tune_obj) %>%
    map(map,
        ~.x$best_tune) %>%
    map(map, map_dfc, function(x) if(is.numeric(x)) signif(x, 2) else x) %>%
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
                                             'minimal model deviance\nrepeated cross-validation',
                                             'minimal classification error\nout-of-bag predictions'),
              `Hyper-parameters` = hyper_params)

  suppl_tabs$tuning <- suppl_tabs$tuning %>%
    as_mdtable(label = 'tuning',
               ref_name = 'tuning',
               caption = paste("Selection of the optimal sets of",
                               "hyper-parameters for Elastic Net and Random",
                               "Forest models of the MIBC consensus molecular",
                               "classes by cross-validation and out-of-bag",
                               "tuning in the training subset of",
                               "the TCGA BLCA cohort."))

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
               caption = paste("Performance of Elastic Net and Random Forest",
                               "models of the luminal papillary (LumP),",
                               "luminal unstable (LumU), stroma-rich,",
                               "and basal/squamous-like (Ba/Sq) consensus",
                               "molecular classes of MIBC.",
                               "The models were trained in a random subset of",
                               "two-thirds of observations of the TCGA BLCA cohort",
                               "with ComBat-processed log2 expression values of",
                               "the consensus class-defining genes by Kamoun",
                               "et al. or the FGFR/FGF/FGFBP genes as",
                               "explanatory variables.",
                               "Minority consensus molecular classes",
                               "(LumU and stroma-rich) were oversampled by the",
                               "SMOTE algorithm to improve their identification",
                               "by the models.",
                               "The models were evaluated in the training",
                               "and test subsets of the TCGA BLCA,",
                               "IMvigor, and BCAN cohorts; the consensus class",
                               "assignment predicted by a nearest centroid",
                               "classifier from consensusMIBC package served as",
                               "the ground truth.",
                               "Metrics of model accuracy, calibration, and",
                               "receiver-operating characteristic are listed."))

# SHAP variable importance -------

  insert_msg('SHAP variable importance')

  ## to meet the Reviewer's requirements, creating two separate tables:
  ## 1) For the FGFR/FGF/FGFBP classifiers
  ## 2) For the classifiers employing the consensus class-defining genes by
  ## Kamoun et al.
  ## in each table we present top five most influential genes

  suppl_tabs[c('shap', 'kamoun_shap')] <-
    list(ml_splots$stats, mibc_splots$stats) %>%
    map(map, mutate, importance = abs(importance)) %>%
    map(map, group_by, class) %>%
    map(map, slice_max, importance, n = 5) %>%
    map(map, ungroup) %>%
    map(compress, names_to = 'algorithm') %>%
    map(left_join, ml_globals$lexicon, by = 'variable') %>%
    map(transmute,
        Algorithm = globals$algo_labs[algorithm],
        `Consensus MIBC class` = factor(class, suppl_globals$class_levels),
        `Gene symbol` = variable,
        `Protein type` = ifelse(is.na(protein_type),
                                'non FGFR/FGF/FGFBP',
                                stri_replace_all(protein_type,
                                                 fixed = '_',
                                                 replacement = ' ')),
        `Importance, mean |SHAP|` = signif(importance, 2))

  suppl_tabs[c('shap', 'kamoun_shap')] <-
    suppl_tabs[c('shap', 'kamoun_shap')] %>%
    map(arrange, Algorithm, `Consensus MIBC class`) %>%
    list(label = names(.),
         ref_name = names(.),
         caption = paste('Importance of explanatory variables for the',
                         'Elastic Net and Random Forest models of',
                         'MIBC consensus molecular classes estimated by',
                         'the SHAP algorithm.',
                         'The models employed',
                         c('the FGFR-, FGF-, and FGFBP-coding genes',
                           'the genuine MIBC consensus class-defining genes by Kamoun et al.'),
                         'as explanatory factors.',
                         'Contribution of the explanatory factors to predictions',
                         'of consensus molecular classes was assessed by mean absolute',
                         'values of SHAP metrics (Shapley additive',
                         'explanations).',
                         'Top five most influential genes for the predictions of',
                         'luminal papillary (LumP), luminal unstable (LumU),',
                         'stroma-rich, and basal/squamous-like (Ba/Sq) classes',
                         'by the Elastic Net and Random Forest algorithms are',
                         'listed.')) %>%
    pmap(as_mdtable)

# Saving the tables on the disc --------

  insert_msg('Saving the tables on the disc')

  suppl_tabs %>%
    save_excel(path = './report/supplementary_tables.xlsx',
               prefix = 'Supplementary Table S')

# END -------

  insert_tail()
