# FGFR-BLCA
Genetic alterations and expression of genes coding for FGF reseptors, ligands and binding proteins in urothelial cancer

## Results

FGFR signaling is of immense importance for urotheliaal cancer development and progression. This is underlined by a recent approval of erdafitinib, a pan-FGFR inhibitor, for treatment of urothelial cancers that carry _FGFR3_ mutations. 

The primary goal of the data mining project was to investigate the genetic and transcriptional background of FGFR signaling in urothelial cancers. The secondary goal was to harness the genetic variability of such malignancies by an analysis of co-occuring and mutually exclusive somatic mutations, amplifications, deletons, and chromosomal aberrations. 

Those goals were addressed by an analysis of genetic aletartions and differential expression of 45 genes coding for FGFR, FGF, and FGF-binding proteins (FGFBP) in five publicly available genetic and transcriptomic collectives: 

* a subset of urothelial cancer collected in the GENIE project [^1]. For this collective, detailed genetic alteration information was available
* MSK IMPACT genetic cohort of urothelical cancers [^2]
* TCGA BLCA cohort which provided both genetic and transcriptomic analysis [^3][^4]
* IMvigor collective of late-stage urothelial cancers treated with anti-PD-L1 immunotherapy [^5]
* BCAN cohort of metastatic neoplasms with predominantly systemic chemotherapy [^6]

Results of differential analyses of genetic alterations and gene expression let us propose FGFR1 and FGFR3 two hubs of FGFR signaling in urothelical cancers. Those proteins get triggered by at least distinct machenisms found to be specific for the consensus molecular subtypes by Kamoun et al [^7]: 

1. In luminal papillary cancers, FGFR3 is the main player of FGFR signaling. Because of the frequent mutations of the transmembrane/hinge region and association with proteoglycans and co-receptors, it may be activated in a ligand-independent way
2. In stroma-rich cancers, FGFR signaling is driven by overexpression of FGFR1, its ligands (FGF2, FGF7, FGF9, FGF10), and co-activating molecules (ANOS1, TGFBR3, HSPG)
3. In basal/squamous cancers, FGFR1 signaling is mediated by overexpression of FGF5 and the binding proteins FGFBP1

![FGFR_signaling_blca](https://github.com/user-attachments/assets/b87b6e82-b0af-4ac7-b496-3fc325ab2176)

By means of machine learning models trained with levels of 38 FGFR-, FGF-, and FGFBP-coding transcript, we could reproduce nearly completely the consensus molecular classification scheme by Kamoun et al [^7]. This suggests that FGFR signaling with its various activation and modulaton mechanisms is the central driver of phynotypic variability of urothelial cancers.

## Analysis report

Figures with the analysis results are available as PDF files [__here__](https://github.com/PiotrTymoszuk/FGFR-BLCA/tree/main/report/figures).
Figure PBG files are available from [__a Dropbox folder__](https://www.dropbox.com/scl/fo/78ufkox671pg26p4ik9gt/AFpKfdibezDv_xxW3WayXyc?rlkey=ou8335apj52d4eadnonyhiksj&dl=0).
Analysis report is available [__as a downloadable docx file__](https://github.com/PiotrTymoszuk/FGFR-BLCA/blob/main/report/report.docx). 
Tables as provided [__as an Excel file__](https://github.com/PiotrTymoszuk/FGFR-BLCA/blob/main/report/tables.xlsx).

## Manuscript parts

Parts/panels of Figure 1 of the manuscript are available as [__PDF files__](https://github.com/PiotrTymoszuk/FGFR-BLCA/tree/main/report/paper%20figure%201) and [__PNG files__](https://github.com/PiotrTymoszuk/FGFR-BLCA/tree/main/report/paper%20PNG%20figure%201).

## Usage

The analysis pipeline requires some development packages, the easiest way to install them is to use `devtools`:

```r

devtools::install_github('PiotrTymoszuk/trafo')
devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/soucer')
devtools::install_github('PiotrTymoszuk/figur')
devtools::install_github('PiotrTymoszuk/microViz')
devtools::install_github('PiotrTymoszuk/coxExtensions')
devtools::install_github('PiotrTymoszuk/htGLMNET')
devtools::install_github('PiotrTymoszuk/ExDA')
devtools::install_github('PiotrTymoszuk/fastTest')
devtools::install_github('PiotrTymoszuk/perich')
devtools::install_github('PiotrTymoszuk/graphExtra')

```
To launch the entire pipeline, source the `exec.R` file:

```r

source('exec.R')

```

## Terms of use

The pipeline results will be included in a future publication. To reference and use analysis results, please cite our GitHub repository and, if available, the publication. 

## Contact

Data and analysis requests should be addressed to [Dr. Renate Pichler](mailto:renate.pichler@i-med.ac.at). [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com) is the maintainer of the repository.

## References

[^1]: Sweeney SM, Cerami E, Baras A, Pugh TJ, Schultz N, Stricker T, Lindsay J, Del Vecchio Fitz C, Kumari P, Micheel C, et al. AACR Project GENIE: Powering Precision Medicine through an International Consortium. Cancer Discov (2017) 7:818–831. doi:10.1158/2159-8290.CD-17-0151
[^2]: Clinton TN, Chen Z, Wise H, Lenis AT, Chavan S, Donoghue MTA, Almassi N, Chu CE, Dason S, Rao P, et al. Genomic heterogeneity as a barrier to precision oncology in urothelial cancer. Cell Rep (2022) 41: doi:10.1016/j.celrep.2022.111859
[^3]: Robertson AG, Kim J, Al-Ahmadie H, Bellmunt J, Guo G, Cherniack AD, Hinoue T, Laird PW, Hoadley KA, Akbani R, et al. Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell (2017) 171:540-556.e25. doi:10.1016/J.CELL.2017.09.007
[^4]: Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, Kovatich AJ, Benz CC, Levine DA, Lee A V., et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell (2018) 173:400-416.e11. doi:10.1016/J.CELL.2018.02.052
[^5]: Balar A V., Galsky MD, Rosenberg JE, Powles T, Petrylak DP, Bellmunt J, Loriot Y, Necchi A, Hoffman-Censits J, Perez-Gracia JL, et al. Atezolizumab as first-line therapy in cisplatin-ineligible patients with locally advanced and metastatic urothelial carcinoma: a single-arm, multicentre, phase 2 trial. Lancet (London, England) (2017) 389:67. doi:10.1016/S0140-6736(16)32455-2
[^6]: Damrauer JS, Beckabir W, Klomp J, Zhou M, Plimack ER, Galsky MD, Grivas P, Hahn NM, O’Donnell PH, Iyer G, et al. Collaborative study from the Bladder Cancer Advocacy Network for the genomic analysis of metastatic urothelial cancer. Nat Commun (2022) 13: doi:10.1038/S41467-022-33980-9
[^7]: Kamoun A, de Reyniès A, Allory Y, Sjödahl G, Robertson AG, Seiler R, Hoadley KA, Groeneveld CS, Al-Ahmadie H, Choi W, et al. A Consensus Molecular Classification of Muscle-invasive Bladder Cancer. Eur Urol (2020) 77:420. doi:10.1016/J.EURURO.2019.09.006

