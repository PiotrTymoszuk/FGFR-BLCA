# FGFR-BLCA
Genetic alterations and expression of genes coding for FGF ligands and FGF reseptors in urothelial cancer

## Results

Figures with the analysis results are available [here](https://github.com/PiotrTymoszuk/FGFR-BLCA/tree/main/report/figures).
Analysis report is available from [a Drobbox folder](https://www.dropbox.com/scl/fo/q79r54u1ng0g31wwam6yl/AJm-THPjz8EuYVpi-y9k1dA?rlkey=lxrqv5vqifspdh27crys7ubpc&dl=0).

## Usage

The analysis pipeline requires some development packages, the easiest way to install them is to use `devtools`:

```r

devtools::install_github('PiotrTymoszuk/trafo')
devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/soucer')
devtools::install_github('PiotrTymoszuk/figur')
devtools::install_github('PiotrTymoszuk/microViz')

```
To launch the entire pipeline, source the `exec.R` file:

```r

source('exec.R')

```

## Terms of use

The pipeline results will be included in a future publication. To reference and use analysis results, please cite our GitHub repository and, if available, the publication. 

## Contact

Data and analysis requests should be addressed to [Dr. Renate Pichler](mailto:renate.pichler@i-med.ac.at). [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com) is the maintainer of the repository.
