# scLASER

[![R-CMD-check](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**scLASER** 

<img src="scLASER_hex_cropped.png" align="right" width="200" alt="scLASER logo"/>



scLASER is an R package that provides a comprehensive framework for analyzing single-cell longitudinal data using neighbor abundance derived principal components (NAM-PCs). It integrates harmonization and mixed-effects modeling to capture dynamic associations between cellular composition, disease status throught time. The package also includes tools to simulate realistic longitudinal single-cell datasets and perform power calculations for study design, enabling users to benchmark analytical approaches and evaluate statistical robustness in complex temporal single-cell omics studies.

<br/>

## Installation

To install the latest development version directly from GitHub:

```r

devtools::install_github("fanzhanglab/scLASER")

```

<br/>

### Dependencies / Other required packages

```r
- R (>= 4.1.0)
- methods
- stats
- utils
- Matrix
- nlme
- lme4
- dplyr
- ggplot2
- ggrepel
- glue
- pbapply
- tidyr
- purrr
- caret
- xgboost
- uwot
- variancePartition
- edgeR
- Seurat
- RSpectra
- harmony
- MOFA2
- MatrixEQTL
- broom.mixed
- foreach
- doParallel
- pheatmap
- presto
- moments
- stevedata
- stevemisc
- shapviz
- rmarkdown
- knitr
```
> **Note:** Only the core dependencies (`Matrix`, `nlme`, `lme4`, `dplyr`, `ggplot2`) are required for standard use.  
> The remaining packages are suggested for optional analyses, visualization, or vignettes.

<br/>

## Tutorials

- [Simulated data (3 timepoints) tutorial](vignettes/Simulated_data_3timepoint_pipeline_tutorial.nb.html)
- [Simulated data (2 timepoints) tutorial](vignettes/Simulated_data_2timepoint_pipeline_tutorial.nb.html)
- [General data simulation overview](vignettes/Data_simulation.nb.html)


<br/>

## Citations

**Authors:**  
Lauren Vanderlinden (lauren.vanderlinden@cuanschutz.edu)  
Juan Vargas (juan.vargas@cuanschutz.edu)  
Fan Zhang (fanzhanglab@gmail.com)


<br/>

## Help, Suggestion and Contribution

Using github [**issues**](https://github.com/fanzhanglab/scLaser/issues)
section, if you have any question, comments, suggestions, or to report
coding related issues of scLASER is highly encouraged than sending
emails.

- Please **check the GitHub
  [issues](https://github.com/fanzhanglab/scLASER/issues)** for similar
  issues that has been reported and resolved. This helps the team to
  focus on adding new features and working on cool projects instead of
  resolving the same issues!
- **Examples** are required when filing a GitHub issue. In certain
  cases, please share your scLASER object and related codes to understand
  the issues.
  
<br/>

## Contact 

Please contact [fanzhanglab@gmail.com](fanzhanglab@gmail.com) for further questions or potential collaborative opportunities!
