# scLASER

[![R-CMD-check](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


**scLASER** 

![scLASER logo](scLASER_hex_cropped.png)

is an R package for modeling, simulating, and visualizing single-cell longitudinal datasets.

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


- [Simulated data (3 timepoints) tutorial](vignettes/Simulated_data_3timepoint_pipeline_tutorial.html)


<br/>

## Citations

Authors


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
