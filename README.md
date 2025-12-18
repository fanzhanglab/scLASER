

[![R-CMD-check](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fanzhanglab/scLASER/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Visitors](https://visitor-badge.laobi.icu/badge?page_id=fanzhanglab.scLASER)




<h1 style="font-size: 38px;">scLASER</h1>

<img src="scLASER_hex.png" align="right" width="220" alt="scLASER logo"/>

<p style="font-size: 18px; line-height: 1.6;">
Longitudinal single-cell clinical studies enable tracking within-individual cellular dynamics, but methods for modeling temporal phenotypic changes and estimating power remain limited. We present <strong>scLASER</strong>, a framework <strong>detecting time-dependent cellular neighborhood dynamics</strong> and <strong>simulating longitudinal single-cell datasets for power estimation</strong>. Across benchmarks, scLASER shows superior sensitivity, particularly for non-linear temporal patterns. Applications to inflammatory bowel disease reveal treatment-responsive <strong>NOTCH3+ stromal trajectories</strong>, while analysis of COVID-19 data identifies distinct axes of <strong>T cell activity</strong> over disease progression. scLASER enables robust longitudinal single-cell analysis and optimization of study design.
</p>

<br/>

<img width="100%" align="center" src="https://github.com/fanzhanglab/scLASER/blob/main/docs/figure.png?raw=true">

**What does scLASER do?**  
- ***Simulate*** multi-timepoint longitudinal single-cell datasets with distinct dynamic patterns for clinical outcome (e.g., treatment response). 
- ***Detect*** time-dependent cellular dynamics (linear and nonlinear) associated with treatment response.  
- ***Generate*** a per-cell association score quantifying each cell's contribution to time 
