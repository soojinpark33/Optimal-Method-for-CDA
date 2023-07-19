# Optimal-Method-for-CDA

Soojin Park<sup>1</sup>, Suyeon Kang<sup>2</sup>, and Chioun Lee<sup>3</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Statistics, University of California, Riverside
<sup>3</sup> Department of Sociology, University of California, Riverside


## Overview

Causal decomposition analysis is among the rapidly growing number of tools for identifying factors (`mediators') that contribute to disparities in outcomes between social groups. An example of such mediators is college completion, which explains later health disparities between Black women and White men. The goal is to quantify how much a disparity would be reduced (or remain) if we hypothetically intervened to set the mediator distribution equal across social groups. Despite increasing interest in estimating disparity reduction and remaining, various estimation procedures are not straightforward and researchers have scant guidance for choosing an optimal method. In this article, we evaluate the performance in terms of bias, variance, and coverage of three approaches that utilize different modeling strategies: 1) regression-based methods that impose restrictive modeling assumptions (e.g., linearity) and 2) weighting-based and 3) imputation-based methods that rely on the observed distribution of variables. We find a trade-off between the modeling assumptions required in the method and its performance. In terms of performance, regression-based methods operate best as long as the restrictive assumption of linearity is met. Methods relying on mediator models without imposing any modeling assumptions are sensitive to the ratio of the group-mediator association to the mediator-outcome association. These results highlight the importance of selecting an appropriate estimation procedure considering the data at hand.

For more details of our proposed methods, see [our paper]([https://journals.sagepub.com/doi/10.1177/00811750231183711]). 
Here, we provide `R` codes to reproduce our simulation study. 


## Simulation Study

* `SM_Section4_Supp.R`  

   This `R` file contains the simulation codes for different estimators of causal decomposition analysis. This code replicates our results in Figures 1 and 2 as well as Supplementary Appendices C and D of our paper.

* `SM_source.R` 
 
   This `R` file includes source functions required to run our simulation codes. 

These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.

