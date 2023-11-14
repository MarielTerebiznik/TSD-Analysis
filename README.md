# TSD-Analysis
Data and analysis of wild, natural sex ratio records of reptiles to investigate the origins of temperature-dependent sex determination. This analysis was conducted in RStudio (version 2022.02.3) using the MCMCglmm package. 


**System Requirements** /n
The only system requirements is the program R

**Installation** /n
To install the MCMCglmm package, run the following code in R. 
```
install.packages(MCMCglmm)
library(MCMCglmm)
```
Instructions on how to use the MCMCglmm package can be found here: https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf

**Instructions for Use** /n
Included in this package is 
1) a unique data set of wild, natural sex ratio records of reptiles
2) a phylogeny of species included in the dataset
3) Code to reformat data, run MCMCglmm models, and model selection for the sex ratio data set

To  run this code, download the data set and phylogeny and save to your working directory. Then run the code in order to first reformat the data into a binary data set, run that data set through each model, identify the best model through model selection.
