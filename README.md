# TSD-Analysis
Data and analysis of wild, natural sex ratio records of reptiles to investigate the origins of temperature-dependent sex determination. This analysis was conducted in RStudio (version 2022.02.3) using the MCMCglmm package. 


**System Requirements** 

The only system requirements is the program R

**Installation**

To install the MCMCglmm package, run the following code in R. 
```
install.packages(MCMCglmm)
library(MCMCglmm)
```
Instructions on how to use the MCMCglmm package can be found here: https://cran.r-project.org/web/packages/MCMCglmm/MCMCglmm.pdf

**Instructions for Use** 

Included in this package is 
1) a unique data set of wild, natural sex ratio records of reptiles
2) a phylogeny of species included in the dataset
3) Code to reformat data, run MCMCglmm models, and model selection for the sex ratio data set
4) Citation Code document that links citations in teh data set to their references

To  run this code, download the data set and phylogeny and save to your working directory. Then run the code in order to first reformat the data into a binary data set, run that data set through each model, identify the best model through model selection. Running each model takes a handful of hours.



**Sex Ratio Dataset Description**

This is a unique data set of 682 sex ratio records representing 182 species of reptiles from 601 populations. The data set was assembled from previously published literature of Bokony et al. (2019), with records from ROSIE, the Reptilian Offspring Sex and Incubation Environment (Krueger and Janzen 2022; ROSIE, 2021; v1.0.0), and personal communications. Inclusion criteria for each sex ratio datum were 1) recorded from a natural, unmanipulated population using 2) non-sex biased capture methods and 3) determined by reliable sexing methods with 4) a recorded sample size and 5) a life stage of sexed individuals.

For each record, the data consists of:

-  a species name 
- a population ID
- taxon: lizard, croc, turtle, snake
- their citation or citation code
- life stage: hatchling, juvenile adult
- sex ratio as proportion male
- supertaxa: whether the species belongs to a sister clade of crocodilians and turtles ‘crocoturtle’ or squamates and tuatara ’squamatara to test phylogenetic history model
- life history: whether the species belongs to crocodilians, turtles, and the tuatara ‘CTT’ that exhibit long-lived, late maturing life history, or squamates ’S’ that have short-lived, early maturing life histories to test the life history model
- sex determination: whether the species has genotypic sex determination (GSD) or temperature-dependent sex determination (TSD)
- SDM type: type of sex determining mechanism of GSD, TSD Ia, TSD Ib, TSD II
- sex determination source: reference for sex determining mechanism for the species
- Notes


**Phylogeny Description**

This phylogeny was created based on the phylogeny used in Bókony et al. (2019). Any species in this sex ratio data set not in the Bókony et al. (2019) were included using phytools package in R with the placement of each species was based on the same phylogenies Bókony et al. (2019) used to construct their phylogeny (Oaks 2011, Guillon et al. 2012, Pyron et al. 2013)


**TSD Analysis Code Description**

Code that takes the sex ratio data set and reformats it into the necessary binary format, then with the phylogeny, runs 4 MCMCglmm models to test competing hypothesis for evolution of temperature-dependent sex determination using model comparison

**Citation Code Description**

In the sex ratio data set, records from the Bokony et al (2019) data have citation codes rather than direct references in the dataset. This spreadsheet matches each citation code to its corresponding reference
