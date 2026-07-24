# Comparing statistical learning models in wastewater-based epidemiology: An application to norovirus

## Overview
This repository contains the data and code used in the paper "Comparing statistical learning models in wastewater-based epidemiology: An application to norovirus". We applied the following statistical and machine learning methods to predict pathogen concentration in wastewater, using norovirus in England as a case study:
- Generalised Additive Models implemented in the [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html)
- Bayesian Generalised Additive Models implemented in [brms](https://cran.r-project.org/web/packages/brms/index.html)
- LASSO regression implemented in [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
- Random Forest implemented in [ranger](https://cran.r-project.org/web/packages/ranger/index.html)

We compared the above approaches to our benchmark model which used the framework outlined in [Li et al. (2023)](https://www.sciencedirect.com/science/article/pii/S0160412023000387?via%3Dihub) using the Integrated Nested Laplace Approximations (INLA) coupled with Stochastic Partial Differential Equations (SPDE) framework, implemented in the [R-INLA](https://www.r-inla.org/) package. This allows probablistic predictions of pathogen concentrations in wastewater

The data used in this research was collected as part of the Environmental Monitoring for Health Protection programme in England. 

This repository contains the data and scripts that can be used to replicate the results of the paper:
-[cross-validation](cross-validation) contains the scripts used for 10-fold sptial cross validation of all the models
- [data](data) contains the data used for model fitting and prediction at the LSOA centroids
- [plots](plots) contains scripts used to generate the plots that were included in the paper
- [prediction](prediction) contains the script to generate teh norovirus wastewater predictions at the LSOA centroids for the INLA-SPDE and Random Forest models


