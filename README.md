# Comparing statistical learning models in wastewater-based epidemiology: An application to norovirus

## Overview
This repository contains the data and code used in the paper "Comparing statistical learning models in wastewater-based epidemiology: An application to norovirus". We applied the following statistical and machine learning methods to predict pathogen concentration in wastewater, using norovirus in England as a case study:
- Generalised Additive Models implemented in the mgcv package
- Bayesian Generalised Additive Models implemented in brms
- LASSO regression implemented in
- Random Forest implemented in ranger

We compared the above approaches to our benchmark model which used the framework outlined in Li et al. (2023) using the Integrated Nested Laplace Approximations (INLA) coupled with Stochastic Partial Differential Equations (SPDE) framework, implemented in the R-INLA package. This allows probablistic predictions of pathogen concentrations in wastewater

The data used in this research was collected as part of the Environmental Monitoring for Health Protection programme in England. 

