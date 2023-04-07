# SPADE
A spatial pattern investigation method to identify spatially variable genes using spatial transcriptomic data.

## Author
Fei Qin, Feifei Xiao, Guoshuai Cai

## Description
To identify spatially varaible (SV) genes with spatial transcriptomic data, the SPADE method was developed based on a Gaussian process regression (GPR) model, which can model the relationship between gene expression and other covariates (i.e., cell groups) incorporating the spatial information of multiple cells. SPADE was able to identify SV genes within a single group and between groups. First, original read counts data were normalized into continuous data. Second, instead of using a fixed length scale hyperparameter in the covariance kernel of GPR, SPADE estimated the optimal hyperparameter for each gene. To identify SV genes within groups, hypothesis testing was conducted based on a quadratic score statistic with a Davies method to compute the P value. With SV gene detection between groups, SPADE exchanged the optimal hyperparameters estimated in two groups and then utilized a crossed likelihood ratio test to calculate the P value for each gene. 

## Installation
```r
library(devtools)
install_github("thecailab/SPADE")
```
