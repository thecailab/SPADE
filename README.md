# SPADE
A spatial pattern investigation method to identify spatially variable genes using spatial transcriptomic data.

## Author
Fei Qin, Feifei Xiao, Guoshuai Cai

## Description
To identify spatially varaible (SV) genes with spatial transcriptomic data, the SPADE method was developed based on a Gaussian process regression (GPR) model, which can model the relationship between gene expression and other covariates (i.e., cell groups) incorporating the spatial information of multiple cells. SPADE was able to identify SV genes within a single group and between groups. First, original read counts data were normalized into continuous data. Second, instead of using a fixed length scale hyperparameter in the covariance kernel of GPR, SPADE estimated the optimal hyperparameter for each gene. To identify SV genes within groups, hypothesis testing was conducted based on a quadratic score statistic with a Davies method to compute the P value. With SV gene detection between groups, SPADE exchanged the optimal hyperparameters estimated in two groups and then utilized a crossed likelihood ratio test to calculate the P value for each gene. 
![Figure 1](https://user-images.githubusercontent.com/68352557/234379105-0d4b88f1-0bf0-4525-9244-610b6cabcc6d.png)

## Installation
```r
library(devtools)
install_github("thecailab/SPADE")
```

## Identifying SV genes within groups
### Data
To help illustrate how SPADE package can be applied, the SeqFISH dataset was provided in the package with 249 genes measured on 131 spots.
```r
library(SPADE)
data(SeqFISH)
data(info)
readcounts <- SeqFISH
dim(readcounts)
dim(info)
```
[1] 249 131
[1] 131   2




