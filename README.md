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

## Identifying SV genes within groups
### Data
To help illustrate how SPADE package can be applied, the SeqFISH dataset was provided in the package with 249 genes measured on 131 spots.
```r
library(SPADE)
data(SeqFISH)
data(info)
readcounts <- SeqFISH
dim(readcounts)
249 131
dim(info)
131   2

readcounts[1:5,1:5]
## C1 C2 C3 C4 C5
## Tal1 3 3 17 20 6
## Dmbx1 11 9 1 4 0
## Emx2 13 14 2 9 0
## Uncx 5 21 3 5 4
## Paxip1 15 9 10 13 7

head(info)
## x y
## C1 686 352
## C2 488 572
## C3 614 698
## C4 516 726
## C5 308 208
## C6 204 225
```
### Normalization
A two-step normalization strategy was implemented to transformed original read counts data into continuous data.
```r
data_norm <- SPADE_norm(readcounts=as.matrix(readcounts), info=info)
```
### Parameter estimation
For each gene, we estimated the optimal length-scale hyperparameter in the Gaussian kernel to increase the accuracy of SV gene identification.
```r
Est <- SPADE_estimate(expr_data=data_norm, info=info)
head(Est)
 
  GeneID theta_Gau  Gamma_hat   Lik_Gau
1      1  11.22829  0.8796637 -91.34812
2      2  32.89034  2.2457761 -87.74524
3      3 130.69611 19.9999373 -83.34457
4      4  29.62392  4.3757445 -77.91255
5      5  44.91683  2.2771921 -67.52400
6      6  58.66151  3.1931255 -63.95050
```
### Testing
After the optimal length-scale hyperparameter was estimated, P-value for each gene was computed based on a quadratic score statistic with a Davies method.
```r
Test_res <- SPADE_test(object=data_norm, location=info, para=Est)
Test_res[c(1, 230),]
## geneid Q Pvalue Adjust.Pvalue
    geneid         Q    Pvalue Adjust.Pvalue
1     Tal1  65.77883 0.4434016     0.5407541
230   lyve 487.63637 0.0000000     0.0000000
```

## Identifying SV genes between groups
### Data
To illustrate how SPADE can be applied to identify SV genes between groups, a real spatial transcriptomic dataset with axolotl telencephalon (i.e., ARRISTA) was provided in the package. Two different post injury stages (i.e., 2DPI, 5DPI) were included in the dataset. This data contain three genes with 1,188 spots and 938 spots in the 2DPI group and 5DPI group, respectively.
```r
data(D2_data)
data(D2_info)
data(D5_data)
data(D5_info)
dim(D2_data)
## [1] 3 1188
dim(D2_info)
## [1] 1188 2
dim(D5_data)
## [1] 3 938
dim(D5_info)
## [1] 938 2
```
### Normalization
Still, the same two-step normalization strategy was implemented for read counts data from each group to transform original data into continuous data.
```r
D2_norm <- SPADE_norm(readcounts=as.matrix(D2_data), info=D2_info)
D5_norm <- SPADE_norm(readcounts=as.matrix(D5_data), info=D5_info)
```
### Parameter estimation and testing
SPADE identifies SV genes between groups based on a crossed likelihood ratio test in spatial transcriptomic data. SPADE first estimates the optimal hyperparameter for kernel matrix in each group, respectively. Thus, for each gene, the log likelihood in each group can be easily calculated with its optimal kernel. Then we exchange the estimated hyperparameters to compute the log likelihoods for both groups, and compare them to their optimal log likelihoods. The likelihood ratio test statistic is calculated to identify SV genes with P-values computed using F test with degree freedom of one.
```r
res <- SPADE_DE(D2_norm, D5_norm, D2_info, D5_info)
res
              geneid theta_Gau1 theta_Gau2           Gamma1           Gamma2   logLik11  logLik21   logLik10
1 AMEX60DDU001010113   49.52013   70.40230  2.0452367735223 1.52640627368895   345.6043  343.4853   279.5261
2 AMEX60DDU001038720   49.63787   81.30095 19.9998220761637 1.33625829767921 -1407.1231 -574.7335 -1563.7369
3 AMEX60DDU001022818   49.63787   50.10993 4.19220647492374   19.99993729622   143.7086  329.6725   143.6045
   logLik20       Diff       Pvalue Adjust.Pvalue
1  327.4874 164.152145 1.401280e-37  2.101920e-37
2 -642.6539 449.068302 1.150481e-99  3.451442e-99
3  329.1845   1.184317 2.764789e-01  2.764789e-01
```

## Citation
 @Article{,
    author = {Fei Qin and Guoshuai Cai},
    title = {Spatial pattern and differential expression analysis with spatial transcriptomic data},
    doi = {10.5281/zenodo.13313887},
  }
