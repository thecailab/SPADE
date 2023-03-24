#' @title SPADE estimation and test between groups
#'
#' @description Estimation of optimal parameter in the kernel function to identify SE genes between groups. P values will be provided based on a crossed likelihood ratio test. 
#'
#' @param readcounts1 Read counts data for identifying spatially expressed genes in group 1. Each row is a gene and each column is a spot.
#' @param readcounts2 Read counts data for identifying spatially expressed genes in group 2. Each row is a gene and each column is a spot.
#' @param location1 Spatial coordinates for all spots in group 1.
#' @param location2 Spatial coordinates for all spots in group 2.
#' 
#'
#' @return This function returns the estimated parameters and some statistics of the SPADE method.
#' \item{geneid}{Gene index.}
#' \item{theta_Gau1}{The estimated optimal length-scale hyperparameter in the Gaussian kernel for group 1.}
#' \item{theta_Gau2}{The estimated optimal length-scale hyperparameter in the Gaussian kernel for group 2.}
#' \item{logLik11}{The corresponding log likelihood calculated with optimal hyperparameter estimated above for group 1.}
#' \item{logLik21}{The corresponding log likelihood calculated with optimal hyperparameter estimated above for group 2.}
#' \item{logLik10}{The log likelihood calculated with expression value from group 1 and optimal hyperparameter estimated for group 2.}
#' \item{logLik20}{The log likelihood calculated with expression value from group 2 and optimal hyperparameter estimated for group 1.}
#' \item{Diff}{The likelihood ratio test statistic is given by 2*(logLik11 + logLik21 - logLik10 - logLik20).}
#' \item{Pvalue}{P values calculated with likelihood ratio test statistic using F test with degree freedom of one.}
#' \item{Adjust.Pvalue}{Adjusted P values calculated with the Benjamini and Hochberg method.}
#'
#' @export
SPADE_DE <- function(readcounts1, readcounts2, location1, location2){ 
  
  ED1 <- as.matrix(dist(location1))
  lrang1 <- SPARK::ComputeGaussianPL(ED1, compute_distance=FALSE)
  
  ED2 <- as.matrix(dist(location2))
  lrang2 <- SPARK::ComputeGaussianPL(ED2, compute_distance=FALSE)
  
  # ED2 = as.matrix(dist(rbind(location, location)))
  
  para <- matrix(NA, nrow(readcounts1), 7)
  colnames(para) <- c("geneid","theta_Gau1", "theta_Gau2",
                      "logLik11","logLik21", "logLik10", "logLik20")
  for (i_gene in 1:nrow(readcounts1)){
    cat(paste("NO. Gene = ",i_gene,"\n"))
    y1 <- readcounts1[i_gene,]
    y2 <- readcounts2[i_gene,]
    
    re1 <- optimize(lengthscale_fit, c(lrang1[3],lrang1[8]), location=location1, y=y1, tol=1)
    re2 <- optimize(lengthscale_fit, c(lrang2[3],lrang2[8]), location=location2, y=y2, tol=1)
    
    Logdelta1 <- Delta_fit(location=location1, y=y1, L=re1$minimum)
    Logdelta2 <- Delta_fit(location=location2, y=y2, L=re2$minimum)
    
    LL10 <- LL_DE(log_delta=Logdelta1, location=location1, L=re2$minimum, y=y1)
    LL20 <- LL_DE(log_delta=Logdelta2, location=location2, L=re1$minimum, y=y2)
    
    para[i_gene,1] <- rownames(readcounts1)[i_gene] 
    
    para[i_gene,2] <- re1$minimum
    para[i_gene,3] <- re2$minimum
    # para[i_gene,4] <- re12$minimum
    
    para[i_gene,4] <- -re1$objective
    para[i_gene,5] <- -re2$objective
    
    para[i_gene,6] <- LL10
    para[i_gene,7] <- LL20
    # para[i_gene,7] <- -re12$objective
    
  }
  
  para <-  as.data.frame(para)
  para$theta_Gau1 <- as.numeric(para$theta_Gau1)
  para$theta_Gau2 <- as.numeric(para$theta_Gau2)
  
  para$logLik11 <- as.numeric(para$logLik11)
  para$logLik21 <- as.numeric(para$logLik21)
  
  para$logLik10 <- as.numeric(para$logLik10)
  para$logLik20 <- as.numeric(para$logLik20)
  
  para$Diff <- 2*(para$logLik11 + para$logLik21 - para$logLik10 - para$logLik20)
  para$Pvalue <- pchisq((para$Diff),df=1,lower.tail=FALSE)
  para$Adjust.Pvalue <- p.adjust(para$Pvalue, method = "BH")
  
  return(para)
}




#' @title Optimization of likelihood function between groups
#'
#' @description Optimization of likelihood function between groups.
#'
#' @param location1 Spatial coordinates for all spots in group 1.
#' @param location2 Spatial coordinates for all spots in group 2.
#' @param y1 Read counts data for the gene in group 1.
#' @param y2 Read counts data for the corresponding gene in group 2.
#' @param L The length-scale hyperparameter in the kernel function.
#' 
#'
#' @return This function returns likelihood at the point of optimal delta.
#' \item{results}{The likelihood with the optimal delta value for identifying SE genes between groups.}
#' 
#'
#' @export
lengthscale_fit_DE <- function(location1, location2, y1, y2, L) {
  R1 <- as.matrix(dist(location1) ** 2)
  K1 <- exp(-R1 / (2 * L ** 2))
  
  R2 <- as.matrix(dist(location2) ** 2)
  K2 <- exp(-R2 / (2 * L ** 2))
  
  US1 <- eigen(K1)
  US1$values[US1$values < 0.] <- 0.
  U1=US1$vectors
  S1=US1$values
  S1[S1 < 0.] <- 0.
  
  US2 <- eigen(K2)
  US2$values[US2$values < 0.] <- 0.
  U2=US2$vectors
  S2=US2$values
  S2[S2 < 0.] <- 0.
  
  UTy1 <- t(y1) %*% U1
  UTy2 <- t(y2) %*% U2
  
  UT1 <- colSums(U1)
  UT2 <- colSums(U2)
  
  n1 <- length(y1)
  n2 <- length(y2)
  ## Find optimal delta value
  o <- optimize(LL_combin, c(-10, 20), UTy1=UTy1, UTy2=UTy2, UT1=UT1, UT2=UT2, S1=S1, S2=S2, n1=n1, n2=n2)$objective
  return(o)
}





#' @title Log likelihood function between groups
#'
#' @description The log transformed multivariate normally distributed marginal likelihood to identify SE genes between groups. 
#'
#' @param log_delta Log transformed delta in the variance function.
#' @param UTy1 The expression UTy will need to be re-computed for the gene in group 1.
#' @param UTy2 The expression UTy will need to be re-computed for the gene in group 2.
#' @param UT1 The expression UT1 only depends on the coordinates X and can be pre-computed and reused for the gene in group 1.
#' @param UT2 The expression UT1 only depends on the coordinates X and can be pre-computed and reused for the gene in group 2.
#' @param S1 A diagonal matrix used in the spectral decomposition of the covariance in group 1.
#' @param S2 A diagonal matrix used in the spectral decomposition of the covariance in group 2.
#' @param n1 The number of spots for each gene in group 1. 
#' @param n2 The number of spots for each gene in group 2. 
#'
#' @return This function returns the log likelihood for each delta.
#' \item{LL_res}{The log likelihood for each delta value.}
#' 
#'
#' @export
LL_combin <- function(log_delta, UTy1, UTy2, UT1, UT2, S1, S2, n1, n2){
  ll_diff <- LL(exp(log_delta), UTy1, UT1, S1, n1)+
    LL(exp(log_delta), UTy2, UT2, S2, n2) 
  return(ll_diff)
}






#' @title Estimation of the optimal delta
#'
#' @description Estimation of the optimal delta in the variance function between groups.
#'
#' @param location Spatial coordinates for all spots.
#' @param y Read counts data for the gene.
#' @param L The length-scale hyperparameter in the kernel function.
#' 
#'
#' @return This function returns the optimal log transformed delta value.
#' \item{results}{The optimal log delta value.}
#' 
#'
#' @export
Delta_fit <- function(location, y, L){
  R2 <- as.matrix(dist(location) ** 2)
  K <- exp(-R2 / (2 * L ** 2))
  
  US <- eigen(K)
  US$values[US$values < 0.] <- 0.
  U=US$vectors
  S=US$values
  S[S < 0.] <- 0.
  
  UTy <- t(y) %*% U
  UT1 <- colSums(U)
  n <- length(y)
  
  ## Find optimal delta value
  results <- optimize(LL, c(-10, 20), UTy=UTy, UT1=UT1, S=S, n=n)$minimum
  return(results)
}



#' @title Log likelihood function
#'
#' @description The log transformed multivariate normally distributed marginal likelihood . 
#'
#' @param log_delta Log transformed delta in the variance function.
#' @param location Spatial coordinates for all spots.
#' @param y Read counts data for the gene.
#' @param L The length-scale hyperparameter in the kernel function.
#'
#' @return This function returns the log likelihood for each delta.
#' \item{LL_res}{The log likelihood for each delta value.}
#' 
#'
#' @export
LL_DE <- function(log_delta, location, y, L) {
  R2 <- as.matrix(dist(location) ** 2)
  K <- exp(-R2 / (2 * L ** 2))
  
  US <- eigen(K)
  US$values[US$values < 0.] <- 0.
  U=US$vectors
  S=US$values
  S[S < 0.] <- 0.
  
  UTy <- t(y) %*% U
  UT1 <- colSums(U)
  n <- length(y)
  
  
  delta <- exp(log_delta)
  mu_h <- mu_hat(delta, UTy, UT1, S)
  
  sum_1 <- sum(((UTy - UT1 * mu_h)^2) / (S + delta))
  sum_2 <- sum(log(S + delta))
  
  LL_res <- -0.5 * (n * log(2 * pi) + n * log(sum_1 / n) + sum_2 + n)
  
  return(LL_res)
}







## Combine mode: Estimating parameters for each gene in SPDE
# SPDE_estimate_DE <- function(readcounts1, readcounts2, location1, location2){ 
#   
#   ED1 <- as.matrix(dist(location1))
#   lrang1 <- SPARK::ComputeGaussianPL(ED1, compute_distance=FALSE)
#   
#   ED2 <- as.matrix(dist(location2))
#   lrang2 <- SPARK::ComputeGaussianPL(ED2, compute_distance=FALSE)
#   
#   para <- matrix(NA, nrow(readcounts1), 7)
#   colnames(para) <- c("geneid","theta_Gau1", "theta_Gau2","theta_Gau12",
#                       "logLik1","logLik2", "logLik12")
#   for (i_gene in 1:nrow(readcounts1)){
#     cat(paste("NO. Gene = ",i_gene,"\n"))
#     y1 <- readcounts1[i_gene,]
#     y2 <- readcounts2[i_gene,]
#     
#     re12 <- optimize(lengthscale_fit_DE, c(lrang1[3],lrang1[8]), location1=location1, location2=location2, y1=y1, y2=y2, tol=1)
#   
#     re1 <- optimize(lengthscale_fit, c(lrang1[3],lrang1[8]), location=location1, y=y1, tol=1)
#     re2 <- optimize(lengthscale_fit, c(lrang2[3],lrang2[8]), location=location2, y=y2, tol=1)
# 
#     para[i_gene,1] <- rownames(readcounts1)[i_gene] 
#     
#     para[i_gene,2] <- re1$minimum
#     para[i_gene,3] <- re2$minimum
#     para[i_gene,4] <- re12$minimum
#     
#     para[i_gene,5] <- -re1$objective
#     para[i_gene,6] <- -re2$objective
#     para[i_gene,7] <- -re12$objective
#     
#   }
#   
#     para <-  as.data.frame(para)
#     para$theta_Gau1 <- as.numeric(para$theta_Gau1)
#     para$theta_Gau2 <- as.numeric(para$theta_Gau2)
#     para$theta_Gau12 <- as.numeric(para$theta_Gau12)
#   
#     para$logLik1 <- as.numeric(para$logLik1)
#     para$logLik2 <- as.numeric(para$logLik2)
#     para$logLik12 <- as.numeric(para$logLik12)
#   
#     para$Diff <- para$logLik1 + para$logLik2 - para$logLik12 
#     para$Pvalue <- pchisq(2*(para$Diff),df=2,lower.tail=FALSE)
#     para$Adjust.Pvalue <- p.adjust(para$Pvalue, method = "BH")
#     
#   return(para)
# }
