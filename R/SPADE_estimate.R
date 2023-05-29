#' @title SPADE estimation for hyperparameter within groups
#'
#' @description Estimation of optimal hyperparameter in the kernel function to identify SE genes within groups.
#'
#' @param expr_data Read counts data for identifying spatially expressed genes. Each row is a gene and each column is a spot.
#' @param info Spatial coordinates for all spots.
#'
#'
#' @return This function returns the estimated parameters and some statistics of the SPADE method.
#' \item{GeneID}{Gene index.}
#' \item{theta_Gau}{The estimated optimal length-scale hyperparameter in the Gaussian kernel}
#' \item{Gamma_hat}{The estimated optimal Gamma_1 in the covariance function.}
#' \item{Lik_Gau}{The log likelihood using optimal hyparameter estimated above.}
#' 
#'
#' @examples
#' Y <- matrix(rnorm(10000, 10, 2),100, 100)
#' info <- matrix(runif(200, 1, 100), 100, 2)
#' output <- SPADE_estimate(expr_data=Y, info=)
#' output
#'
#' @export
SPADE_estimate <- function(expr_data, info){
  ED <- as.matrix(dist(info))
  lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)
  final <- NULL
  for (i in 1:nrow(expr_data)){
    print(paste0("Gene ", i))
    y = expr_data[i,]
    re_Gau <- optimize(lengthscale_fit, c(lrang[3],lrang[8]), location=info, y=y, tol=1)
    
    Gamma_hat <- Delta_fit(location=info, y=y, L=re_Gau$minimum)$Tao_hat
    
    # re_Per <- optimize(lengthscale_fit_Per, c(lrang[3],lrang[8]), location=info, y=y, tol=1)
    
    final1 <- data.frame(GeneID=i, theta_Gau=re_Gau$minimum, Gamma_hat=Gamma_hat, Lik_Gau=-re_Gau$objective)
    # final1 <- data.frame(GeneID=i, theta_Gau=re_Gau$minimum, Lik_Gau=-re_Gau$objective,
    #                                theta_Per=re_Per$minimum, Lik_Per=-re_Per$objective)
    final <- rbind(final, final1)
  }
  return(final)
}



#' @title Optimization of likelihood function within groups
#'
#' @description Optimization of likelihood function within groups.
#'
#' @param y Read counts data for a gene.
#' @param location Spatial coordinates for all spots.
#' @param L The length-scale hyperparameter in the kernel function.
#' 
#'
#' @return This function returns likelihood at the point of optimal delta.
#' \item{results}{The likelihood with the optimal delta value.}
#' 
#'
#' @export
lengthscale_fit <- function(location, y, L){
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
  results <- optimize(LL, c(-10, 20), UTy=UTy, UT1=UT1, S=S, n=n)$objective
  return(results)
}


#' @title Log likelihood function within groups
#'
#' @description The log transformed multivariate normally distributed marginal likelihood to identify SE genes within groups. 
#'
#' @param log_delta Log transformed delta in the variance function.
#' @param UTy The expression UTy will need to be re-computed for each gene.
#' @param UT1 The expression UT1 only depends on the coordinates X and can be pre-computed and reused for every gene.
#' @param S A diagonal matrix used in the spectral decomposition of the covariance.
#' @param n The number of spots for each gene. 
#'
#' @return This function returns the log likelihood for each delta.
#' \item{LL_res}{The log likelihood for each delta value.}
#' 
#'
#' @export
LL <- function(log_delta, UTy, UT1, S, n) {
  delta <- exp(log_delta)
  mu_h <- mu_hat(delta, UTy, UT1, S)
  
  sum_1 <- sum(((UTy - UT1 * mu_h)^2) / (S + delta))
  sum_2 <- sum(log(S + delta))
  
  LL_res <- -0.5 * (n * log(2 * pi) + n * log(sum_1 / n) + sum_2 + n)
  
  return(-LL_res)
}




#' @title The estimate of mu
#'
#' @description The estimate of mean expression level for each gene. 
#'
#' @param delta Delta value in the variance function.
#' @param UTy The expression UTy will need to be re-computed for each gene,
#' @param UT1 The expression UT1 only depends on the coordinates X and can be pre-computed and reused for each gene.
#' @param S A diagonal matrix used in the spectral decomposition of the covariance.
#'
#' @return This function returns the estimate of mean expression level for each gene.
#' \item{mu_h}{The estimate of mean expression level for each gene.}
#' 
#'
#' @export
mu_hat <- function(delta, UTy, UT1, S) {
  UT1_scaled <- UT1 / (S + delta)
  sum1 <- UT1_scaled %*% t(UTy)
  sum2 <- UT1_scaled %*% UT1
  mu_h <- sum1 / sum2
  return(mu_h)
}




#' @title Hyperparameter range.
#' @description The range for the hyperparameter in the Gaussian Kernel, which was calculated from the coordinate of all spots.
#' 
#' @param X Cell coordinates matrix n x 2 or kernel matrix computed already
#' @param compute_distance Compute the distance matrix using generic function dist, default=TRUE
#' 
#' @return This function returns the range for the hyperparameter in the Gaussian Kernel.
#' \item{lrang}{The range for the hyperparameter in the Gaussian kernel}
#' 
#' @export
ComputeGaussianPL <- function(X, compute_distance=TRUE){
  if(compute_distance){
    if(ncol(X)<2){stop("X has to be a coordinate matrix with number of column greater than 1")}
    D <- dist(X)
  }else{
    D <- X
  }# end fi
  
  #Dval <- unique(as.vector(D))
  Dval <- D
  Dval.nz <- Dval[Dval>1e-8]
  lmin <- min(Dval.nz)/2
  lmax <- max(Dval.nz)*2
  lrang <- 10^(seq(log10(lmin),log10(lmax),length.out=10))
  return(lrang)
}# end func




# lengthscale_fit_Gau <- function(location, y, L){
#   R2 <- as.matrix(dist(location) ** 2)
#   K <- exp(-R2 / (2 * L ** 2))
#   
#   US <- eigen(K)
#   US$values[US$values < 0.] <- 0.
#   U=US$vectors
#   S=US$values
#   S[S < 0.] <- 0.
#   
#   UTy <- t(y) %*% U
#   UT1 <- colSums(U)
#   n <- length(y)
#   
#   ## Find optimal delta value
#   results <- optimize(LL, c(-10, 20), UTy=UTy, UT1=UT1, S=S, n=n)$objective
#   return(results)
# }


# 
# lengthscale_fit_Per <- function(location, y, L){
#   R <- as.matrix(dist(location))
#   # K <- exp(-R2 / (2 * L ** 2))
#   K <- cos(2*pi*R/L)
#     
#   US <- eigen(K)
#   US$values[US$values < 0.] <- 0.
#   U=US$vectors
#   S=US$values
#   S[S < 0.] <- 0.
#   
#   UTy <- t(y) %*% U
#   UT1 <- colSums(U)
#   n <- length(y)
#   
#   ## Find optimal delta value
#   results <- optimize(LL, c(-10, 20), UTy=UTy, UT1=UT1, S=S, n=n)$objective
#   return(results)
# }
# 

