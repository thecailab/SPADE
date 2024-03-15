#' @title SPADE test within groups
#'
#' @description Test for identifying SE genes within groups
#'
#' @param object Read counts data for identifying spatially expressed genes. Each row is a gene and each column is a spot.
#' @param location Spatial coordinates for all spots.
#' @param para Parameter data frame estimated from the SPADE_estimate function.
#'
#' 
#'
#' @return This function returns the statistics and P values for the test
#' \item{GeneID}{Gene index.}
#' \item{Q}{The Q statistics calculated using the SKAT method.}
#' \item{Pvalue}{P value calculated from the Q statistics based on the Davies method.}
#' \item{Adjust.Pvalue}{Adjusted P values calculated with the Benjamini and Hochberg method.}
#' 
#'
#'
#' @export
#' @importFrom SKAT SKAT_Null_Model
#' @importFrom CompQuadForm davies
SPADE_test <- function(object, location, para){
  
  ED <- as.matrix(dist(location))
  R2 <- as.matrix(dist(location) ** 2)

  num_gene_test <- nrow(object)
  
  res_all <- NULL
  for (ig in 1:num_gene_test){
    if (ig %% 20==0) cat(ig,"..")

    y.c <- as.numeric(object[ig, ])
    
    K <- exp(-R2 / (2 * para$theta[ig] ** 2))
   
    obj<-SKAT::SKAT_Null_Model(y.c ~ 1, out_type="C")
    res <- obj$res
    s2 <- obj$s2
    X <- obj$X1
    
    Q = t(res) %*% K %*% res / (2*s2)
    W = K - X%*%solve( t(X)%*%X)%*%( t(X) %*% K)	# W = P0 K
    W1 = W - (W %*% X) %*% solve(t(X)%*%X) %*% t(X)
    K1 <- W1/2
    lambda <- Get_Lambda(K1)
    
    Pvalue <- CompQuadForm::davies(Q, lambda, h=rep(1, length(lambda)))$Qq
    Pvalue

    cat(paste("NO. Gene = ",ig,"\n"))
    
    res_each <- data.frame(geneid = rownames(object)[ig], 
                           Q=Q,
                           Pvalue=Pvalue)
    res_all <- rbind(res_all, res_each)
  }
  res_all$Adjust.Pvalue <- p.adjust(res_all$Pvalue, method = "BH")
  return(res_all)
}



#' @title Eigenvalue for Kernel matrix
#'
#' @description Calculating the eigenvalue for the kernel matrix K.
#'
#' @param K Kernel matrix.
#'
#' 
#'
#' @return This function returns eigenvalue for the Kernel matrix.
#' \item{lambda}{Eigenvalue for Kernel matrix.}
#' 
#'
#'
#' @export
Get_Lambda<-function(K){
  
  lambda<-NULL
  
  out.s<-eigen(K,symmetric=TRUE, only.values = TRUE)
  
  lambda1<-out.s$values
  IDX1<-which(lambda1 >= 0)
  
  # eigenvalue bigger than sum(eigenvalues)/1000
  IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
  
  if(length(IDX2) == 0){
    stop("No Eigenvalue is bigger than 0!!")
  }
  lambda<-lambda1[IDX2]
  
  return(lambda)
  
}




# 
# SPDE_test2 <- function(object, location, para){
#   
#   ED <- as.matrix(dist(location))
#   R2 <- as.matrix(dist(location) ** 2)
#   R <- as.matrix(dist(location))
#   
#   num_gene_test <- nrow(object)
#   
#   res_all <- NULL
#   for (ig in 1:num_gene_test){
#     
#     y.c <- as.numeric(object[ig, ])
#     
#     # if (para$Lik_Gau[ig] > para$Lik_Per[ig]){
#     K <- exp(-R2 / (2 * para$theta_Gau[ig] ** 2))
#     # } else{
#     # }
#     # K= exp(-ED^2/(2*(140^2)))
#     
#     obj<-SKAT::SKAT_Null_Model(y.c ~ 1, out_type="C")
#     res <- obj$res
#     s2 <- obj$s2
#     X <- obj$X1
#     
#     Q = t(res) %*% K %*% res / (2*s2)
#     W = K - X%*%solve( t(X)%*%X)%*%( t(X) %*% K)	# W = P0 K
#     W1 = W - (W %*% X) %*% solve(t(X)%*%X) %*% t(X)
#     K1 <- W1/2
#     lambda <- Get_Lambda(K1)
#     
#     Pvalue_Gau <- CompQuadForm::davies(Q, lambda, h=rep(1, length(lambda)))$Qq
#     Pvalue_Gau
#     
#     
#     K <- cos(2*pi*R / para$theta_Per[ig])
#     obj<-SKAT::SKAT_Null_Model(y.c ~ 1, out_type="C")
#     res <- obj$res
#     s2 <- obj$s2
#     X <- obj$X1
#     
#     Q = t(res) %*% K %*% res / (2*s2)
#     W = K - X%*%solve( t(X)%*%X)%*%( t(X) %*% K)	# W = P0 K
#     W1 = W - (W %*% X) %*% solve(t(X)%*%X) %*% t(X)
#     K1 <- W1/2
#     lambda <- Get_Lambda(K1)
#     
#     Pvalue_Per <- CompQuadForm::davies(Q, lambda, h=rep(1, length(lambda)))$Qq
#     Pvalue_Per
#     
#     if (para$Lik_Gau[ig] > para$Lik_Per[ig]){
#       Pvalue <- Pvalue_Gau
#       Index_Gau=1
#     } else{
#       Pvalue <- Pvalue_Per
#       Index_Gau=0
#     }
#     
#     cat(paste("NO. Gene = ",ig,"\n"))
#     
#     res_each <- data.frame(geneid = rownames(object)[ig], 
#                            Pvalue_Gau=Pvalue_Gau,
#                            Pvalue_Per=Pvalue_Per, 
#                            Pvalue=Pvalue, 
#                            Index_Gau=Index_Gau)
#     res_all <- rbind(res_all, res_each)
#   }
#   res_all$Adjust.Pvalue <- p.adjust(res_all$Pvalue, method = "BH")
#   return(res_all)
# }
