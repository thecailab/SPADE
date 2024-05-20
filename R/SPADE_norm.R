#' @title SPADE normalization
#'
#' @description Normlization of read counts data into continuous data.
#'
#' @param readcounts Read counts data for identifying spatially expressed genes. Each row is a gene and each column is a spot.
#' @param info Spatial coordinates for all spots.
#'
#'
#' @return This function returns the normalized continuous data.
#' \item{regdata}{The normlized continuous data}
#' 
#'
#' @export
SPADE_norm <- function(readcounts, info){
  vst_counts <- stabilize(readcounts)
  info1$log_total_counts <- log(info1$total_counts)
  regdata <- regress_out(info1, vst_counts, "~ log_total_counts")
  return(regdata)
}



#' @title SPADE stabilize variance
#'
#' @description Use Anscombe's approximation to variance stabilize Negative Binomial data.
#'
#' @param expression_matrix Read counts data for identifying spatially expressed genes. Each row is a gene and each column is a spot.
#'
#' @return This function returns the stabilized data.
#' \item{stabilized_matrix}{The stabilized data}
#' 
#'
#' @importFrom stats lm nls
stabilize <- function(expression_matrix) {
  # Assumes columns are samples, and rows are genes
  
  # Calculate the mean and variance of each row (gene)
  gene_means <- rowMeans(expression_matrix)
  gene_vars <- apply(expression_matrix, 1, var)
  
  # Define the function to fit
  model <- function(mu, phi) {
    mu + phi * mu^2
  }
  
  # Perform the curve fitting to estimate phi_hat
  fit <- nls(gene_vars ~ model(gene_means, phi), start = list(phi = 0.1))
  phi_hat <- coef(fit)["phi"]
  
  # Apply the Anscombe transformation
  stabilized_matrix <- log(expression_matrix + 1 / (2 * phi_hat))
  
  return(stabilized_matrix)
}



#' @title SPADE residualize data
#'
#' @description Implementation of limma's removeBatchEffect function to residualize data.
#'
#' @param sample_info Spatial coordinates for all spots.
#' @param expression_matrix Stabilized data from above step. Each row is a gene and each column is a spot.
#' @param covariate_formula Specify covariates to regress out.
#' @param design_formula Design formula.
#'
#' @return This function returns the residualize data.
#' \item{regressed}{The residualized data}
#' 
#' @importFrom stats lmFit
regress_out <- function(sample_info, expression_matrix, covariate_formula, design_formula = "~ 1") {
  # Ensure intercept is not part of covariates
  covariate_formula <- paste0(covariate_formula, " - 1")
  
  # Create design and covariate matrices
  covariate_matrix <- model.matrix(as.formula(covariate_formula), data = sample_info)
  design_matrix <- model.matrix(as.formula(design_formula), data = sample_info)
  
  # Combine design and covariate matrices
  design_batch <- cbind(design_matrix, covariate_matrix)
  
  # Perform linear regression to find coefficients
  fit <- lmFit(expression_matrix, design_batch)
  coefficients <- fit$coefficients
  
  # Extract beta (coefficients of covariates)
  beta <- coefficients[, (ncol(design_matrix) + 1):ncol(design_batch)]
  
  # Regress out the covariates
  regressed <- expression_matrix - t(covariate_matrix %*% t(beta))
  
  return(regressed)
}
