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
#' @importFrom spatialDE stabilize regress_out
SPADE_norm <- function(readcounts, info){
  stabilized <- spatialDE::stabilize(readcounts)
  info$total_counts <- colSums(readcounts)
  regdata <- spatialDE::regress_out(counts=stabilized, sample_info=info)
  return(regdata)
}

