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
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom stats lm resid
SPADE_norm <- function(readcounts, info){
  NGS <- round(apply(readcounts, 2, as.numeric))
  vst_counts <- DESeq2::varianceStabilizingTransformation(NGS)
  info1$total_counts <- colSums(vst_counts)
	regdata <- t(apply(vst_counts, 1, function(x){resid(lm(x ~ log(info1$total_counts)))} ))
  rownames(regdata) <- rownames(readcounts)
  return(regdata)
}

