#' Adjusted Rand Index distance 
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns (group labels), as used in \code{\link{total_dissim}}.
#' @param ind1 the first individual
#' @param ind2 the second individual
#' @return The ARI distance between the two individuals
#'
#' @export
ari_distance<-function(data,ind1,ind2){
  return(sqrt(1-adjusted_rand_index(data[,ind1],data[,ind2])))
}
