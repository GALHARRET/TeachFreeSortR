#' Mean Adjusted Rand Index between a partition and a set of partitions
#'
#' @param partition A vector of group labels for the candidate (consensus)
#'   partition.
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns (individual partitions / group labels).
#'
#' @return The mean ARI between \code{partition} and each column of
#'   \code{data}.
#'
#' @export
mean_ari <- function(partition, data) {
  mean(vapply(seq_len(ncol(data)), function(j) {
    adjusted_rand_index(partition, data[[j]])
  }, numeric(1)))
}
