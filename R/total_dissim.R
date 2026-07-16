#' Compute the total dissimilarity matrix across individuals
#'
#' Aggregates individual dissimilarity matrices (as produced by
#' \code{\link{dissim_partition}}) into a single total dissimilarity matrix
#' across all individuals, either as a sum or an average.
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns. Each cell contains the group label (numeric, character, or
#'   factor) assigned by the individual to that product (typical output of a
#'   free-sorting task). Row names (if present) are used as product
#'   identifiers.
#' @param aggregate Aggregation method: \code{"sum"} (default) returns, for
#'   each pair of products, the number of individuals who sorted them into
#'   different groups. \code{"mean"} returns the same quantity divided by the
#'   number of individuals, i.e. the proportion of individuals who separated
#'   the two products (a value between 0 and 1).
#'
#' @return A square symmetric matrix of total dissimilarity between products,
#'   with 0 on the diagonal.
#'
#' @examples
#' df <- data.frame(
#'   ind1 = c("A", "A", "B", "B"),
#'   ind2 = c(1, 2, 2, 1),
#'   ind3 = c("X", "X", "X", "Y")
#' )
#' rownames(df) <- paste0("Product", 1:4)
#'
#' total_dissim(df)
#' total_dissim(df, aggregate = "mean")
#'
#' @export
total_dissim <- function(data, aggregate = c("sum", "mean")) {

  aggregate <- match.arg(aggregate)

  # --- Input checks ---
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a matrix or a data.frame.")
  }

  data <- as.data.frame(data)

  if (nrow(data) < 2) {
    stop("`data` must contain at least two products (rows).")
  }

  if (anyNA(data)) {
    warning("`data` contains missing values; affected pairs will be NA in the output.")
  }

  products <- rownames(data)
  if (is.null(products)) {
    products <- paste0("Product", seq_len(nrow(data)))
  }

  n_products <- nrow(data)
  n_individuals <- ncol(data)

  # --- Sum dissimilarity matrices across individuals ---
  total_mat <- matrix(0, nrow = n_products, ncol = n_products)

  for (j in seq_len(n_individuals)) {
    groups <- data[[j]]
    mat_j <- outer(groups, groups, FUN = function(a, b) as.integer(a != b))
    total_mat <- total_mat + mat_j
  }

  if (aggregate == "mean") {
    total_mat <- total_mat / n_individuals
  }

  rownames(total_mat) <- products
  colnames(total_mat) <- products

  return(total_mat)
}
