#' Compute individual dissimilarity (partition) matrices
#'
#' For each individual (column), computes a binary dissimilarity matrix
#' indicating whether pairs of products were sorted into different groups
#' by that individual.
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns. Each cell contains the group label (numeric, character, or
#'   factor) assigned by the individual to that product (typical output of a
#'   free-sorting task). Row names (if present) are used as product
#'   identifiers.
#'
#' @return A named list (named after the column names of \code{data}),
#'   containing for each individual a square symmetric binary matrix with
#'   0 if the two products were sorted into the same group, and 1 otherwise.
#'   The diagonal is always 0 (a product is never dissimilar to itself).
#'
#' @examples
#' df <- data.frame(
#'   ind1 = c("A", "A", "B", "B"),
#'   ind2 = c(1, 2, 2, 1),
#'   ind3 = c("X", "X", "X", "Y")
#' )
#' rownames(df) <- paste0("Product", 1:4)
#'
#' res <- dissim_partition(df)
#' res$ind1
#'
#' @export
dissim_partition <- function(data) {

  # --- Input checks ---
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("`data` must be a matrix or a data.frame.")
  }

  data <- as.data.frame(data)

  if (nrow(data) < 2) {
    stop("`data` must contain at least two products (rows).")
  }

  if (anyNA(data)) {
    warning("`data` contains missing values; corresponding cells will be NA in the output.")
  }

  products <- rownames(data)
  if (is.null(products)) {
    products <- paste0("Product", seq_len(nrow(data)))
  }

  individuals <- colnames(data)
  if (is.null(individuals)) {
    individuals <- paste0("Individual", seq_len(ncol(data)))
  }

  # --- Compute dissimilarity matrix for each individual ---
  dissim_list <- lapply(seq_len(ncol(data)), function(j) {
    groups <- data[[j]]

    mat <- outer(groups, groups, FUN = function(a, b) as.integer(a != b))

    rownames(mat) <- products
    colnames(mat) <- products

    mat
  })

  names(dissim_list) <- individuals

  return(dissim_list)
}