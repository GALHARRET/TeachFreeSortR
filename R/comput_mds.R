#' Compute multidimensional scaling (MDS) on a dissimilarity matrix
#'
#' Performs classical (metric) or non-metric multidimensional scaling on a
#' total dissimilarity matrix, such as the one produced by
#' \code{\link{total_dissim}}, and returns product coordinates in a
#' low-dimensional space.
#'
#' @param dissim_matrix A square symmetric dissimilarity matrix (numeric),
#'   typically the output of \code{\link{total_dissim}}.
#' @param k Number of dimensions to retain. Defaults to \code{2}.
#' @param method MDS method: \code{"classical"} (default) uses
#'   \code{\link[stats]{cmdscale}} (metric MDS based on eigen-decomposition).
#'   \code{"nonmetric"} uses \code{\link[MASS]{isoMDS}} (Kruskal's non-metric
#'   MDS), which requires the \pkg{MASS} package and a strictly positive
#'   dissimilarity matrix (zeros on the diagonal are handled automatically,
#'   but off-diagonal zeros should be avoided or jittered).
#' @param sc Scaling or not the configuration 
#' @return A list with the following elements:
#' \describe{
#'   \item{points}{A matrix of product coordinates, with \code{k} columns.}
#'   \item{eig}{Eigenvalues (classical MDS only, \code{NULL} otherwise).}
#'   \item{goodness_of_fit}{Goodness-of-fit statistic (classical MDS) or
#'     final stress value (non-metric MDS).}
#'   \item{method}{The method used.}
#' }
#'
#' @examples
#' df <- data.frame(
#'   ind1 = c("A", "A", "B", "B"),
#'   ind2 = c(1, 2, 2, 1),
#'   ind3 = c("X", "X", "X", "Y")
#' )
#' rownames(df) <- paste0("Product", 1:4)
#'
#' d <- total_dissim(df)
#' res <- compute_mds(d)
#' res$points
#'
#' @importFrom stats cmdscale as.dist
#' @export
compute_mds <- function(dissim_matrix, k = 2, method = c("classical", "nonmetric"), sc = FALSE) {

  method <- match.arg(method)

  if (!requireNamespace("smacof", quietly = TRUE)) {
    stop("Package 'smacof' is required. Please install it.")
  }

  # --- Input checks ---
  if (!is.matrix(dissim_matrix)) {
    dissim_matrix <- as.matrix(dissim_matrix)
  }

  if (nrow(dissim_matrix) != ncol(dissim_matrix)) {
    stop("`dissim_matrix` must be a square matrix.")
  }

  if (!isSymmetric(unname(dissim_matrix))) {
    warning("`dissim_matrix` is not exactly symmetric; symmetrizing it.")
    dissim_matrix <- (dissim_matrix + t(dissim_matrix)) / 2
  }

  if (nrow(dissim_matrix) <= k) {
    stop("`k` must be smaller than the number of products.")
  }

  labels <- rownames(dissim_matrix)

  # --- Helper: rotate configuration to orthogonal, variance-ordered axes ---
  rotate_config <- function(points, k, sc) {
    Config <- scale(points, center = TRUE, scale = sc)
    W <- Config %*% t(Config)
    bid <- svd(W)

    if (k > 1) {
      Config <- bid$u[, 1:k] %*% sqrt(diag(bid$d[1:k]))
    } else {
      Config <- as.matrix(bid$u[, 1] * sqrt(bid$d[1]))
    }

    Percent <- bid$d[1:k] / sum(bid$d[1:k])
    colnames(Config) <- paste0("Dim", seq_len(k))
    rownames(Config) <- labels

    list(config = Config, percent = Percent, eig = sqrt(bid$d[1:k]))
  }

  # --- Compute MDS ---
  if (method == "classical") {

    fit <- smacof::smacofSym(dissim_matrix, ndim = k, type = "interval", verbose = FALSE)
    Stress <- sqrt(sum((fit$dhat - fit$confdist)^2) / sum(fit$confdist^2))

    rot <- rotate_config(fit$conf, k, sc)

    result <- list(
      points      = rot$config,
      eig         = rot$eig,
      stress      = Stress,
      explain_var = rot$percent,
      method      = method
    )

  } else {

    if (any(dissim_matrix[lower.tri(dissim_matrix)] == 0)) {
      warning("Zero dissimilarities detected; adding a small jitter before smacofSym().")
      dissim_matrix[lower.tri(dissim_matrix)][dissim_matrix[lower.tri(dissim_matrix)] == 0] <- 1e-6
      dissim_matrix <- (dissim_matrix + t(dissim_matrix)) / 2  # re-symmetrize after jitter
    }

    fit <- smacof::smacofSym(dissim_matrix, ndim = k, type = "ordinal", ties = "primary", verbose = FALSE)
    Stress <- sqrt(sum((fit$dhat - fit$confdist)^2) / sum(fit$confdist^2))

    rot <- rotate_config(fit$conf, k, sc)

    result <- list(
      points      = rot$config,
      eig         = rot$eig,
      stress      = Stress,
      explain_var = rot$percent,
      method      = method
    )
  }

  return(result)
}