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
compute_mds <- function(dissim_matrix, k = 2, method = c("classical", "nonmetric"),sc=FALSE) {

  method <- match.arg(method)

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

  d <- stats::as.dist(dissim_matrix)

  # --- Compute MDS ---
  if (method == "classical") {

    fit <- smacof::smacofSym(dissim_matrix, ndim = k, type = "interval", 
        verbose = FALSE)

    points <- fit$conf
    colnames(points) <- paste0("Dim", seq_len(k))
    Stress <- sqrt(sum((fit$dhat - fit$confdist)^2)/sum(fit$confdist^2))
    Config <- scale(points, center = TRUE, scale = sc)
    W <- Config %*% t(Config)
    bid <- svd(W)
    Config <- bid$u[, 1:k] %*% sqrt(diag(bid$d[1:k]))
    Percent <- bid$d[1:k]/sum(bid$d[1:k])  
    colnames(Config)<-paste("Dim",1:k,sep="")
    rownames(Config)<-rownames(D)  
    result <- list(
      points = Config,
      stress = Stress,
      explain_var=Percent,
      method = method
    )

  } else {

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for non-metric MDS. Please install it.")
    }

    if (any(d == 0)) {
      warning("Zero dissimilarities detected; adding a small jitter for isoMDS().")
      d[d == 0] <- d[d == 0] + 1e-6
    }

    fit <- smacof::smacofSym(dissim_matrix, ndim = k, type = "ordinal",ties="primary", 
        verbose = FALSE)

    points <- fit$conf
    #colnames(points) <- paste0("Dim", seq_len(k))

    Config <- scale(points, center = TRUE, scale = sc)
    W <- Config %*% t(Config)
    bid <- svd(W)
    Config <- bid$u[, 1:k] %*% sqrt(diag(bid$d[1:k]))
    Eig<-sqrt(bid$d[1:k])
    Stress <- sqrt(sum((fit$dhat - fit$confdist)^2)/sum(fit$confdist^2))
    Percent <- bid$d[1:k]/sum(bid$d[1:k])
    colnames(Config)<-paste("Dim",1:k,sep="")
    rownames(Config)<-rownames(D)  
    result <- list(
      points = Config,
      eig=Eig,
      stress = Stress,
      explain_var=Percent,
      method = method
    )
  }

  return(result)
}
