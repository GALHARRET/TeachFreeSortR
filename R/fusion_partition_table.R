#' Fusion algorithm for consensus partition (internal)
#'
#' Builds the full agglomeration table by iteratively merging, at each step,
#' the pair of groups that maximizes the mean Adjusted Rand Index with the
#' individual partitions. This follows the "fusion" algorithm from
#' Courcoux, Faye & Qannari (2014), as implemented in
#' \code{FreeSortR::ConsensusPartition(type = "fusion")}.
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns.
#' @param ngroups Target number of groups in the consensus. If \code{0}, the
#'   fusion is carried out all the way down to 2 groups, and results are kept
#'   for every intermediate number of groups.
#'
#' @return A list with \code{tabsubjopt} (matrix of partitions, one column
#'   per number of groups tested) and \code{tabcritopt} (named vector of the
#'   corresponding mean ARI criterion).
#'
#' @keywords internal
fusion_partition_table <- function(data, ngroups = 0) {

  n_products <- nrow(data)

  ngrclass <- if (ngroups != 0) n_products - ngroups + 1 else n_products - 1

  tabsubjopt <- matrix(0L, nrow = n_products, ncol = ngrclass - 1)
  tabcritopt <- numeric(ngrclass - 1)

  Cons <- seq_len(n_products)

  for (iter in seq_len(ngrclass - 1)) {

    nclasses <- n_products - iter + 1
    maxcrit <- -Inf
    maxgr1 <- NA_integer_
    maxgr2 <- NA_integer_

    for (gr1 in seq_len(nclasses - 1)) {
      for (gr2 in (gr1 + 1):nclasses) {

        # candidate fusion of gr1 and gr2
        prov <- Cons
        prov[prov == gr2] <- gr1
        prov[prov > gr2] <- prov[prov > gr2] - 1L

        crit <- mean_ari(prov, data)

        if (crit > maxcrit) {
          maxcrit <- crit
          maxgr1 <- gr1
          maxgr2 <- gr2
        }
      }
    }

    # apply the best fusion found
    Cons[Cons == maxgr2] <- maxgr1
    Cons[Cons > maxgr2] <- Cons[Cons > maxgr2] - 1L

    tabsubjopt[, ngrclass - iter] <- Cons
    tabcritopt[ngrclass - iter] <- maxcrit
  }

  group_counts <- (n_products - ngrclass + 1):(n_products - 1)
  colnames(tabsubjopt) <- group_counts
  names(tabcritopt) <- group_counts

  list(tabsubjopt = tabsubjopt, tabcritopt = tabcritopt)
}