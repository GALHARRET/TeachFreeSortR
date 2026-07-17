#' Compute the consensus partition from free-sorting data
#'
#' Derives the consensus partition of products from individual free-sorting
#' data, following Courcoux, Faye & Qannari (2014), replicating
#' \code{FreeSortR::ConsensusPartition()} (including its \code{optim}
#' local-search refinement with random permutation and singleton handling).
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns (group labels), as used in \code{\link{total_dissim}}.
#' @param ngroups Number of groups in the consensus partition. If \code{0}
#'   (default), the consensus is computed for every possible number of
#'   groups, and the best one (highest mean ARI) is returned.
#' @param type \code{"cutree"} (default), \code{"fusion"}, or
#'   \code{"medoid"}. See \code{FreeSortR::ConsensusPartition()} for details.
#' @param optim Logical. If \code{TRUE}, refines the partition with the
#'   exact local-search algorithm from \code{FreeSortR}: random-order
#'   single-product reassignment, repeated until the criterion no longer
#'   improves, followed (in case of unresolved singletons) by an attempt to
#'   swap two products' group memberships.
#' @param maxiter Maximum number of outer iterations for \code{optim}.
#'   Defaults to \code{100}.
#' @param plotDendrogram Logical. If \code{TRUE} and \code{type = "cutree"},
#'   plots the dendrogram used for initialization.
#' @param hclust_method Agglomeration method for \code{"cutree"}. Defaults
#'   to \code{"ward.D2"} (as in \code{FreeSortR}).
#' @param verbose Logical. Print progress information. Defaults to
#'   \code{FALSE}.
#'
#' @return A list with \code{consensus} (named group-assignment vector) and
#'   \code{crit} (mean ARI of the consensus with the individual partitions).
#'
#' @references
#' Courcoux, Ph., Faye, P., Qannari, E.M. (2014). Determination of the
#' consensus partition and cluster analysis of subjects in a free sorting
#' task experiment. \emph{Food Quality and Preference}, 32, 107-112.
#'
#' @importFrom stats hclust as.dist cutree
#' @export
consensus_partition <- function(data,
                                 ngroups = 0,
                                 type = c("cutree", "fusion", "medoid"),
                                 optim = FALSE,
                                 maxiter = 100,
                                 plotDendrogram = FALSE,
                                 hclust_method = "ward.D2",
                                 verbose = FALSE) {

  type <- match.arg(type)
  data <- as.data.frame(data)

  n_products <- nrow(data)
  npart <- ncol(data)
  products <- rownames(data)
  if (is.null(products)) products <- paste0("Product", seq_len(n_products))

  if (type == "fusion") {
    message("Fusion algorithm. May be time consuming.")
  }

  # ============================================================
  # type = "cutree" : build the hierarchical clustering
  # ============================================================
  if (type == "cutree") {
    d_mat <- total_dissim(data)
    hres <- stats::hclust(stats::as.dist(d_mat), method = hclust_method)
    if (plotDendrogram) {
      plot(hres, labels = products, hang = -1)
    }
  }

  # ============================================================
  # type = "fusion" : build the full agglomeration table
  # ============================================================
  if (type == "fusion") {
    fusion_res <- fusion_partition_table(data, ngroups = ngroups)
  }

  # ============================================================
  # type = "medoid" : closest partition to all others
  # ============================================================
  if (type == "medoid") {

    maxcrit <- 0
    pmax <- 0

    for (p in seq_len(npart)) {
      crit <- 0
      for (p1 in seq_len(npart)) {
        if (p != p1) {
          crit <- crit + adjusted_rand_index(data[[p]], data[[p1]])
        }
      }
      crit <- crit / (npart - 1)
      if (crit > maxcrit) {
        maxcrit <- crit
        pmax <- p
      }
    }

    consensus <- as.integer(factor(data[[pmax]]))
    names(consensus) <- products

    if (verbose) {
      message("Medoid partition")
      message(sprintf(
        "The consensus is the partition of subject '%s'.",
        colnames(data)[pmax]
      ))
    }
  }

  # ============================================================
  # searching for optimum (cutree / fusion only)
  # ============================================================
  if (type != "medoid") {

    if (ngroups != 0) {

      if (type == "cutree") {
        consensus <- stats::cutree(hres, k = ngroups)
        names(consensus) <- products
        maxcrit <- mean_ari(consensus, data)
      } else { # fusion
        consensus <- fusion_res$tabsubjopt[, 1]
        names(consensus) <- products
        maxcrit <- fusion_res$tabcritopt[1]
      }

      if (verbose) {
        message(sprintf("Criterion: %f", maxcrit))
        message("Consensus:")
        print(consensus)
      }

    } else {

      if (type == "cutree") {

        tabsubjopt <- matrix(0L, n_products, n_products - 2)
        tabcritopt <- matrix(0, 1, n_products - 2)
        colnames(tabsubjopt) <- 2:(n_products - 1)
        colnames(tabcritopt) <- 2:(n_products - 1)

        for (groups in 2:(n_products - 1)) {
          Cons <- stats::cutree(hres, k = groups)
          crit <- mean_ari(Cons, data)
          tabsubjopt[, groups - 1] <- Cons
          tabcritopt[1, groups - 1] <- crit
        }

      } else { # fusion
        tabsubjopt <- fusion_res$tabsubjopt
        tabcritopt <- matrix(fusion_res$tabcritopt, nrow = 1)
        colnames(tabsubjopt) <- names(fusion_res$tabcritopt)
        colnames(tabcritopt) <- names(fusion_res$tabcritopt)
      }

      groupmax <- colnames(tabcritopt)[which.max(tabcritopt)]
      maxcrit <- max(tabcritopt)
      consensus <- tabsubjopt[, colnames(tabsubjopt) == groupmax]
      names(consensus) <- products

      if (verbose) {
        message("Table of criterion as a function of the number of groups:")
        print(tabcritopt)
        message(sprintf(
          "\nOptimal consensus with %s groups and criterion %f:", groupmax, maxcrit
        ))
        print(consensus)
      }
    }
  }

  # ============================================================
  # optim = TRUE : exact local-search refinement (as in FreeSortR)
  # ============================================================
  if (type != "medoid" && optim) {

    converge <- FALSE
    iter <- 0
    obj <- 0
    nclass <- sum(unique(consensus) != 0)

    while (iter < maxiter && !converge) {

      change <- FALSE
      permut <- sample(seq_len(n_products))
      Singleton <- FALSE

      for (i in permut) {
        prov <- consensus

        if (sum(prov == prov[i]) > 1) {   # i is not a singleton

          maxcrit_i <- 0
          cmax <- 1

          for (cl in seq_len(nclass)) {
            prov[i] <- cl
            crit <- mean_ari(prov, data)
            if (crit > maxcrit_i) {
              maxcrit_i <- crit
              cmax <- cl
            }
          }

          if (consensus[i] != cmax) {
            consensus[i] <- cmax
            maxcrit <- maxcrit_i
            change <- TRUE
          }

        } else {
          Singleton <- TRUE
        }
      }

      converge <- abs(maxcrit - obj) < 1e-06

      if (converge && Singleton) {
        # attempt to exchange two products when blocked by a singleton
        for (i in permut) {
          if (sum(consensus == consensus[i]) == 1) {  # i is a singleton
            for (j in permut) {
              prov <- consensus
              gri <- prov[i]
              prov[i] <- prov[j]
              prov[j] <- gri

              crit <- mean_ari(prov, data)

              if (crit > maxcrit) {
                change <- TRUE
                maxcrit <- crit
                candi <- i
                candj <- j
              }
            }
          }
        }

        if (change) {
          grj <- consensus[candj]
          consensus[candj] <- consensus[candi]
          consensus[candi] <- grj
        }
      }

      obj <- maxcrit
      iter <- iter + 1
    }

    message(sprintf(
      "Final consensus with %d groups and criterion %f:",
      length(unique(consensus)), maxcrit
    ))
    print(consensus)
  }

  return(list(consensus = consensus, crit = maxcrit))
}
