#' Bootstrap MDS configurations with Procrustes alignment
#'
#' Performs bootstrap resampling over individuals to assess the stability of
#' an MDS configuration. Each bootstrap replicate's configuration is aligned
#' to the reference configuration via Procrustes analysis before being
#' collected, allowing confidence ellipses to be drawn per product.
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns (group labels), as used in \code{\link{total_dissim}}.
#' @param k Number of MDS dimensions. Defaults to \code{2}.
#' @param n_boot Number of bootstrap replicates. Defaults to \code{200}.
#' @param method MDS method, passed to \code{\link{compute_mds}}.
#'   Defaults to \code{"classical"}.
#'
#' @param sc Scaling or not the configuration. 
#' @return A data.frame with columns \code{product}, \code{Dim1}, ...,
#'   \code{DimK}, \code{replicate}, suitable for plotting with
#'   \code{ggplot2::stat_ellipse()} grouped by \code{product}.
#'
#' @importFrom stats cmdscale
#' @export
bootstrap_mds <- function(data, k = 2, n_boot = 200, method = "classical",sc=FALSE) {

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for Procrustes alignment. Please install it.")
  }

  data <- as.data.frame(data)
  products <- rownames(data)
  if (is.null(products)) products <- paste0("Product", seq_len(nrow(data)))
  n_individuals <- ncol(data)

  # --- Reference configuration (full data) ---
  d_ref <- total_dissim(data)
  ref_mds <- compute_mds(d_ref, k = k, method = method,sc=sc)
  ref_points <- ref_mds$points

  # --- Bootstrap replicates ---
  boot_list <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {

    idx <- sample(seq_len(n_individuals), n_individuals, replace = TRUE)
    data_boot <- data[, idx, drop = FALSE]

    d_boot <- total_dissim(data_boot)

    # Skip degenerate replicates (e.g. all individuals identical -> zero matrix)
    fit_boot <- tryCatch(
      compute_mds(d_boot, k = k, method = method,sc=sc),
      error = function(e) NULL
    )
    if (is.null(fit_boot)) next

    # Procrustes alignment onto the reference configuration
    proc <- vegan::procrustes(X = ref_points, Y = fit_boot$points, symmetric = FALSE)
    aligned_points <- proc$Yrot

    df_b <- as.data.frame(aligned_points)
    colnames(df_b) <- paste0("Dim", seq_len(k))
    df_b$product <- products
    df_b$replicate <- b

    boot_list[[b]] <- df_b
  }

  result <- do.call(rbind, boot_list)
  rownames(result) <- NULL

  return(result)
}
