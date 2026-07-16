#' Plot MDS configurations with Procrustes-aligned bootstrap confidence ellipses
#'
#' Performs bootstrap resampling over individuals to assess the stability of
#' an MDS configuration. Each bootstrap replicate's configuration is aligned
#' to the reference configuration via Procrustes analysis before being
#' collected, allowing confidence ellipses to be drawn per product. If
#' \code{n_boot = 0}, only the reference configuration is plotted, without
#' ellipses.
#'
#' @param data A matrix or data.frame with products in rows and individuals
#'   in columns (group labels), as used in \code{\link{total_dissim}}.
#' @param k Number of MDS dimensions. Defaults to \code{2}.
#' @param n_boot Number of bootstrap replicates. Defaults to \code{500}. If
#'   \code{n_boot = 0}, the reference configuration is plotted without
#'   confidence ellipses.
#' @param method MDS method, passed to \code{\link{compute_mds}}.
#'   Defaults to \code{"classical"}.
#' @param sc Logical. Whether to scale the configuration.
#' @param dim A vector of length 2 giving the indices of the dimensions to
#'   plot. Defaults to \code{c(1, 2)}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom stats cmdscale
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_bw labs
#'   ggtitle theme coord_fixed
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#' @export
mds_plot <- function(data, k = 2, n_boot = 500, sc = FALSE, method = "classical", dim = 1:2) {

  # --- Input checks ---
  if (length(dim) != 2) {
    stop("`dim` must be a vector of length 2.")
  }
  if (any(dim > k)) {
    stop(sprintf("`dim` requests dimension(s) up to %d, but k = %d.", max(dim), k))
  }
  if (n_boot < 0) {
    stop("`n_boot` must be >= 0.")
  }

  Diss <- total_dissim(data)
  config <- compute_mds(dissim_matrix = Diss, k = k, method = method, sc = sc)

  dim1_name <- paste0("Dim", dim[1])
  dim2_name <- paste0("Dim", dim[2])

  # --- Axis labels (guarded: explain_var may be absent for nonmetric MDS) ---
  if (!is.null(config$explain_var)) {
    x_lab <- paste0(dim1_name, " (", round(100 * config$explain_var[dim[1]], 1), "%)")
    y_lab <- paste0(dim2_name, " (", round(100 * config$explain_var[dim[2]], 1), "%)")
  } else {
    x_lab <- dim1_name
    y_lab <- dim2_name
  }

  # --- Title (guarded: stress may be absent) ---
  plot_title <- if (!is.null(config$stress)) {
    paste("Stress =", round(config$stress, 3))
  } else {
    NULL
  }

  if (n_boot > 0) {

    boot_df <- bootstrap_mds(data, k = k, n_boot = n_boot, method = method, sc = sc)

    centres_boot <- boot_df |>
      dplyr::group_by(.data$product) |>
      dplyr::summarise(
        !!dim1_name := mean(.data[[dim1_name]]),
        !!dim2_name := mean(.data[[dim2_name]]),
        .groups = "drop"
      )

    p <- ggplot2::ggplot(boot_df, ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                                color = .data$product, fill = .data$product)) +
      ggplot2::stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.90) +
      ggplot2::geom_point(data = centres_boot,
                           ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]]),
                           color = "black", size = 2, inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = centres_boot,
                                ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                             label = .data$product),
                                color = "black")

  } else {

    df_mds <- as.data.frame(config$points)
    df_mds$product <- rownames(config$points)

    p <- ggplot2::ggplot(df_mds, ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                               label = .data$product)) +
      ggplot2::geom_point() +
      ggrepel::geom_text_repel()
  }

  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme(legend.position = "none")

  if (!is.null(plot_title)) {
    p <- p + ggplot2::ggtitle(plot_title)
  }

  return(p)
}
