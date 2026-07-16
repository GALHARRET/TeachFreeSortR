#' Plot MDS configurations colored by consensus partition
#'
#' Same as \code{\link{mds_plot}}, but points and confidence ellipses are
#' colored according to the consensus partition of products (as returned by
#' \code{\link{consensus_partition}}), instead of one color per product.
#' Product labels are still displayed individually via \code{ggrepel}.
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
#' @param consensus Optional. A pre-computed result from
#'   \code{\link{consensus_partition}} (a list with a \code{consensus}
#'   element). If \code{NULL} (default), it is computed internally using
#'   \code{ngroups}, \code{partition_type}, and \code{partition_optim}.
#' @param ngroups Number of groups for the consensus partition, passed to
#'   \code{\link{consensus_partition}} if \code{consensus} is not supplied.
#'   Defaults to \code{0} (automatic selection).
#' @param partition_type Type of consensus algorithm, passed to
#'   \code{\link{consensus_partition}}. Defaults to \code{"cutree"}.
#' @param partition_optim Logical, passed to \code{\link{consensus_partition}}
#'   as \code{optim}. Defaults to \code{FALSE}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse theme_bw labs
#'   ggtitle theme coord_fixed scale_color_discrete scale_fill_discrete
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by summarise left_join
#' @importFrom rlang .data
#' @export
mds_plot_consensus <- function(data,
                                k = 2,
                                n_boot = 500,
                                method = "classical",
                                sc = FALSE,
                                dim = 1:2,
                                consensus = NULL,
                                ngroups = 0,
                                partition_type = "cutree",
                                partition_optim = FALSE) {

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

  # --- Consensus partition ---
  if (is.null(consensus)) {
    consensus <- consensus_partition(
      data,
      ngroups = ngroups,
      type    = partition_type,
      optim   = partition_optim
    )
  }

  if (!"consensus" %in% names(consensus)) {
    stop("`consensus` must be a list containing a `consensus` element (as returned by consensus_partition()).")
  }

  group_df <- data.frame(
    product = names(consensus$consensus),
    group   = factor(consensus$consensus)
  )

  # --- MDS configuration ---
  Diss <- total_dissim(data)
  config <- compute_mds(dissim_matrix = Diss, k = k, method = method, sc = sc)

  dim1_name <- paste0("Dim", dim[1])
  dim2_name <- paste0("Dim", dim[2])

  # --- Axis labels ---
  if (!is.null(config$explain_var)) {
    x_lab <- paste0(dim1_name, " (", round(100 * config$explain_var[dim[1]], 1), "%)")
    y_lab <- paste0(dim2_name, " (", round(100 * config$explain_var[dim[2]], 1), "%)")
  } else {
    x_lab <- dim1_name
    y_lab <- dim2_name
  }

  # --- Title ---
  plot_title <- if (!is.null(config$stress)) {
    paste("Stress =", round(config$stress, 3))
  } else {
    NULL
  }
  ngroups_used <- length(unique(consensus$consensus))
  plot_subtitle <- paste("Consensus partition:", ngroups_used, "groups")

  if (n_boot > 0) {

    boot_df <- bootstrap_mds(data, k = k, n_boot = n_boot, method = method, sc = sc)
    boot_df <- dplyr::left_join(boot_df, group_df, by = "product")

    centres_boot <- boot_df |>
      dplyr::group_by(.data$product, .data$group) |>
      dplyr::summarise(
        !!dim1_name := mean(.data[[dim1_name]]),
        !!dim2_name := mean(.data[[dim2_name]]),
        .groups = "drop"
      )

    p <- ggplot2::ggplot(boot_df, ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                                color = .data$group, fill = .data$group)) +
      ggplot2::stat_ellipse(ggplot2::aes(group = .data$product),
                             geom = "polygon", alpha = 0.15, level = 0.90) +
      ggplot2::geom_point(data = centres_boot,
                           ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                        color = .data$group),
                           size = 2, inherit.aes = FALSE) +
      ggrepel::geom_text_repel(data = centres_boot,
                                ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                             label = .data$product, color = .data$group),
                                show.legend = FALSE, inherit.aes = FALSE)

  } else {

    df_mds <- as.data.frame(config$points)
    df_mds$product <- rownames(config$points)
    df_mds <- dplyr::left_join(df_mds, group_df, by = "product")

    p <- ggplot2::ggplot(df_mds, ggplot2::aes(x = .data[[dim1_name]], y = .data[[dim2_name]],
                                               label = .data$product, color = .data$group)) +
      ggplot2::geom_point(size = 2) +
      ggrepel::geom_text_repel(show.legend = FALSE)
  }

  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::labs(x = x_lab, y = y_lab, color = "Consensus group", fill = "Consensus group") +
    ggplot2::theme(legend.position = "right")

  if (!is.null(plot_title)) {
    p <- p + ggplot2::ggtitle(plot_title, subtitle = plot_subtitle)
  }

  return(p)
}