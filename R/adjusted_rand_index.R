#' Adjusted Rand Index between two partitions
#'
#' Computes the Adjusted Rand Index (ARI) between two partitions of the same
#' set of objects, as defined by Hubert & Arabie (1985).
#'
#' @param part1,part2 Two vectors of group labels (numeric, character, or
#'   factor), of the same length, giving the partition of each object under
#'   two different classifications.
#'
#' @return A single numeric value (ARI), equal to 1 for identical partitions
#'   and close to 0 for partitions no more similar than chance.
#'
#' @keywords internal

adjusted_rand_index <- function(part1, part2) {

  tab <- table(part1, part2)

  n <- sum(tab)
  sum_comb_tab <- sum(choose(tab, 2))
  sum_comb_rows <- sum(choose(rowSums(tab), 2))
  sum_comb_cols <- sum(choose(colSums(tab), 2))
  comb_n <- choose(n, 2)

  expected_index <- (sum_comb_rows * sum_comb_cols) / comb_n
  max_index <- (sum_comb_rows + sum_comb_cols) / 2

  if (max_index == expected_index) return(1)

  (sum_comb_tab - expected_index) / (max_index - expected_index)
}