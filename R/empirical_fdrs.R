#' Empirical FDR estimates from observed and null scores
#'
#' Internal helper to compute empirical false discovery rates (FDR) across a set
#' of score thresholds by comparing tail counts in observed scores to tail counts
#' from a null distribution.
#'
#' @param z Numeric vector of observed scores (e.g., test statistics).
#' @param znull Numeric vector (or matrix) of null scores used to estimate the
#'   expected number of false positives at each threshold.
#' @param thresholds Numeric vector of bin/threshold cutpoints used by
#'   \code{tail_counts()}.
#'
#' @return Numeric vector of empirical FDR estimates (one per threshold interval).
#'
#' @keywords internal
#' @noRd
.empirical_fdrs <- function(z, znull, thresholds) {
  N <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:N, ])
  ranks <- t(tail_counts(thresholds, z)[1:N, ])

  fdp <- sweep(tails, 2, ranks, "/")
  fdr <- Matrix::colMeans(fdp)

  fdr
}
