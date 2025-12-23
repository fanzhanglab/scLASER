#' Tail counts of null statistics beyond thresholds
#'
#' Internal helper used for empirical FDR estimation. For each column in a null
#' statistic matrix, counts (in a cumulative/tail sense) how many null values
#' exceed a set of threshold cutpoints based on squared magnitude.
#'
#' @param thresholds Numeric vector of threshold cutpoints. These define bins on
#'   \eqn{z^2} via \code{cut()} using breaks \code{c(0, thresholds^2)}.
#' @param znull Numeric matrix (or data.frame coercible to matrix) of null scores.
#'   Each column is treated as one null replicate.
#'
#' @return Numeric matrix of tail counts with one column per null replicate.
#'
#' @keywords internal
#' @noRd
.tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
  })
}
