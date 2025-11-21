#' Title
#'
#' @noRd
#' @export
#'
empirical_fdrs <- function(z, znull, thresholds) {
  N <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:N, ])
  ranks <- t(tail_counts(thresholds, z)[1:N, ])
  
  
  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)
  
  return(fdr)
}