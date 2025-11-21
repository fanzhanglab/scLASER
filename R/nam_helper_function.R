#' Title (optional, can be brief)
#' @noRd
#' @import internal
conditional_permutation <- function(B, Y, num) {
  purrr::map(seq_len(num), function(i) {
    split(seq_len(length(Y)), B) %>%
      purrr::map(function(idx) {
        data.frame(idx, val=sample(Y[idx], replace = TRUE))
      }) %>% dplyr::bind_rows() %>%
      dplyr::arrange(idx) %>%
      with(val)
  }) %>%
    purrr::reduce(Matrix::cbind2)
}




#' Title (optional, can be brief)
#' @noRd
#' @import internal
empirical_fdrs <- function(z, znull, thresholds) {
  N <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:N, ])
  ranks <- t(tail_counts(thresholds, z)[1:N, ])


  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)

  return(fdr)
}

#' Title (optional, can be brief)
#' @noRd
#' @import
tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
  })
}
