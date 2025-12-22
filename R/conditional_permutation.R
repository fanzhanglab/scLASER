#' Conditional permutation within blocks
#'
#' Internal helper used to generate block-wise resampled permutations of a
#' response vector. Values in \code{Y} are resampled with replacement within
#' levels of the blocking factor \code{B}.
#'
#' @param B A vector defining permutation blocks (e.g., cluster or subject IDs).
#' @param Y A numeric or factor vector to be permuted within blocks.
#' @param num Integer. Number of permutation replicates to generate.
#'
#' @return A matrix where each column corresponds to a block-wise permuted
#'   version of \code{Y}.
#'
#' @keywords internal
#' @noRd
.conditional_permutation <- function(B, Y, num) {
  purrr::map(seq_len(num), function(i) {
    split(seq_len(length(Y)), B) %>%
      purrr::map(function(idx) {
        data.frame(idx, val = sample(Y[idx], replace = TRUE))
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(idx) %>%
      with(val)
  }) %>%
    purrr::reduce(Matrix::cbind2)
}

