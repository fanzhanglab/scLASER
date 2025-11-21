#' Title
#'
#' @noRd
#' @export
#'
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
