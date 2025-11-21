#' Title
#'
#' @noRd
#' @export
#'
tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
  })
}
