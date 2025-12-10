#' Show method for scLASER objects
#'
#' Custom display when a scLASER object is printed.
#'
#' @param object A scLASER object.
#' @export
setMethod("show", "scLASER", function(object) {
  cat("scLASER object\n")
  cat("---------------\n")

})

