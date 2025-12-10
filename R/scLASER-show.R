#' Show method for scLASER objects
#'
#' Custom display when a scLASER object is printed.
#'
#' @param object A scLASER object.
#' @export
setMethod("show", "scLASER", function(object) {
  cat("scLASER object\n")
  cat("---------------\n")
  cat("Slots:","\n" )
  cat("\t metadata: ",dim(object@metadata),"\n" )
  cat("\t pcs: ",dim(object@pcs ),"\n" )
  cat("\t umap: ",dim(object@umap ),"\n" )
  cat("\t harmony: ",dim(object@harmony ),"\n" )
  cat("\t nam_pcs: ",dim(object@nam_pcs ),"\n" )
  cat("\t NAM_matrix: ",dim(object@NAM_matrix ),"\n" )
  cat("\t plsda_LV: ",dim(object@plsda_LV ),"\n" )
  cat("\t pipeline_output: ",dim(object@pipeline_output ),"\n" )
})

