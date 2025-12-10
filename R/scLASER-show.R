#' Show method for scLASER objects
#'
#' Custom display when a scLASER object is printed.
#'
#' @param object A scLASER object.
#' @export
setMethod("show", "scLASER", function(object) {
  cat("scLASER object\n")
  cat("---------------\n")
  cat("Slots:")
  cat("\t metadata: ",dim(object@metadata) )
  cat("\t pcs: ",dim(object@pcs ))
  cat("\t umap: ",dim(object@umap ))
  cat("\t harmony: ",dim(object@harmony ))
  cat("\t nam_pcs: ",dim(object@nam_pcs ))
  cat("\t NAM_matrix: ",dim(object@NAM_matrix ))
  cat("\t plsda_LV: ",dim(object@plsda_LV ))
  cat("\t pipeline_output: ",dim(object@pipeline_output ))
})

