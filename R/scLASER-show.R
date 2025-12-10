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
  cat("\t metadata: ",dim(metadata) )
  cat("\t pcs: ",dim(pcs ))
  cat("\t umap: ",dim(umap ))
  cat("\t harmony: ",dim(harmony ))
  cat("\t nam_pcs: ",dim(nam_pcs ))
  cat("\t NAM_matrix: ",dim(NAM_matrix ))
  cat("\t plsda_LV: ",dim(plsda_LV ))
  cat("\t pipeline_output: ",dim(pipeline_output ))
})

