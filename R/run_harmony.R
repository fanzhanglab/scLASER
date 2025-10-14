#' Title
#'
#' @param object 
#' @param batch_col 
#' @param theta 
#' @param do_pca 
#' @param epsilon.cluster 
#' @param epsilon.harmony 
#' @param max.iter.cluster 
#' @param max.iter.harmony 
#' @param plot_convergence 
#'
#' @return
#' @export
#'
#' @examples
run_harmony <- function(object,
                        batch_col = "sample_id",
                        theta = 2,
                        do_pca = FALSE,
                        epsilon.cluster = -Inf,
                        epsilon.harmony = -Inf,
                        max.iter.cluster = 30,
                        max.iter.harmony = 10,
                        plot_convergence = TRUE) {
  stopifnot(is(object, "scLASER"))
  if (is.null(object@pcs)) stop("@pcs is NULL in scLASER object.")
  if (!is.data.frame(object@metadata)) stop("@metadata must be a data.frame.")
  if (!batch_col %in% names(object@metadata))
    stop(sprintf("Column '%s' not found in @metadata.", batch_col))
  if (nrow(object@pcs) != nrow(object@metadata))
    stop("Row counts of @pcs and @metadata must match for Harmony.")

  meta <- object@metadata
  meta[[batch_col]] <- as.factor(meta[[batch_col]])

  harmonized <- harmony::HarmonyMatrix(
    data_mat         = object@pcs,
    meta_data        = meta,
    do_pca           = do_pca,
    vars_use         = batch_col,
    theta            = theta,
    epsilon.cluster  = epsilon.cluster,
    epsilon.harmony  = epsilon.harmony,
    max.iter.cluster = max.iter.cluster,
    max.iter.harmony = max.iter.harmony,
    plot_convergence = plot_convergence
  )

  harmonized <- as.matrix(harmonized)
  rownames(harmonized) <- rownames(object@pcs)
  colnames(harmonized) <- colnames(object@pcs)

  object@harmony <- harmonized
  object
}
