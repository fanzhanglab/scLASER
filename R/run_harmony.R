#' Run Harmony batch correction on stored NAM PCs
#'
#' Applies Harmony integration to the principal components stored in
#' \code{object@pcs} using batch information from \code{object@metadata}.
#' The harmonized embeddings are stored in \code{object@harmony} as a matrix
#' with the same row/column names as \code{object@pcs}.
#'
#' This is a convenience wrapper around \code{harmony::HarmonyMatrix()}.
#'
#' @param object A \code{\linkS4class{scLASER}} object. Must contain a non-\code{NULL}
#'   matrix in \code{object@pcs} and a data.frame in \code{object@metadata}.
#' @param batch_col Character scalar. Column name in \code{object@metadata} that
#'   encodes batch / sample / donor labels for integration (default \code{"sample_id"}).
#'   This column is coerced to a factor internally.
#' @param theta Numeric. Harmony diversity clustering penalty parameter
#'   (passed to \code{harmony::HarmonyMatrix()}).
#' @param do_pca Logical. If \code{TRUE}, Harmony will perform PCA internally.
#'   If \code{FALSE} (default), Harmony uses \code{object@pcs} as the input
#'   embedding.
#' @param epsilon.cluster Numeric. Convergence tolerance for the clustering step.
#'   Passed to \code{harmony::HarmonyMatrix()}.
#' @param epsilon.harmony Numeric. Convergence tolerance for Harmony correction.
#'   Passed to \code{harmony::HarmonyMatrix()}.
#' @param max.iter.cluster Integer. Maximum number of clustering iterations.
#' @param max.iter.harmony Integer. Maximum number of Harmony iterations.
#' @param plot_convergence Logical. If \code{TRUE}, Harmony may display convergence
#'   diagnostics/plots depending on the Harmony implementation.
#'
#' @return The input \code{scLASER} object with \code{object@harmony} populated as a
#'   numeric matrix of harmonized embeddings.
#'
#' @examples
#' \donttest{
#' ## Requires an scLASER object with PCs and metadata:
#' ## obj <- readRDS("path/to/your_sclaser_object.rds")
#' ## obj <- run_harmony(obj, batch_col = "sample_id", theta = 2, plot_convergence = FALSE)
#' ## dim(obj@harmony)
#' }
#'
#' @seealso \code{\link[harmony:HarmonyMatrix]{harmony::HarmonyMatrix}}
#'
#' @export
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
