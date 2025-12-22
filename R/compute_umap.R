#' Compute a UMAP embedding from stored PCs
#'
#' Computes a two-dimensional UMAP embedding using the principal components
#' stored in \code{object@pcs} and saves the result in \code{object@umap}
#' as a data.frame with columns \code{UMAP1} and \code{UMAP2}.
#'
#' This is a convenience wrapper around \code{uwot::umap()}.
#'
#' @param object A \code{\linkS4class{scLASER}} object. Must contain a
#'   non-\code{NULL} numeric matrix in \code{object@pcs}.
#' @param n_neighbors Integer. Number of nearest neighbors used in UMAP.
#'   Larger values preserve more global structure. Passed to
#'   \code{uwot::umap()}.
#' @param metric Character. Distance metric for UMAP. Passed to
#'   \code{uwot::umap()} (e.g., \code{"cosine"}, \code{"euclidean"}).
#' @param min_dist Numeric. Minimum distance between points in the embedding.
#'   Smaller values produce tighter clusters. Passed to \code{uwot::umap()}.
#' @param ... Additional arguments passed to \code{uwot::umap()} (e.g.,
#'   \code{set_op_mix_ratio}, \code{spread}, \code{n_threads}, \code{seed}).
#'
#' @return The input \code{scLASER} object with \code{object@umap} populated as a
#'   data.frame containing \code{UMAP1} and \code{UMAP2}.
#'
#' @examples
#' \donttest{
#' ## Example with a precomputed PCs matrix stored in an scLASER object:
#' ## obj <- readRDS("path/to/your_sclaser_object.rds")
#' ## obj <- compute_umap(obj, n_neighbors = 30, metric = "cosine", min_dist = 0.01)
#' ## head(obj@umap)
#' }
#'
#' @export
compute_umap <- function(object,
                         n_neighbors = 30,
                         metric = "cosine",
                         min_dist = 0.01,
                         ...) {
  if (!inherits(object, "scLASER")) stop("object must be an 'scLASER'")
  if (is.null(object@pcs)) stop("object@pcs is NULL; provide a PCs matrix.")
  pcs_matrix <- object@pcs
  if (!is.matrix(pcs_matrix)) stop("object@pcs must be a matrix.")

  umap_res <- uwot::umap(
    pcs_matrix,
    n_neighbors = n_neighbors,
    metric = metric,
    min_dist = min_dist,
    ...
  )
  umap_res <- as.data.frame(umap_res)
  names(umap_res) <- c("UMAP1", "UMAP2")

  object@umap <- umap_res
  object
}
