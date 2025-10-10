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
