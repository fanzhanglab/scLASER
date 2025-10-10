setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL",     c("data.frame", "NULL"))

setClass(
  "scLASER",
  slots = list(
    metadata = "data.frame",   # sample/cell metadata
    pcs      = "matrixOrNULL", # principal components
    masc    = "dfOrNULL", # adjacency or kNN graph (optional)
    umap     = "dfOrNULL"      # UMAP coordinates (optional)
  ),
  prototype = list(
    metadata = data.frame(),
    pcs      = NULL,
    masc    = NULL,
    umap     = NULL
  )
)

#wrapper
scLASER <- function(metadata = data.frame(),
                    pcs = NULL,
                    masc = NULL,
                    umap = NULL) {
  new("scLASER", metadata = metadata, pcs = pcs, masc = masc, umap = umap)
}
