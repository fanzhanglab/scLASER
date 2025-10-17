setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL",     c("data.frame", "NULL"))

#' Title
#' @export
setClass(
  "scLASER",
  slots = list(
    metadata        = "data.frame",
    pcs             = "matrixOrNULL",
    masc            = "dfOrNULL",
    umap            = "dfOrNULL",
    harmony         = "matrixOrNULL",
    nam_pcs         = "matrixOrNULL",
    pipeline_output = "dfOrNULL"
  ),
  prototype = list(
    metadata        = data.frame(),
    pcs             = NULL,
    masc            = NULL,
    umap            = NULL,
    harmony         = NULL,
    nam_pcs         = NULL,
    pipeline_output = NULL
  )
)


scLASER <- function(metadata = data.frame(),
                    pcs = NULL,
                    masc = NULL,
                    umap = NULL,
                    harmony = NULL,
                    nam_pcs = NULL,
                    pipeline_output = NULL) {
  new("scLASER",
      metadata        = metadata,
      pcs             = pcs,
      masc            = masc,
      umap            = umap,
      harmony         = harmony,
      nam_pcs         = nam_pcs,
      pipeline_output = pipeline_output)
}
