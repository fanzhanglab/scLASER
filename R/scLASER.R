setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL",     c("data.frame", "NULL"))

setClass(
  "scLASER",
  slots = list(
    metadata = "data.frame",
    pcs      = "matrixOrNULL",
    masc     = "dfOrNULL",
    umap     = "dfOrNULL",
    harmony  = "matrixOrNULL",
    nam_pcs  = "matrixOrNULL"
  ),
  prototype = list(
    metadata = data.frame(),
    pcs      = NULL,
    masc     = NULL,
    umap     = NULL,
    harmony  = NULL,
    nam_pcs  = NULL
  )
)

scLASER <- function(metadata = data.frame(),
                    pcs = NULL,
                    masc = NULL,
                    umap = NULL,
                    harmony = NULL,
                    nam_pcs = NULL) {
  new("scLASER",
      metadata = metadata,
      pcs      = pcs,
      masc     = masc,
      umap     = umap,
      harmony  = harmony,
      nam_pcs  = nam_pcs)
}
