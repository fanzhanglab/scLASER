
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL",     c("data.frame", "NULL"))

#' scLASER S4 class
#'
#' @slot metadata data.frame
#' @slot pcs matrix or NULL
#' @slot masc data.frame or NULL
#' @slot umap data.frame or NULL
#' @slot harmony matrix or NULL
#' @slot nam_pcs matrix or NULL
#' @slot NAM_matrix matrix or NULL
#' @slot pipeline_output data.frame or NULL
#'
#' @exportClass scLASER
#' @import methods
setClass(
  "scLASER",
  slots = list(
    metadata        = "data.frame",
    pcs             = "matrixOrNULL",
    masc            = "dfOrNULL",
    umap            = "dfOrNULL",
    harmony         = "matrixOrNULL",
    nam_pcs         = "matrixOrNULL",
    NAM_matrix      = "matrixOrNULL",
    pipeline_output = "dfOrNULL"
  ),
  prototype = list(
    metadata        = data.frame(),
    pcs             = NULL,
    masc            = NULL,
    umap            = NULL,
    harmony         = NULL,
    nam_pcs         = NULL,
    NAM_matrix      = NULL,
    pipeline_output = NULL
  )
)

#' Create a scLASER object
#' @export
scLASER <- function(metadata = data.frame(),
                    pcs = NULL,
                    masc = NULL,
                    umap = NULL,
                    harmony = NULL,
                    nam_pcs = NULL,
                    NAM_matrix = NULL,
                    pipeline_output = NULL) {
  methods::new("scLASER",
               metadata        = metadata,
               pcs             = pcs,
               masc            = masc,
               umap            = umap,
               harmony         = harmony,
               nam_pcs         = nam_pcs,
               NAM_matrix      = NAM_matrix,
               pipeline_output = pipeline_output)
}

