setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("dfOrNULL",     c("data.frame", "NULL"))

#' scLASER S4 class
#'
#' @slot metadata        data.frame
#' @slot pcs             matrix or NULL
#' @slot masc            data.frame or NULL
#' @slot umap            data.frame or NULL
#' @slot harmony         matrix or NULL
#' @slot nam_pcs         matrix or NULL
#' @slot NAM_matrix      matrix or NULL
#' @slot plsda_LV        matrix or NULL (latent variables from PLS-DA)
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
    plsda_LV        = "matrixOrNULL",
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
    plsda_LV        = NULL,
    pipeline_output = NULL
  )
)

#' Create a scLASER object
#'
#' Construct a \code{\linkS4class{scLASER}} object. This is a lightweight
#' constructor that stores provided components into the corresponding slots.
#'
#' @param metadata A data.frame of cell/sample metadata.
#' @param pcs A numeric matrix of principal components (rows aligned to
#'   \code{metadata}). Or \code{NULL}.
#' @param masc A data.frame of MASC/GLM results. Or \code{NULL}.
#' @param umap A data.frame with columns \code{UMAP1} and \code{UMAP2}. Or \code{NULL}.
#' @param harmony A numeric matrix of Harmony-corrected embeddings. Or \code{NULL}.
#' @param nam_pcs A numeric matrix of NAM PC embeddings. Or \code{NULL}.
#' @param NAM_matrix A numeric matrix storing the full NAM embedding matrix. Or \code{NULL}.
#' @param plsda_LV A numeric matrix of PLS-DA latent variables. Or \code{NULL}.
#' @param pipeline_output A data.frame storing tidy outputs from pipeline steps.
#'   Or \code{NULL}.
#'
#' @return A \code{\linkS4class{scLASER}} object.
#'
#' @examples
#' \donttest{
#' ## Minimal constructor usage:
#' ## obj <- scLASER(metadata = data.frame())
#' ## obj
#' }
#'
#' @export
scLASER <- function(metadata        = data.frame(),
                    pcs             = NULL,
                    masc            = NULL,
                    umap            = NULL,
                    harmony         = NULL,
                    nam_pcs         = NULL,
                    NAM_matrix      = NULL,
                    plsda_LV        = NULL,
                    pipeline_output = NULL) {
  methods::new(
    "scLASER",
    metadata        = metadata,
    pcs             = pcs,
    masc            = masc,
    umap            = umap,
    harmony         = harmony,
    nam_pcs         = nam_pcs,
    NAM_matrix      = NAM_matrix,
    plsda_LV        = plsda_LV,
    pipeline_output = pipeline_output
  )
}
