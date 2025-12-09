#' PLS-DA step for scLASER objects using mixOmics
#'
#' Convenience wrapper around [mixOmics::plsda()] that pulls the
#' feature matrix from a `scLASER` object, constructs the class labels
#' from a metadata column, optionally performs multilevel decomposition
#' with [mixOmics::withinVariation()], and returns the fitted PLS-DA model.
#'
#' @param object        A `scLASER` object.
#' @param response_var  Character. Name of the column in `object@metadata`
#'   to use as the class label (default `"cell_type"`).
#' @param use_nam       Logical. If `TRUE` (default), use `object@nam_pcs`;
#'   if `FALSE`, use `object@pcs` instead.
#' @param multilevel_var Optional character. Name of the metadata column
#'   identifying the unit for repeated measures (e.g. `"sample_id"`).
#'   If not `NULL`, multilevel PLS-DA is performed by applying
#'   [mixOmics::withinVariation()] before `plsda()`.
#' @param ncomp         Integer. Number of latent components to compute.
#'   Passed to [mixOmics::plsda()].
#' @param ...           Additional arguments passed on to
#'   [mixOmics::plsda()].
#'
#' @return A `mixOmics` `plsda` object.
#'
#' @export
#' @importFrom mixOmics plsda withinVariation
mixOmics_plsda <- function(object,
                          response_var   = "cell_type",
                          use_nam        = TRUE,
                          multilevel_var = NULL,
                          ncomp          = 2,
                          ...) {

  # Basic sanity check
  stopifnot(inherits(object, "scLASER"))

  ## 1. Choose feature matrix: NAM PCs or regular PCs ----
  if (use_nam) {
    X <- object@nam_pcs
    if (is.null(X)) {
      stop(
        "object@nam_pcs is NULL but `use_nam = TRUE`. ",
        "Run the NAM step first or set `use_nam = FALSE`."
      )
    }
  } else {
    X <- object@pcs
    if (is.null(X)) {
      stop(
        "object@pcs is NULL and `use_nam = FALSE`. ",
        "Provide PCs or use_nam = TRUE."
      )
    }
  }

  X <- as.matrix(X)

  ## 2. Build response vector from metadata ----
  meta <- object@metadata
  if (!response_var %in% colnames(meta)) {
    stop("response_var '", response_var, "' not found in object@metadata.")
  }

  Y <- factor(meta[[response_var]])

  ## 3. Optional multilevel PLS-DA (within-subject variation) ----
  if (!is.null(multilevel_var)) {
    if (!multilevel_var %in% colnames(meta)) {
      stop("multilevel_var '", multilevel_var,
           "' not found in object@metadata.")
    }

    design <- data.frame(sample = meta[[multilevel_var]])
    # mixOmics::withinVariation
    X <- withinVariation(X, design = design)
  }

  ## 4. Fit PLS-DA via mixOmics ----
  # mixOmics::plsda
  fit <- plsda(X, Y, ncomp = ncomp, ...)

  fit$variates$X

}
