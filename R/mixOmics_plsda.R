#' PLS-DA step for scLASER objects using mixOmics
#'
#' Wrapper around [mixOmics::plsda()] that pulls the feature matrix from a
#' `scLASER` object, constructs class labels from metadata, optionally performs
#' multilevel decomposition, fits PLS-DA, and saves the latent variables
#' (component scores) into `object@plsda_LV`.
#'
#' @param object        A `scLASER` object.
#' @param response_var  Character. Column in `object@metadata` used as class
#'   label (default `"cell_type"`).
#' @param multilevel_var Optional character. Column in `object@metadata`
#'   identifying repeated-measures units (e.g. `"sample_id"`). If not `NULL`,
#'   multilevel PLS-DA is performed via [mixOmics::withinVariation()].
#' @param ncomp         Integer. Number of latent components to compute.
#' @param ...           Additional arguments passed to [mixOmics::plsda()].
#'
#' @return A `scLASER` object with the `plsda_LV` slot populated with the PLS-DA
#'   latent variables (`n_cells x ncomp` matrix).
#'
#' @export
mixOmics_plsda <- function(object,
                          response_var   = "cell_type",
                          multilevel_var = NULL,
                          ncomp          = 2,
                          ...) {

  stopifnot(inherits(object, "scLASER"))

  ## 1. Feature matrix from NAM_matrix ----
  X <- object@NAM_matrix
  if (is.null(X)) {
    stop(
      "object@NAM_matrix is NULL. ",
      "Please compute NAM_matrix before calling scLASER_plsda()."
    )
  }
  X <- as.matrix(X)

  ## 2. Response from metadata ----
  meta <- object@metadata
  if (!response_var %in% colnames(meta)) {
    stop("response_var '", response_var, "' not found in object@metadata.")
  }
  Y <- factor(meta[[response_var]])

  if (nrow(X) != nrow(meta)) {
    stop(
      "Number of rows in NAM_matrix (", nrow(X),
      ") does not match number of rows in metadata (", nrow(meta), ")."
    )
  }

  ## 3. Optional multilevel design ----
  if (!is.null(multilevel_var)) {
    if (!multilevel_var %in% colnames(meta)) {
      stop("multilevel_var '", multilevel_var,
           "' not found in object@metadata.")
    }
    design <- data.frame(sample = meta[[multilevel_var]])
    X <- withinVariation(X, design = design)
  }

  ## 4. Fit PLS-DA ----
  fit <- plsda(X, Y, ncomp = ncomp, ...)

  ## 5. Save latent variables (scores) into LV slot ----
  # mixOmics convention: fit$variates$X is an n_cells x ncomp matrix
  LV <- fit$variates$X
  object@plsda_LV <- as.matrix(LV)


  object
}
