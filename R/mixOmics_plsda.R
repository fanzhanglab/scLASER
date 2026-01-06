#' PLS-DA step for scLASER objects using mixOmics
#'
#' Wrapper around [mixOmics::plsda()] that pulls the feature matrix from a
#' `scLASER` object, constructs class labels from metadata, optionally performs
#' multilevel decomposition, fits PLS-DA, and saves the latent variables
#' (component scores) into `object@plsda_LV`.
#'
#' Optionally, a subject-wise longitudinal cross-validation (CV) is performed,
#' where folds are created at the subject level (all visits for a subject are
#' kept in the same fold) and the Balanced Error Rate (BER) is computed for
#' 1..`cv_max_ncomp` components. CV results are attached as
#' `attr(object, "plsda_cv")`.
#'
#' @param object         A `scLASER` object.
#' @param response_var   Character. Column in `object@metadata` used as class
#'   label (default `"cell_type"`).
#' @param multilevel_var Optional character. Column in `object@metadata`
#'   identifying repeated-measures units (e.g. `"sample_id"`). If not `NULL`,
#'   multilevel PLS-DA is performed via [mixOmics::withinVariation()] prior
#'   to fitting.
#' @param ncomp          Integer. Number of latent components to compute for
#'   the main PLS-DA fit and for CV (upper bound).
#' @param cv_longitudinal Logical. If `TRUE`, perform subject-wise longitudinal
#'   CV as described above. Default is `FALSE`.
#' @param cv_subject_var Character. Column in `object@metadata` that defines
#'   the subject ID for CV folds (default `"subject_id"`).
#' @param cv_K           Integer. Number of folds for longitudinal CV
#'   (default `5`).
#' @param cv_max_ncomp   Integer. Maximum number of components to evaluate in
#'   longitudinal CV. Default is `ncomp`. Values greater than `ncomp` are
#'   truncated to `ncomp`.
#' @param cv_seed        Optional integer seed for reproducible fold creation.
#'   If `NULL`, the current RNG state is used.
#' @param ...            Additional arguments passed to [mixOmics::plsda()].
#'
#' @return A `scLASER` object with:
#'   \itemize{
#'     \item `@plsda_LV` populated with the PLS-DA latent variables
#'       (scores; `n_cells x ncomp`).
#'     \item If `cv_longitudinal = TRUE`, an attribute `plsda_cv` attached:
#'       `attr(object, "plsda_cv")`, a list with elements
#'       `ber` (numeric BER per component),
#'       `opt_ncomp` (optimal number of components),
#'       and `ncomp_seq` (sequence of components evaluated).
#'   }
#'
#' @export
#' @importFrom mixOmics plsda withinVariation
mixOmics_plsda <- function(object,
                           response_var    = "cell_type",
                           multilevel_var  = NULL,
                           ncomp           = 2,
                           cv_longitudinal = FALSE,
                           cv_subject_var  = "subject_id",
                           cv_K            = 5,
                           cv_max_ncomp    = ncomp,
                           cv_seed         = NULL,
                           ...) {

  stopifnot(inherits(object, "scLASER"))

  ## 1. Feature matrix from NAM_matrix ----
  X <- object@NAM_matrix
  if (is.null(X)) {
    stop(
      "object@NAM_matrix is NULL. ",
      "Please compute NAM_matrix before calling mixOmics_plsda()."
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

  ## 3. Optional multilevel design (withinVariation) ----
  if (!is.null(multilevel_var)) {
    if (!multilevel_var %in% colnames(meta)) {
      stop("multilevel_var '", multilevel_var,
           "' not found in object@metadata.")
    }
    design <- data.frame(sample = meta[[multilevel_var]])
    X <- withinVariation(X, design = design)
  }

  ## 4. Fit main PLS-DA model ----
  # Use mixOmics::plsda; ncomp is the max number of components we keep
  # fit <- plsda(X, Y, ncomp = ncomp, ...)

  # Latent variables (scores): n_cells x ncomp
  #LV <- fit$variates$X
  #object@plsda_LV <- as.matrix(LV)

  ## 5. Optional subject-wise longitudinal CV ----
  if (isTRUE(cv_longitudinal)) {
    if (!cv_subject_var %in% colnames(meta)) {
      stop("cv_subject_var '", cv_subject_var,
           "' not found in object@metadata.")
    }

    # Prepare subject-level table (unique (subject, class) combinations)
    subj_df <- meta[, c(cv_subject_var, response_var)]
    subj_df <- subj_df[!duplicated(subj_df), , drop = FALSE]

    # Randomize order of subjects for fold assignment
    if (!is.null(cv_seed)) {
      set.seed(cv_seed)
    }
    subj_df <- subj_df[sample(nrow(subj_df)), , drop = FALSE]

    # Build fold IDs, stratified by response class
    K <- cv_K
    fold_id <- integer(nrow(subj_df))
    for (cls in unique(subj_df[[response_var]])) {
      idx_cls <- which(subj_df[[response_var]] == cls)
      fold_id[idx_cls] <- rep(seq_len(K), length.out = length(idx_cls))
    }

    # List of subject IDs per fold
    subject_folds <- split(subj_df[[cv_subject_var]], fold_id)

    # Convert subject folds to row index folds (longitudinal, all visits)
    index_folds <- lapply(subject_folds, function(subj_grp) {
      which(meta[[cv_subject_var]] %in% subj_grp)
    })

    # Balanced Error Rate CV, similar to your manual code
    max_ncomp <- min(cv_max_ncomp, ncomp)
    n         <- nrow(X)
    classes   <- levels(Y)
    ber_cv    <- numeric(max_ncomp)

    for (nc in seq_len(max_ncomp)) {
      class_err_sum <- setNames(numeric(length(classes)), classes)
      class_n_sum   <- setNames(numeric(length(classes)), classes)

      for (fold in seq_along(index_folds)) {
        test_idx  <- index_folds[[fold]]
        train_idx <- setdiff(seq_len(n), test_idx)

        fit_fold <- plsda(
          X[train_idx, , drop = FALSE],
          Y[train_idx],
          ncomp = nc
        )

        pred <- predict(
          fit_fold,
          X[test_idx, , drop = FALSE],
          dist = "max.dist"
        )

        # Safely select available component
        ncomp_fit <- length(pred$class)
        ncomp_use <- min(nc, ncomp_fit)

        y_pred <- pred$class[[ncomp_use]]
        y_true <- Y[test_idx]

        for (cls in classes) {
          idx_cls_fold <- which(y_true == cls)
          if (length(idx_cls_fold) == 0) next

          err <- mean(y_pred[idx_cls_fold] != y_true[idx_cls_fold])
          class_err_sum[cls] <- class_err_sum[cls] + err * length(idx_cls_fold)
          class_n_sum[cls]   <- class_n_sum[cls]   + length(idx_cls_fold)
        }
      }

      class_err   <- class_err_sum / class_n_sum
      ber_cv[nc]  <- mean(class_err, na.rm = TRUE)
    }

    opt_ncomp <- which.min(ber_cv)

    fit <- plsda(X, Y, ncomp = opt_ncomp, ...)

    # Latent variables (scores): n_cells x ncomp
    LV <- fit$variates$X
    object@plsda_LV <- as.matrix(LV)
    return( object)
  } else{
     fit <- plsda(X, Y, ncomp = ncomp, ...)

    # Latent variables (scores): n_cells x ncomp
    LV <- fit$variates$X
    object@plsda_LV <- as.matrix(LV)
    return(object)
  }


}

