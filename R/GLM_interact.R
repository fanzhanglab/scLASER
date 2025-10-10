GLM_interact <- function(object,
                         cluster_col,            # e.g. "cell_type"
                         contrast1, contrast2,   # e.g. "visit", "disease"
                         random_effects = NULL,  # e.g. "sample_id"
                         fixed_effects  = NULL,
                         verbose = FALSE,
                         save_models = FALSE,
                         save_model_dir = NULL,
                         save_name = NULL) {

  stopifnot(methods::is(object, "scLASER"))
  if (!requireNamespace("lme4", quietly = TRUE))
    stop("Package 'lme4' is required.")

  ## -- copy metadata and coerce ----------------------------------------------
  dataset2 <- as.data.frame(object@metadata, stringsAsFactors = FALSE)

  # make IDs factors (prevents lme4/Matrix crashes)
  if ("sample_id"  %in% names(dataset2)) dataset2$sample_id  <- factor(dataset2$sample_id)
  if ("subject_id" %in% names(dataset2)) dataset2$subject_id <- factor(dataset2$subject_id)

  need_cols <- c(cluster_col, contrast1, contrast2, if (!is.null(random_effects)) random_effects)
  miss <- setdiff(need_cols, names(dataset2))
  if (length(miss)) stop("Missing required columns in metadata: ", paste(miss, collapse = ", "))

  # make contrasts factors with treatment coding
  dataset2[[contrast1]] <- droplevels(factor(dataset2[[contrast1]]))
  dataset2[[contrast2]] <- droplevels(factor(dataset2[[contrast2]]))
  stats::contrasts(dataset2[[contrast1]]) <- stats::contr.treatment(nlevels(dataset2[[contrast1]]))
  stats::contrasts(dataset2[[contrast2]]) <- stats::contr.treatment(nlevels(dataset2[[contrast2]]))

  # drop incomplete rows
  dataset2 <- dataset2[stats::complete.cases(dataset2[, need_cols, drop = FALSE]), , drop = FALSE]

  # random effect as factor with >=2 levels (if used)
  if (!is.null(random_effects)) {
    dataset2[[random_effects]] <- droplevels(factor(dataset2[[random_effects]]))
    if (nlevels(dataset2[[random_effects]]) < 2)
      stop("random_effects '", random_effects, "' has <2 levels after filtering")
  }

  ## -- build cluster 1-vs-rest responses -------------------------------------
  cl_fac     <- droplevels(factor(dataset2[[cluster_col]]))
  designmat  <- stats::model.matrix(~ cl_fac + 0)
  colnames(designmat) <- sub("^cl_fac", "cluster", colnames(designmat))
  clust_cols <- colnames(designmat)
  dataset2   <- cbind.data.frame(as.data.frame(designmat, check.names = FALSE), dataset2)

  ## -- model RHS like your original ------------------------------------------
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + "),
                          contrast1, contrast2), collapse = " + ")
    if (verbose) message("Using null model: cluster ~ ", model_rhs)
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(paste0(fixed_effects, collapse = " + "), " + ",
                        contrast1, " + ", contrast2)
    if (verbose) { message("Using null model: cluster ~ ", model_rhs); stop("No random effects specified") }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(paste0("(1|", random_effects, ")", collapse = " + "), " + ",
                        contrast1, " + ", contrast2)
    if (verbose) message("Using null model: cluster ~ ", model_rhs)
  } else {
    model_rhs <- paste0(contrast1, " + ", contrast2)
    if (verbose) { message("Using null model: cluster ~ ", model_rhs); stop("No random or fixed effects specified") }
  }
  message("Using full model: cluster ~ ", model_rhs, " + ", contrast1, ":", contrast2)

  ## -- containers -------------------------------------------------------------
  cluster_models <- vector("list", length(clust_cols))
  names(cluster_models) <- clust_cols

  coef_name <- paste(
    paste0(contrast1, levels(dataset2[[contrast1]])[min(2, nlevels(dataset2[[contrast1]]))]),
    paste0(contrast2, levels(dataset2[[contrast2]])[min(2, nlevels(dataset2[[contrast2]]))]),
    sep = ":"
  )

  fit_glmer <- function(form) {
    lme4::glmer(formula = form, data = dataset2,
                family = stats::binomial(), nAGQ = 1, verbose = 0,
                control = lme4::glmerControl(optimizer = "bobyqa"))
  }

  ## -- fit per cluster with fallback to glm ----------------------------------
  for (i in seq_along(clust_cols)) {
    test_cluster <- clust_cols[i]
    y <- dataset2[[test_cluster]]

    # skip degenerate responses
    if (all(y %in% c(0, 1)) && (all(y == 0) || all(y == 1))) {
      if (verbose) message("Skipping ", test_cluster, " (all zeros or all ones).")
      cluster_models[[i]] <- list(null_model = NULL, full_model = NULL,
                                  model_lrt = NA, confint = NA, method = NA)
      next
    }

    if (verbose) message("Fitting GLMMs for ", test_cluster)

    null_fm <- stats::as.formula(paste0(test_cluster, " ~ 1 + ", model_rhs))
    full_fm <- stats::as.formula(paste0(test_cluster, " ~ ", contrast1, ":", contrast2, " + ", model_rhs))

    res <- try({
      null_model <- fit_glmer(null_fm)
      full_model <- fit_glmer(full_fm)
      model_lrt  <- stats::anova(null_model, full_model)
      ci         <- try(lme4::confint.merMod(full_model, method = "Wald", parm = coef_name), silent = TRUE)
      list(null_model = null_model, full_model = full_model,
           model_lrt = model_lrt, confint = ci, method = "glmer")
    }, silent = TRUE)

    # Fallback: plain logistic regression (no random effects)
    if (inherits(res, "try-error")) {
      if (verbose) message("GLMM failed for ", test_cluster, " — trying glm fallback.")
      fe_rhs   <- gsub("\\(1\\|[^)]+\\) \\+ ?", "", model_rhs) # strip (1|re) if present
      null_fe  <- stats::as.formula(paste0(test_cluster, " ~ 1 + ", fe_rhs))
      full_fe  <- stats::as.formula(paste0(test_cluster, " ~ ", contrast1, ":", contrast2, " + ", fe_rhs))
      res2 <- try({
        null_m <- stats::glm(null_fe, data = dataset2, family = stats::binomial())
        full_m <- stats::glm(full_fe, data = dataset2, family = stats::binomial())
        lr     <- stats::anova(null_m, full_m, test = "LRT")
        bb     <- stats::coef(full_m)
        ci     <- try(stats::confint.default(full_m), silent = TRUE)  # Wald-style CI
        list(null_model = null_m, full_model = full_m,
             model_lrt = lr, confint = ci, method = "glm")
      }, silent = TRUE)

      if (inherits(res2, "try-error")) {
        if (verbose) message("Fallback glm also failed for ", test_cluster, ": ", as.character(res2))
        cluster_models[[i]] <- list(null_model = NULL, full_model = NULL,
                                    model_lrt = NA, confint = NA, method = NA)
      } else {
        cluster_models[[i]] <- res2
      }
    } else {
      cluster_models[[i]] <- res
    }
  }

  ## -- assemble MASC-like output ---------------------------------------------
  out <- data.frame(cluster = clust_cols,
                    size    = colSums(designmat),
                    stringsAsFactors = FALSE)

  out$model.pvalue <- vapply(cluster_models, function(x) {
    if (is.null(x$model_lrt) || all(is.na(x$model_lrt))) return(NA_real_)
    # glmer anova or glm anova both put the test in row 2
    suppressWarnings(x$model_lrt[["Pr(>Chisq)"]][2])
  }, numeric(1))

  or_name <- paste(coef_name, "OR", sep = ".")
  lo_name <- paste(coef_name, "OR", "95pct.ci.lower", sep = ".")
  hi_name <- paste(coef_name, "OR", "95pct.ci.upper", sep = ".")

  out[[or_name]] <- vapply(cluster_models, function(x) {
    if (is.null(x$full_model)) return(NA_real_)
    if (identical(x$method, "glmer")) {
      bb <- lme4::fixef(x$full_model)
      if (!coef_name %in% names(bb)) return(NA_real_)
      return(exp(bb[[coef_name]]))
    } else if (identical(x$method, "glm")) {
      bb <- stats::coef(x$full_model)
      if (!coef_name %in% names(bb)) return(NA_real_)
      return(exp(bb[[coef_name]]))
    }
    NA_real_
  }, numeric(1))

  out[[lo_name]] <- vapply(cluster_models, function(x) {
    if (is.null(x$confint) || inherits(x$confint, "try-error")) return(NA_real_)
    if (identical(x$method, "glmer")) {
      rn <- rownames(x$confint); if (is.null(rn) || !(coef_name %in% rn)) return(NA_real_)
      return(exp(x$confint[coef_name, "2.5 %"]))
    } else if (identical(x$method, "glm")) {
      cn <- colnames(x$confint); rn <- rownames(x$confint)
      if (is.null(rn) || is.null(cn) || !(coef_name %in% rn)) return(NA_real_)
      return(exp(x$confint[coef_name, "2.5 %"]))
    }
    NA_real_
  }, numeric(1))

  out[[hi_name]] <- vapply(cluster_models, function(x) {
    if (is.null(x$confint) || inherits(x$confint, "try-error")) return(NA_real_)
    if (identical(x$method, "glmer")) {
      rn <- rownames(x$confint); if (is.null(rn) || !(coef_name %in% rn)) return(NA_real_)
      return(exp(x$confint[coef_name, "97.5 %"]))
    } else if (identical(x$method, "glm")) {
      cn <- colnames(x$confint); rn <- rownames(x$confint)
      if (is.null(rn) || is.null(cn) || !(coef_name %in% rn)) return(NA_real_)
      return(exp(x$confint[coef_name, "97.5 %"]))
    }
    NA_real_
  }, numeric(1))

  if (isTRUE(save_models)) {
    if (is.null(save_model_dir)) stop("save_model_dir must be provided when save_models=TRUE")
    if (!dir.exists(save_model_dir)) dir.create(save_model_dir, recursive = TRUE)
    saveRDS(cluster_models, file = file.path(
      save_model_dir, paste0(ifelse(is.null(save_name), "masc_models", save_name), ".rds")
    ))
  }

  object@masc <- out
  object
}

