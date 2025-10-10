GLM_interact <- function(object,                      # scLASER object
                         cluster_col,                 # column name in metadata
                         contrast1, contrast2,
                         random_effects = NULL, fixed_effects = NULL,
                         verbose = FALSE, save_models = FALSE,
                         save_model_dir = NULL, save_name = NULL) {

  # libs as in your original
  library(tidyverse)
  library(stevemisc)
  library(stevedata)
  library(lme4)
  library(broom.mixed)

  # pull metadata
  dataset <- object@metadata
  if (!is.data.frame(dataset) || nrow(dataset) == 0)
    stop("object@metadata is empty or not a data.frame")
  if (!all(c(cluster_col, contrast1, contrast2) %in% names(dataset))) {
    stop("Missing required columns in metadata: ",
         paste(setdiff(c(cluster_col, contrast1, contrast2), names(dataset)), collapse = ", "))
  }

  # build 1-vs-rest design per cluster
  cluster   <- as.character(dataset[[cluster_col]])
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset   <- cbind(designmat, dataset)

  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- colnames(designmat)

  # RHS (kept)
  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + "),
                          contrast1, contrast2), collapse = " + ")
    if (verbose) message(paste("Using null model:", "cluster ~", model_rhs))
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(paste0(fixed_effects, collapse = " + "), " + ",
                        contrast1, " + ", contrast2)
    if (verbose) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(paste0("(1|", random_effects, ")", collapse = " + "), " + ",
                        contrast1, " + ", contrast2)
    if (verbose) message(paste("Using null model:", "cluster ~", model_rhs))
  } else {
    model_rhs <- paste0(contrast1, " + ", contrast2)
    if (verbose) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  message(paste0("Using full model: cluster ~ ", model_rhs, " + ", contrast1, ":", contrast2))

  cluster_models <- vector("list", length = ncol(designmat))
  names(cluster_models) <- colnames(designmat)

  for (i in seq_along(colnames(designmat))) {
    test_cluster <- colnames(designmat)[i]
    if (verbose) message(paste("Creating logistic mixed models for", test_cluster))

    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast1, ":", contrast2 , " + "),
                                   model_rhs), collapse = ""))

    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = lme4::glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)

    contrast_lvl2 <- paste(paste0(contrast1, levels(dataset[[contrast1]])[2]),
                           paste0(contrast2, levels(dataset[[contrast2]])[2]),
                           sep=":")
    contrast_ci <- confint.merMod(full_model, method = "Wald", parm = contrast_lvl2)

    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt  <- model_lrt
    cluster_models[[i]]$confint    <- contrast_ci
  }

  output <- data.frame(cluster = colnames(designmat),
                       size    = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])

  # use full_model (your original had x$full which errors)
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <-
    sapply(cluster_models, function(x) exp(lme4::fixef(x$full_model)[[contrast_lvl2]]))

  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <-
    sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))

  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <-
    sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))

  if (save_models) {
    if (is.null(save_model_dir)) stop("save_model_dir must be provided when save_models=TRUE")
    if (!dir.exists(save_model_dir)) dir.create(save_model_dir, recursive = TRUE)
    saveRDS(cluster_models, file = file.path(save_model_dir, paste0(ifelse(is.null(save_name), "masc_models", save_name), ".rds")))
  }

  # save into the object and return it
  object@masc <- output
  object
}

