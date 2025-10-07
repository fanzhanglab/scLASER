GLM_interact <- function(dataset, cluster, contrast1, contrast2, random_effects = NULL, fixed_effects = NULL,
         verbose = FALSE, save_models = FALSE, save_model_dir = NULL, save_name = NULL) {
  library(tidyverse)
  library(stevemisc)
  library(stevedata)
  library(lme4)
  library(broom.mixed)






  cluster <- as.character(cluster)
  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)


  cluster <- as.character(cluster)

  designmat <- model.matrix(~ cluster + 0, data.frame(cluster = cluster))
  dataset <- cbind(designmat, dataset)

  res <- vector(mode = "list", length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]


  if (!is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(c(paste0(fixed_effects, collapse = " + "),
                          paste0("(1|", random_effects, ")", collapse = " + "),
                          contrast1,
                          contrast2
    ),
    collapse = " + ")
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else if (!is.null(fixed_effects) && is.null(random_effects)) {
    model_rhs <- paste0(paste0(fixed_effects, collapse = " + "), " + ",
                        contrast1, " + ",
                        contrast2)
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))

      stop("No random effects specified")
    }
  } else if (is.null(fixed_effects) && !is.null(random_effects)) {
    model_rhs <- paste0(paste0("(1|", random_effects, ")", collapse = " + "), " + ",
                        contrast1, " + ",
                        contrast2)
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
    }
  } else {
    model_rhs <- paste0(contrast1, " + ",
                        contrast2)
    if (verbose == TRUE) {
      message(paste("Using null model:", "cluster ~", model_rhs))
      stop("No random or fixed effects specified")
    }
  }
  message(paste0("Using full model: cluster ~ ", model_rhs, " + ", contrast1, ":", contrast2 ))


  cluster_models <- vector(mode = "list",
                           length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]


  for (i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if (verbose == TRUE) {
      message(paste("Creating logistic mixed models for", test_cluster))
    }
    null_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ 1 + "),
                                   model_rhs), collapse = ""))
    full_fm <- as.formula(paste0(c(paste0(test_cluster, " ~ ", contrast1, ":", contrast2 , " + "),
                                   model_rhs), collapse = ""))

    null_model <- lme4::glmer(formula = null_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    full_model <- lme4::glmer(formula = full_fm, data = dataset,
                              family = binomial, nAGQ = 1, verbose = 0,
                              control = glmerControl(optimizer = "bobyqa"))
    model_lrt <- anova(null_model, full_model)

    contrast_lvl2 <- paste(paste0(contrast1, levels(dataset[[contrast1]])[2]),
                           paste0(contrast2, levels(dataset[[contrast2]])[2]),
                           sep=":")
    contrast_ci <- confint.merMod(full_model, method = "Wald",
                                  parm = contrast_lvl2)

    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt <- model_lrt
    cluster_models[[i]]$confint <- contrast_ci
  }


  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]],
                       size = colSums(designmat))
  output$model.pvalue <- sapply(cluster_models, function(x) x$model_lrt[["Pr(>Chisq)"]][2])
  output[[paste(contrast_lvl2, "OR", sep = ".")]] <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.lower", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "2.5 %"]))
  output[[paste(contrast_lvl2, "OR", "95pct.ci.upper", sep = ".")]] <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, "97.5 %"]))


  if (save_models == TRUE) {
    saveModelObj(cluster_models, save_dir = save_model_dir, save_name = save_name)
    return(output)
  } else {
    return(output)
  }
}
