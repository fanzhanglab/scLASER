longAnalyses_plusModelSelection <- function(
    forAnalysis,
    outcome_var, # e.g., "NAM_PC1"
    disease_var, # e.g., "disease"
    time_var, # e.g., "visit"
    covariates = NULL,  # e.g., c("age","sex")
    subject_var = "subject_id",
    sample_var  = "sample_id"
){

  #check that data is a data.frame
  stopifnot(is.data.frame(forAnalysis))
  # make a "tmp" dataset and ensure types for disease and time factor
  tmp <- forAnalysis %>%
    mutate(
      !!disease_var := factor(.data[[disease_var]]),
      !!time_var    := .data[[time_var]]
    )
  # make a factor of time
  time_fac_var <- paste0(time_var, "_fac__tmp")
  tmp[[time_fac_var]] <- factor(tmp[[time_var]])

  # utility: collapse covariates into formula string
  covar_str <- if (length(covariates) && !is.null(covariates)) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }
  rhs_plus_covars <- function(core_rhs) {
    if (is.null(covar_str)) core_rhs else paste(core_rhs, "+", covar_str)
  }

  # Fixed-effect RHS for each model
  rhs_m1 <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_var))
  rhs_m2 <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_fac_var))
  rhs_m3 <- rhs_plus_covars(sprintf("%s * %s + %s * I(%s^2)",
                                    disease_var, time_var, disease_var, time_var))

  # Build full fixed formulas
  f_m1 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m1))
  f_m2 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m2))
  f_m3 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m3))

  # Random-effects structure (nested: subject/sample)
  rand_f <- as.formula(sprintf("~ 1 | %s/%s", subject_var, sample_var))

  # Fit all three with ML for fair comparison
  m1 <- lme(fixed = f_m1, random = rand_f, data = tmp, method = "ML", control = lmeControl(opt = "optim"))
  m2 <- lme(fixed = f_m2, random = rand_f, data = tmp, method = "ML", control = lmeControl(opt = "optim"))
  m3 <- lme(fixed = f_m3, random = rand_f, data = tmp, method = "ML", control = lmeControl(opt = "optim"))

  # Compare fits
  fits <- list(m1 = m1, m2 = m2, m3 = m3)
  fit_tbl <- dplyr::tibble(
    model = names(fits),
    k     = vapply(fits, function(x) attr(logLik(x), "df"), numeric(1)),
    logLik= vapply(fits, function(x) as.numeric(logLik(x)), numeric(1)),
    AIC   = vapply(fits, AIC, numeric(1)),
    BIC   = vapply(fits, BIC, numeric(1))
  ) %>% arrange(AIC)

  best_name <- fit_tbl$model[1]
  best_fit  <- fits[[best_name]]

  # Build null models (remove the disease:time interaction but keep main effects)
  # Note: we build null formulas explicitly for clarity/stability.
  if (best_name == "m1") {
    rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_var))
  } else if (best_name == "m2") {
    rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_fac_var))
  } else { # m3
    rhs_null <- rhs_plus_covars(sprintf("%s + %s + I(%s^2)", disease_var, time_var, time_var))
  }
  f_null <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_null))
  null_fit <- lme(fixed = f_null, random = rand_f, data = tmp, method = "ML", control = lmeControl(opt = "optim"))

  # LRT for interaction
  lrt <- anova(null_fit, best_fit)

  want <- list(
    data_used   = tmp,
    formulas    = list(m1 = f_m1, m2 = f_m2, m3 = f_m3, random = rand_f, null_for_best = f_null),
    fits        = fits,
    fit_stats   = fit_tbl,
    best_model  = list(name = best_name, fit = best_fit),
    lrt_best_vs_null = lrt
  )

  #### Subset to data we want for ####
  p <- want$lrt_best_vs_null$`p-value`[2]
  # Fit statistics (this table is sorted by AIC, lowest value (i.e. best model) at top, but contains other info like BIC, logLikelihood and the number of coefficients)
  fitStats <- want$fit_stats
  # Fixed Effects T-table for Effect Sizes (the structure of this may change depending on the best fit model)
  tTab <- summary(want$best_model$fit)$tTable
  # Description of the best fit model
  best_model_desc <- case_when(
    want$best_model$name == "m1" ~ "time as continuous",
    want$best_model$name == "m2" ~ "time as categorical",
    want$best_model$name == "m3" ~ "quadratic time",
    TRUE ~ NA_character_
  )

  want_red <- data.frame(NAMscore = outcome_var, best_model = best_model_desc,coefficient = rownames(tTab), tTab, lrt_pval = want$lrt_best_vs_null$`p-value`[2], AIC = fit_tbl$AIC[1])
  rownames(want_red) <- c(1:nrow(want_red))
  return(want_red)
}
