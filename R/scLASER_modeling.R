#' Fit longitudinal mixed models across all NAM PCs
#'
#' Fits mixed-effects models for each NAM PC stored in \code{object@nam_pcs}
#' using metadata from \code{object@metadata}. The function adapts to the number
#' of unique time points in \code{time_var}:
#'
#' \itemize{
#'   \item For \strong{2 time points}, fits two \code{nlme::lme()} models per PC:
#'   a full model with \code{disease_var * time_var} and a matched null without
#'   the interaction, then computes an LRT p-value.
#'   \item For \strong{3+ time points}, fits and compares three candidate models
#'   (continuous time, categorical time, quadratic time) by AIC, selects the best,
#'   and computes an LRT versus the matched null (no interaction for the chosen
#'   time specification).
#' }
#'
#' Results are saved to \code{object@pipeline_output} as a tidy data.frame.
#'
#' @param object A \code{\linkS4class{scLASER}} object with \code{@metadata} and
#'   \code{@nam_pcs}. Each NAM PC is modeled as an outcome.
#' @param disease_var Character scalar. Column name in \code{object@metadata}
#'   indicating group/disease status (coerced to factor internally).
#' @param time_var Character scalar. Column name in \code{object@metadata}
#'   indicating time/visit.
#' @param covariates Optional character vector of additional covariate column
#'   names in \code{object@metadata} to include as fixed effects (e.g.
#'   \code{c("age", "sex")}). Use \code{NULL} or \code{character(0)} to omit.
#' @param subject_var Character scalar. Column name identifying subjects for the
#'   level-1 random intercept.
#' @param sample_var Character scalar. Column name identifying samples nested
#'   within subjects for the level-2 random intercept.
#'
#' @return The input \code{scLASER} object with \code{object@pipeline_output}
#'   populated. The table includes, per PC: the selected model (3+ time points),
#'   fixed-effect coefficient rows, test statistics, and an LRT p-value for the
#'   interaction.
#'
#' @examples
#' \donttest{
#' ## Requires an scLASER object with metadata and NAM PCs:
#' ## obj <- readRDS("path/to/your_sclaser_object.rds")
#' ## obj <- scLASER_modeling(
#' ##   obj,
#' ##   disease_var = "disease",
#' ##   time_var    = "visit",
#' ##   covariates  = c("age", "sex"),
#' ##   subject_var = "subject_id",
#' ##   sample_var  = "sample_id"
#' ## )
#' ## head(obj@pipeline_output)
#' }
#'
#' @export
scLASER_modeling <- function(
  object,
  disease_var = "disease",
  time_var    = "visit",
  covariates  = c("age", "sex"),
  subject_var = "subject_id",
  sample_var  = "sample_id"
) {
  forAnalysis <- cbind(object@metadata, as.data.frame(object@nam_pcs))

  getLMM <- function(pc) {
    tmp <- forAnalysis
    tmp$nam_PC <- tmp[, which(colnames(tmp) == pc)]

    # --- build fixed effects using same arguments as 3+ timepoints ---
    covar_str <- if (length(covariates) && !is.null(covariates)) {
      paste(covariates, collapse = " + ")
    } else {
      NULL
    }
    rhs_plus_covars <- function(core_rhs) {
      if (is.null(covar_str)) core_rhs else paste(core_rhs, "+", covar_str)
    }

    rhs_full <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_var))
    rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_var))

    f_full <- as.formula(sprintf("nam_PC ~ %s", rhs_full))
    f_null <- as.formula(sprintf("nam_PC ~ %s", rhs_null))

    rand_f <- as.formula(sprintf("~ 1 | %s/%s", subject_var, sample_var))

    mod <- lme(
      fixed  = f_full,
      random = rand_f,
      data   = tmp,
      method = "ML"
    )

    null <- lme(
      fixed  = f_null,
      random = rand_f,
      data   = tmp,
      method = "ML"
    )

    full_tTab <- summary(mod)$tTable
    lrt <- anova(null, mod)
    lrt_pval <- lrt$`p-value`[2]
    want <- data.frame(
      NAMscore   = pc,
      coefficient = rownames(full_tTab),
      full_tTab,
      lrt_pval   = lrt_pval
    )
    return(want)
  }

  longAnalyses_plusModelSelection <- function(
    forAnalysis,
    outcome_var,
    disease_var,
    time_var,
    covariates = NULL,
    subject_var = "subject_id",
    sample_var = "sample_id"
  ) {
    stopifnot(is.data.frame(forAnalysis))
    tmp <- forAnalysis %>%
      dplyr::mutate(
        `:=`(!!disease_var, factor(.data[[disease_var]])),
        `:=`(!!time_var, .data[[time_var]])
      )

    time_fac_var <- paste0(time_var, "_fac__tmp")
    tmp[[time_fac_var]] <- factor(tmp[[time_var]])

    covar_str <- if (length(covariates) && !is.null(covariates)) {
      paste(covariates, collapse = " + ")
    } else {
      NULL
    }

    rhs_plus_covars <- function(core_rhs) {
      if (is.null(covar_str)) core_rhs else paste(core_rhs, "+", covar_str)
    }

    rhs_m1 <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_var))
    rhs_m2 <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_fac_var))
    rhs_m3 <- rhs_plus_covars(sprintf("%s * %s + %s * I(%s^2)",
                                      disease_var, time_var, disease_var, time_var))

    f_m1 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m1))
    f_m2 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m2))
    f_m3 <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_m3))

    rand_f <- as.formula(sprintf("~ 1 | %s/%s", subject_var, sample_var))

    m1 <- lme(fixed = f_m1, random = rand_f, data = tmp, method = "ML",
              control = lmeControl(opt = "optim"))
    m2 <- lme(fixed = f_m2, random = rand_f, data = tmp, method = "ML",
              control = lmeControl(opt = "optim"))
    m3 <- lme(fixed = f_m3, random = rand_f, data = tmp, method = "ML",
              control = lmeControl(opt = "optim"))

    fits <- list(m1 = m1, m2 = m2, m3 = m3)

    fit_tbl <- dplyr::tibble(
      model = names(fits),
      k     = vapply(fits, function(x) attr(logLik(x), "df"), numeric(1)),
      logLik = vapply(fits, function(x) as.numeric(logLik(x)), numeric(1)),
      AIC   = vapply(fits, AIC, numeric(1)),
      BIC   = vapply(fits, BIC, numeric(1))
    ) %>% dplyr::arrange(AIC)

    best_name <- fit_tbl$model[1]
    best_fit  <- fits[[best_name]]

    if (best_name == "m1") {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_var))
    } else if (best_name == "m2") {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_fac_var))
    } else {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s + I(%s^2)", disease_var, time_var, time_var))
    }

    f_null <- as.formula(sprintf("%s ~ %s", outcome_var, rhs_null))

    null_fit <- lme(fixed = f_null, random = rand_f, data = tmp, method = "ML",
                    control = lmeControl(opt = "optim"))

    lrt <- anova(null_fit, best_fit)

    want <- list(
      data_used = tmp,
      formulas  = list(m1 = f_m1, m2 = f_m2, m3 = f_m3, random = rand_f, null_for_best = f_null),
      fits      = fits,
      fit_stats = fit_tbl,
      best_model = list(name = best_name, fit = best_fit),
      lrt_best_vs_null = lrt
    )

    p <- want$lrt_best_vs_null$`p-value`[2]
    fitStats <- want$fit_stats
    tTab <- summary(want$best_model$fit)$tTable

    best_model_desc <- dplyr::case_when(
      want$best_model$name == "m1" ~ "time as continuous",
      want$best_model$name == "m2" ~ "time as categorical",
      want$best_model$name == "m3" ~ "quadratic time",
      TRUE ~ NA_character_
    )

    want_red <- data.frame(
      NAMscore   = outcome_var,
      best_model = best_model_desc,
      coefficient = rownames(tTab),
      tTab,
      lrt_pval   = want$lrt_best_vs_null$`p-value`[2],
      AIC        = fit_tbl$AIC[1]
    )

    rownames(want_red) <- seq_len(nrow(want_red))
    return(want_red)
  }

  n_time   <- length(unique(forAnalysis[[time_var]]))
  pc_names <- colnames(object@nam_pcs)

  if (n_time == 2) {
    out_list <- lapply(pc_names, getLMM)
  } else {
    out_list <- lapply(pc_names, function(pc) {
      longAnalyses_plusModelSelection(
        forAnalysis = forAnalysis,
        outcome_var = pc,
        disease_var = disease_var,
        time_var    = time_var,
        covariates  = covariates,
        subject_var = subject_var,
        sample_var  = sample_var
      )
    })
  }

  object@pipeline_output <- dplyr::bind_rows(out_list)
  return(object)
}
