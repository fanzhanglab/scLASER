#' Run longitudinal mixed-model selection for all NAM PCs in a scLASER object
#'
#' Pulls metadata and @nam_pcs from `object`, fits three mixed models per PC
#' (continuous time, categorical time, quadratic time), selects best by AIC,
#' computes an LRT vs the matched null (no interaction), and saves a tidy table
#' into `object@pipeline_output`.
#'
#' @param object      scLASER object (requires @metadata and @nam_pcs)
#' @param disease_var Column name in metadata for disease/group (will be factored)
#' @param time_var    Column name in metadata for time/visit
#' @param covariates  Optional character vector of additional covariate names
#' @param subject_var Column name for subject ID (random intercept level 1)
#' @param sample_var  Column name for sample ID (random intercept level 2)
#' @param na_action   One of "omit" (default) or "fail" for handling NA rows
#'
#' @return The modified scLASER object with `@pipeline_output` set to a data.frame.
#' @export
setGeneric("longAnalyses_plusModelSelection",
           function(object, ...) standardGeneric("longAnalyses_plusModelSelection"))

setMethod(
  "longAnalyses_plusModelSelection",
  signature(object = "scLASER"),
  function(
    object,
    disease_var,
    time_var,
    covariates = NULL,
    subject_var = "subject_id",
    sample_var  = "sample_id",
    na_action   = c("omit", "fail")
  ){
    na_action <- match.arg(na_action)

    ## --- sanity checks ---
    if (!is.data.frame(object@metadata) || nrow(object@metadata) == 0)
      stop("@metadata is missing or empty in scLASER object.")

    if (is.null(object@nam_pcs) || !is.matrix(object@nam_pcs))
      stop("@nam_pcs is NULL or not a matrix; run association_nam_scLASER() first.")

    meta <- object@metadata
    pcs  <- object@nam_pcs

    need <- c(disease_var, time_var, subject_var, sample_var, covariates)
    need <- need[!is.null(need)]
    miss <- setdiff(need, names(meta))
    if (length(miss))
      stop("Missing columns in metadata: ", paste(miss, collapse = ", "))

    if (is.null(rownames(pcs)))
      stop("rownames(@nam_pcs) are NULL; they must be sample IDs matching '", sample_var, "'.")

    if (!all(meta[[sample_var]] %in% rownames(pcs))) {
      bad <- setdiff(unique(meta[[sample_var]]), rownames(pcs))
      stop("Some ", sample_var, " not found in rownames(@nam_pcs). Example: ",
           paste(head(bad, 5), collapse = ", "))
    }

    ## --- base data prep ---
    tmp_base <- meta

    # disease as factor
    tmp_base[[disease_var]] <- factor(tmp_base[[disease_var]])

    # time: ensure numeric for continuous model; also keep a factor copy
    # if not numeric, map ordered/factor levels to 0..K-1
    if (!is.numeric(tmp_base[[time_var]])) {
      tf <- factor(tmp_base[[time_var]])
      tmp_base[[time_var]] <- as.numeric(tf) - 1L
    }
    time_fac_var <- paste0(time_var, "_fac__tmp")
    tmp_base[[time_fac_var]] <- factor(meta[[time_var]])

    # assemble covariates
    covar_str <- if (length(covariates) && !is.null(covariates)) paste(covariates, collapse = " + ") else NULL
    rhs_plus_covars <- function(core_rhs) if (is.null(covar_str)) core_rhs else paste(core_rhs, "+", covar_str)

    # three fixed-effects templates
    rhs_m1_tpl <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_var))                 # continuous time
    rhs_m2_tpl <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_fac_var))            # categorical time
    rhs_m3_tpl <- rhs_plus_covars(sprintf("%s * %s + %s * I(%s^2)",
                                          disease_var, time_var, disease_var, time_var))    # quadratic time

    # random effects (nested: subject/sample)
    rand_f <- stats::as.formula(sprintf("~ 1 | %s/%s", subject_var, sample_var))
    ctrl   <- nlme::lmeControl(opt = "optim", returnObject = TRUE, msMaxIter = 200)





    # small helper: safe lme fit
    fit_safe <- function(fml, data) {
      tryCatch({
        fit <- nlme::lme(
          fixed   = fml,
          random  = rand_f,
          data    = data,
          method  = "ML",
          control = ctrl
        )
        ## IMPORTANT: make sure the stored call holds the *formula*,
        ## not the symbol "fml"
        fit$call$fixed <- fml
        fit
      }, error = function(e) NULL)
    }

    # optionally drop rows with NAs across needed columns for modeling
    needed_for_model <- unique(c(disease_var, time_var, time_fac_var, subject_var, sample_var, covariates))
    if (na_action == "omit") {
      keep <- stats::complete.cases(tmp_base[, needed_for_model, drop = FALSE])
      tmp_base <- tmp_base[keep, , drop = FALSE]
    } else if (anyNA(tmp_base[, needed_for_model, drop = FALSE])) {
      stop("NA present in required columns and na_action='fail'.")
    }

    pc_names <- colnames(pcs)

    fit_one_pc <- function(pc_name) {
      # align outcome vector by sample_id order in tmp_base
      out_vec <- pcs[ tmp_base[[sample_var]], pc_name ]
      out_col <- ".tmp_outcome__NAMPC"
      tmp <- tmp_base
      tmp[[out_col]] <- out_vec

      # formulas
      f_m1 <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m1_tpl))
      f_m2 <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m2_tpl))
      f_m3 <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m3_tpl))

      m1 <- fit_safe(f_m1, tmp)
      m2 <- fit_safe(f_m2, tmp)
      m3 <- fit_safe(f_m3, tmp)

      fits <- Filter(Negate(is.null), list(m1 = m1, m2 = m2, m3 = m3))
      if (length(fits) == 0L) {
        return(data.frame(
          NAMscore    = pc_name,
          best_model  = NA_character_,
          coefficient = NA_character_,
          t.value     = NA_real_,
          `p-value`   = NA_real_,
          lrt_pval    = NA_real_,
          AIC         = NA_real_,
          row.names   = NULL,
          check.names = FALSE
        ))
      }

      fit_tbl <- dplyr::tibble(
        model  = names(fits),
        k      = vapply(fits, function(x) attr(stats::logLik(x), "df"), numeric(1)),
        logLik = vapply(fits, function(x) as.numeric(stats::logLik(x)), numeric(1)),
        AIC    = vapply(fits, stats::AIC, numeric(1)),
        BIC    = vapply(fits, stats::BIC, numeric(1))
      ) |> dplyr::arrange(.data$AIC)

      best_name <- fit_tbl$model[1]
      best_fit  <- fits[[best_name]]

      # matched null: drop interaction, keep main effects (+ quadratic term for m3)
      if (best_name == "m1") {
        rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_var))
      } else if (best_name == "m2") {
        rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_fac_var))
      } else { # m3
        rhs_null <- rhs_plus_covars(sprintf("%s + %s + I(%s^2)", disease_var, time_var, time_var))
      }
      f_null   <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_null))
      null_fit <- fit_safe(f_null, tmp)

      lrt <- if (!is.null(null_fit)) stats::anova(null_fit, best_fit) else NA
      lrt_p <- if (is.data.frame(lrt) && nrow(lrt) >= 2 && "p-value" %in% colnames(lrt)) lrt$`p-value`[2] else NA_real_

      best_model_desc <- dplyr::case_when(
        best_name == "m1" ~ "time as continuous",
        best_name == "m2" ~ "time as categorical",
        best_name == "m3" ~ "quadratic time",
        TRUE ~ NA_character_
      )

      tTab <- summary(best_fit)$tTable

      out <- data.frame(
        NAMscore    = pc_name,
        best_model  = best_model_desc,
        coefficient = rownames(tTab),
        tTab,
        lrt_pval    = lrt_p,
        AIC         = fit_tbl$AIC[1],
        row.names   = NULL,
        check.names = FALSE
      )

      # add a few breadcrumbs to help downstream filtering
      out$subject_var <- subject_var
      out$sample_var  <- sample_var
      out$disease_var <- disease_var
      out$time_var    <- time_var

      out
    }

    results_list <- lapply(pc_names, fit_one_pc)
    results_all  <- do.call(rbind, results_list)

    # save into the object
    object@pipeline_output <- results_all
    return(object)
  }
)
