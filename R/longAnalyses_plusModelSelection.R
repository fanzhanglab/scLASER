#' Title
#'
#' @param object 
#' @param disease_var 
#' @param time_var 
#' @param covariates 
#' @param subject_var 
#' @param sample_var 
#'
#' @return
#' @export
#'
#' @examples
longAnalyses_plusModelSelection_scLASER <- function(
  object,
  disease_var,        # e.g., "disease"
  time_var,           # e.g., "visit"
  covariates = NULL,  # e.g., c("age","sex")
  subject_var = "subject_id",
  sample_var  = "sample_id"
){
  stopifnot(inherits(object, "scLASER"))
  if (is.null(object@nam_pcs) || !is.matrix(object@nam_pcs)) {
    stop("@nam_pcs is NULL or not a matrix; run association_nam_scLASER() first.")
  }
  if (!is.data.frame(object@metadata) || nrow(object@metadata) == 0) {
    stop("@metadata is missing or empty in scLASER object.")
  }
  
  meta <- object@metadata
  pcs  <- object@nam_pcs
  
  need <- c(disease_var, time_var, subject_var, sample_var, covariates)
  miss <- setdiff(need[!is.null(need)], names(meta))
  if (length(miss)) stop("Missing columns in metadata: ", paste(miss, collapse = ", "))
  
  if (is.null(rownames(pcs))) {
    stop("rownames(@nam_pcs) are NULL; they must be sample IDs matching '", sample_var, "'.")
  }
  if (!all(meta[[sample_var]] %in% rownames(pcs))) {
    bad <- setdiff(unique(meta[[sample_var]]), rownames(pcs))
    stop("Some ", sample_var, " not found in rownames(@nam_pcs). Example: ", paste(head(bad, 5), collapse = ", "))
  }
  
  # --- prep types (key tweak: continuous time must be numeric; for 2TP ensure 0/1) ---
  tmp_base <- meta
  tmp_base[[disease_var]] <- factor(tmp_base[[disease_var]])
  # continuous time (numeric)
  if (!is.numeric(tmp_base[[time_var]])) {
    tf <- factor(tmp_base[[time_var]])
    # map to 0,1,... (for 2TP this becomes 0/1 automatically)
    tmp_base[[time_var]] <- as.numeric(tf) - 1L
  }
  # categorical time (factor copy)
  time_fac_var <- paste0(time_var, "_fac__tmp")
  tmp_base[[time_fac_var]] <- factor(meta[[time_var]])
  
  # covariates
  covar_str <- if (length(covariates) && !is.null(covariates)) paste(covariates, collapse = " + ") else NULL
  rhs_plus_covars <- function(core_rhs) if (is.null(covar_str)) core_rhs else paste(core_rhs, "+", covar_str)
  
  # 3TP logic (works for 2TP as well)
  rhs_m1_tpl <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_var))                 # continuous time
  rhs_m2_tpl <- rhs_plus_covars(sprintf("%s * %s", disease_var, time_fac_var))            # categorical time
  rhs_m3_tpl <- rhs_plus_covars(sprintf("%s * %s + %s * I(%s^2)",                         # quadratic time
                                        disease_var, time_var, disease_var, time_var))
  
  # random effects (same as yours)
  rand_f <- stats::as.formula(sprintf("~ 1 | %s/%s", subject_var, sample_var))
  ctrl   <- nlme::lmeControl(opt = "optim", msMaxIter = 200, returnObject = TRUE)
  
  # safe fitter that returns NULL on any fitting/singularity issue
  fit_safe <- function(fml, data) {
    tryCatch(
      nlme::lme(fixed = fml, random = rand_f, data = data, method = "ML", control = ctrl),
      error = function(e) NULL
    )
  }
  
  pc_names <- colnames(pcs)
  
  fit_one_pc <- function(pc_name) {
    out_vec <- pcs[ tmp_base[[sample_var]], pc_name ]
    out_col <- ".tmp_outcome__NAMPC"
    tmp <- tmp_base
    tmp[[out_col]] <- out_vec
    
    # build formulas
    f_m1   <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m1_tpl))
    f_m2   <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m2_tpl))
    f_m3   <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_m3_tpl))
    
    m1 <- fit_safe(f_m1, tmp)
    m2 <- fit_safe(f_m2, tmp)
    m3 <- fit_safe(f_m3, tmp)
    
    fits <- Filter(Negate(is.null), list(m1 = m1, m2 = m2, m3 = m3))
    if (length(fits) == 0L) {
      return(data.frame(
        NAMscore = pc_name, best_model = NA, coefficient = NA,
        t.value = NA, `p-value` = NA, lrt_pval = NA, AIC = NA,
        row.names = NULL, check.names = FALSE
      ))
    }
    
    fit_tbl <- dplyr::tibble(
      model  = names(fits),
      k      = vapply(fits, function(x) attr(stats::logLik(x), "df"), numeric(1)),
      logLik = vapply(fits, function(x) as.numeric(stats::logLik(x)), numeric(1)),
      AIC    = vapply(fits, stats::AIC, numeric(1)),
      BIC    = vapply(fits, stats::BIC, numeric(1))
    ) |>
      dplyr::arrange(.data$AIC)
    
    best_name <- fit_tbl$model[1]
    best_fit  <- fits[[best_name]]
    
    # build the matched null (drop interaction, keep main effects)
    if (best_name == "m1") {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_var))
    } else if (best_name == "m2") {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s", disease_var, time_fac_var))
    } else {
      rhs_null <- rhs_plus_covars(sprintf("%s + %s + I(%s^2)", disease_var, time_var, time_var))
    }
    f_null   <- stats::as.formula(sprintf("%s ~ %s", out_col, rhs_null))
    null_fit <- fit_safe(f_null, tmp)
    
    lrt <- if (!is.null(null_fit)) stats::anova(null_fit, best_fit) else NA
    lrt_p <- if (is.data.frame(lrt) && nrow(lrt) >= 2 && "p-value" %in% colnames(lrt)) lrt$`p-value`[2] else NA_real_
    
    # label best model
    best_model_desc <- dplyr::case_when(
      best_name == "m1" ~ "time as continuous",
      best_name == "m2" ~ "time as categorical",
      best_name == "m3" ~ "quadratic time",
      TRUE ~ NA_character_
    )
    
    tTab <- summary(best_fit)$tTable
    data.frame(
      NAMscore    = pc_name,
      best_model  = best_model_desc,
      coefficient = rownames(tTab),
      tTab,
      lrt_pval    = lrt_p,
      AIC         = fit_tbl$AIC[1],
      row.names   = NULL,
      check.names = FALSE
    )
  }
  
  results_list <- lapply(pc_names, fit_one_pc)
  results_all  <- do.call(rbind, results_list)
  
  object@pipeline_output <- results_all
  object
}