#' Generate pseudo principal components with time-dependent structure
#'
#' Simulates a matrix of pseudo principal components (PCs) by combining
#' Gaussian noise with structured effects derived from metadata variables
#' (e.g., cell type, disease status, visit, batch). Each PC can be assigned
#' to reflect one or more biological or technical factors with configurable
#' signal-to-noise ratios.
#'
#' This function is primarily intended for simulation and benchmarking of
#' longitudinal single-cell analysis workflows.
#'
#' @param data A data.frame containing per-cell metadata. Each row represents
#'   one cell/sample and must include the columns specified by the `*_col`
#'   arguments.
#' @param n_pcs Integer. Number of pseudo PCs to generate.
#'
#' @param cluster_pcs,disease_pcs,sex_pcs,age_pcs,bmi_pcs,batch_pcs,
#'   interaction_pcs,visit_pcs,subject_pcs Integer vectors indicating which
#'   PC indices are influenced by the corresponding metadata variable.
#'   Use `0` or an empty vector to disable an effect.
#'
#' @param scale_factor Numeric. Controls how variance decays across PCs.
#'   Higher values produce lower variance for higher-index PCs.
#'
#' @param cluster_ratio,disease_ratio,sex_ratio,age_ratio,bmi_ratio,
#'   batch_ratio,visit_ratio,interaction_ratio,subject_ratio Numeric values
#'   between 0 and 1 specifying the proportion of variance attributable to
#'   each factor for PCs listed in the corresponding `*_pcs` argument.
#'
#' @param cluster_col,disease_col,sex_col,age_col,bmi_col,batch_col,
#'   visit_col,interact_term_col Character scalars giving column names in
#'   `data` used to derive structured effects.
#'
#' @param seed Integer. Random seed for reproducibility of the simulated PCs.
#'   Note that the seed is applied implicitly through repeated calls to
#'   `rnorm()` and `sample()` inside the function.
#' @param interaction_pcs Integer vector of PC indices influenced by the
#'   interaction term (from `interact_term_col`). Use `0` or `integer(0)` to
#'   disable.
#' @param visit_pcs Integer vector of PC indices influenced by the visit/time
#'   variable (from `visit_col`). Use `0` or `integer(0)` to disable.
#' @param subject_pcs Integer vector of PC indices influenced by subject-level
#'   structure. (Currently reserved for future extension; kept for API
#'   completeness.)
#'
#' @param batch_ratio Numeric in [0, 1]. Proportion of variance attributable to
#'   batch effects for PCs listed in `batch_pcs`.
#' @param visit_ratio Numeric in [0, 1]. Proportion of variance attributable to
#'   visit/time effects for PCs listed in `visit_pcs`.
#' @param interaction_ratio Numeric in [0, 1]. Proportion of variance
#'   attributable to the interaction term for PCs listed in `interaction_pcs`.
#' @param subject_ratio Numeric in [0, 1]. Proportion of variance attributable
#'   to subject effects for PCs listed in `subject_pcs`.
#'
#' @param visit_col Character scalar. Column name in `data` encoding visit/time.
#' @param interact_term_col Character scalar. Column name in `data` encoding the
#'   interaction term used for `interaction_pcs` (e.g., a visit-by-disease term).
#'
#' @return A numeric matrix of dimension `nrow(data) x n_pcs`, where each column
#'   corresponds to a simulated pseudo principal component.
#'
#' @examples
#' \donttest{
#' ## Example using pre-loaded metadata:
#' ## meta <- readRDS("path/to/metadata.rds")
#' ## pcs <- generate_pseudo_pcs_time(
#' ##   data = meta,
#' ##   n_pcs = 10,
#' ##   cluster_pcs = 1:5,
#' ##   cluster_ratio = 0.3,
#' ##   visit_pcs = 6:10,
#' ##   visit_ratio = 0.2,
#' ##   seed = 123
#' ## )
#' ## dim(pcs)
#' }
#'
#' @export
generate_pseudo_pcs_time <- function(data,

                                     n_pcs = 20,

                                     cluster_pcs = 1:20,
                                     disease_pcs = 0,
                                     sex_pcs = 0,
                                     age_pcs = 0,
                                     bmi_pcs = 0,
                                     batch_pcs = 0,
                                     interaction_pcs = 0,
                                     visit_pcs = 0,
                                     subject_pcs = 0,

                                     scale_factor = 2,

                                     cluster_ratio = 0.25,
                                     disease_ratio = 0,
                                     sex_ratio = 0,
                                     age_ratio = 0,
                                     bmi_ratio = 0,
                                     batch_ratio = 0,
                                     visit_ratio = 0,
                                     interaction_ratio = 0,
                                     subject_ratio = 0,

                                     cluster_col = "cell_type",
                                     disease_col = "disease",
                                     sex_col = "sex",
                                     age_col = "age",
                                     bmi_col = "bmi",
                                     batch_col = "batch",
                                     visit_col = "visit",
                                     interact_term_col = "interact_term",
                                     seed = 1234
) {
  # set.seed(seed)
  disease_pcs = sample(disease_pcs)
  # set.seed(seed*2)
  sex_pcs = sample(sex_pcs)
  # set.seed(seed*3)
  age_pcs = sample(age_pcs)
  # set.seed(seed*4)
  bmi_pcs = sample(bmi_pcs)
  # set.seed(seed*5)
  batch_pcs = sample(batch_pcs)

  pcs_matrix <- sapply(1:n_pcs, function(x){

    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluster_col]))
    n_pcs <- n_pcs

    cell_clusters <- as.integer(factor(data[, cluster_col]))

    cell_clusters_list <- vector("list", 3)
    for (j in 1:3) {
      cell_clusters_tmp <- cell_clusters
      shufle <- sample(unique(cell_clusters_tmp)) * 10L
      for (i in seq_along(unique(cell_clusters))) {
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(shufle[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      cell_clusters_list[[j]] <- cell_clusters_tmp
    }

    cell_clusters1 <- cell_clusters_list[[1]]
    # (optional: you can also use [[2]] / [[3]] if you want later)

    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    cell_interaction <- as.integer(factor(data[,interact_term_col]))

    variance <- 1 / (x * scale_factor)
    # set.seed(seed*x)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))
    # set.seed(seed*x)
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters1), sd = sqrt(variance))
    # set.seed(seed*x)
    pc_disease = rnorm(n_cells, mean = scale(cell_diseases), sd = sqrt(variance))
    # set.seed(seed*x*2)
    pc_sex = rnorm(n_cells, mean = scale(cell_sex), sd = sqrt(variance))
    # set.seed(seed*x*3)
    pc_age = rnorm(n_cells, mean = scale(cell_age), sd = sqrt(variance))
    # set.seed(seed*x*4)
    pc_bmi = rnorm(n_cells, mean = scale(cell_bmi), sd = sqrt(variance))
    # set.seed(seed*x*5)
    pc_batch = rnorm(n_cells, mean = scale(cell_batch), sd = sqrt(variance))
    # set.seed(seed*x*6)
    pc_interact = rnorm(n_cells, mean = scale(cell_interaction), sd = sqrt(variance))
    # set.seed(seed*x*7)
    pc_visit = rnorm(n_cells, mean = scale(data[,visit_col]), sd = sqrt(variance))

    cluster_ratio_tmp     = cluster_ratio     / charmatch(x, cluster_pcs)
    disease_ratio_tmp     = disease_ratio     / charmatch(x, disease_pcs)
    sex_ratio_tmp         = sex_ratio         / charmatch(x, sex_pcs)
    age_ratio_tmp         = age_ratio         / charmatch(x, age_pcs)
    bmi_ratio_tmp         = bmi_ratio         / charmatch(x, bmi_pcs)
    batch_ratio_tmp       = batch_ratio       / charmatch(x, batch_pcs)
    interaction_ratio_tmp = interaction_ratio / charmatch(x, interaction_pcs)
    visit_ratio_tmp       = visit_ratio       / charmatch(x, visit_pcs)

    if (!(x %in% cluster_pcs))     cluster_ratio_tmp <- 0
    if (!(x %in% disease_pcs))     disease_ratio_tmp <- 0
    if (!(x %in% sex_pcs))         sex_ratio_tmp     <- 0
    if (!(x %in% age_pcs))         age_ratio_tmp     <- 0
    if (!(x %in% bmi_pcs))         bmi_ratio_tmp     <- 0
    if (!(x %in% batch_pcs))       batch_ratio_tmp   <- 0
    if (!(x %in% interaction_pcs)) interaction_ratio_tmp <- 0
    if (!(x %in% visit_pcs))       visit_ratio_tmp   <- 0

    pc * (1
          - cluster_ratio_tmp
          - disease_ratio_tmp
          - sex_ratio_tmp
          - age_ratio_tmp
          - bmi_ratio_tmp
          - batch_ratio_tmp
          - interaction_ratio_tmp
          - visit_ratio_tmp) +
      pc_cluster  * cluster_ratio_tmp +
      pc_disease  * disease_ratio_tmp +
      pc_sex      * sex_ratio_tmp +
      pc_age      * sex_ratio_tmp +
      pc_bmi      * sex_ratio_tmp +
      pc_batch    * sex_ratio_tmp +
      pc_interact * interaction_ratio_tmp +
      pc_visit    * visit_ratio_tmp
  })

  colnames(pcs_matrix) <- paste0("PC", seq_len(ncol(pcs_matrix)))
  return(pcs_matrix)
}
