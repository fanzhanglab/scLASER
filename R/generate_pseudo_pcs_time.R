#' Title
#'
#' @param data
#' @param n_pcs
#' @param cluster_pcs
#' @param disease_pcs
#' @param sex_pcs
#' @param age_pcs
#' @param bmi_pcs
#' @param batch_pcs
#' @param interaction_pcs
#' @param visit_pcs
#' @param subject_pcs
#' @param scale_factor
#' @param cluster_ratio
#' @param disease_ratio
#' @param sex_ratio
#' @param age_ratio
#' @param bmi_ratio
#' @param batch_ratio
#' @param visit_ratio
#' @param interaction_ratio
#' @param subject_ratio
#' @param cluster_col
#' @param disease_col
#' @param sex_col
#' @param age_col
#' @param bmi_col
#' @param batch_col
#' @param visit_col
#' @param interact_term_col
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
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

    cell_clusters <- as.integer(factor(data[,cluster_col]))
    for (j in 1:3){
      cell_clusters_tmp = cell_clusters
      # set.seed(x*j)
      shufle = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(shufle[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters",j),cell_clusters_tmp)
    }
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
