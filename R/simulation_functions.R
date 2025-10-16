# Required libraries

# remotes::install_version("Matrix", version = "1.6-3")
# suppressWarnings({
#   library(dplyr)
#   library(tidyr)
#  # library(tidyverse)
#   library(MASS) # Used for sampling from normal distribution
#   library(caret)
#
#   library(Seurat)
#   library(glue)
#
#   library(stevemisc)
#   library(stevedata)
#   library(lme4)
#   library(broom.mixed)
#
#   # With parallelism
#   library(doParallel)
#   library(foreach)
#   library(pbapply)
#
#   library(variancePartition)
#
#   library(Seurat)
# })

# Functions
## generate meta data

#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
generate_dummy_data_woInteraction <- function(n_cells = 3000, # cells of major cell types per individual
                                              sd_celltypes = 0.1,  # standard deviation of number of cells
                                              n_major_cell_types = 7,
                                              n_minor_cell_types = 3,
                                              relative_abundance = 0.1, # ratio between major and rare
                                              n_major_diff_celltypes = 1,
                                              n_minor_diff_celltypes = 1,
                                              n_individuals = 30, # total individuals
                                              n_batchs = 4,
                                              prop_sex = 0.5,
                                              prop_disease = 0.5,
                                              fc_interact = 0.1, # additional proportion of interacted cells which are from people with interacted group compared to non-interacted cell types
                                              seed = 1234
) {
  n_cell_types = n_major_cell_types + n_minor_cell_types

  # Generate unique subject IDs
  subject_id <- paste0("SUB_", 1:n_individuals)
  set.seed(seed)
  sex = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_sex)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_sex)))
  )
  if (length(sex) != length(subject_id)) {
    sex = c(sex, rep(1,length(subject_id)-length(sex)))
  }
  set.seed(seed*2)
  disease = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_disease)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_disease))
    )
  )
  if (length(disease) != length(subject_id)) {
    disease = c(disease, rep(1,length(subject_id)-length(disease)))
  }
  set.seed(seed*3)
  age <- sample(18:60, n_individuals, replace = TRUE)
  set.seed(seed*4)
  bmi <- sample(15:35, n_individuals, replace = TRUE)
  batch <- rep(rep(1:n_batchs), length(subject_id))[1:length(subject_id)]
  dummy_data = data.frame(subject_id = subject_id,
                          sex = factor(sex, levels = c(0, 1)),
                          disease = factor(disease, levels = c(0, 1)),
                          age = age,
                          batch = factor(batch, levels = c(0:max(batch))),
                          bmi = bmi
  )

  # Major and rare cell type counts
  ## major_cell_types <- ceiling(n_cell_types / 2)
  major_cell_types <- n_major_cell_types
  ## rare_cell_types <- n_cell_types - major_cell_types
  rare_cell_types <- n_minor_cell_types

  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "subject_id")

  # Generate baseline of cell type data.frame
  for (id in dummy_data$subject_id){
    set.seed(seed*5*grep(id,dummy_data$subject_id)*10)
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
    set.seed(seed*6*grep(id,dummy_data$subject_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
      celltype_df = rbind(celltype_df,
                          data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], n),
                                     subject_id = id)
      )
    }
  }

  # check
  ## celltype_df %>%
  ##   dplyr::group_by(subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  ## celltype_df %>%
  ##   dplyr::group_by(cell_type,subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()

  diff_clusters = c(1:n_major_diff_celltypes, (n_cell_types-n_minor_diff_celltypes+1):n_cell_types)
  diff_cell_types = LETTERS[diff_clusters]
  prop_control_types <- 0.1

  for (i in LETTERS[1:n_cell_types]) {
    abundance = dplyr::left_join(celltype_df,
                                 dummy_data,
                                 by="subject_id") %>%
      dplyr::filter(cell_type == i) %>%
      dplyr::group_by(subject_id) %>%
      dplyr::summarise(pro = dplyr::n()/n_cells,
                       count = dplyr::n())
    if(i %in% diff_cell_types){

      # abundance = dplyr::left_join(celltype_df,
      #                              dummy_data,
      #                              by="subject_id") %>%
      #   dplyr::filter(cell_type == i) %>%
      #   dplyr::group_by(subject_id) %>%
      #   dplyr::summarise(pro = dplyr::n()/n_cells,
      #                    count = dplyr::n())

      # prop(ortion) = cells for each subject that are this cell type
      print(head(abundance))
      # diff is the number of cells that the differential cell types should have (use fc_interact to scale/increase)
      diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) - n_cells*median(abundance$pro))
      print(paste("Inside diff cell types. diff = ", diff))
      print(paste("n_cells = ", n_cells, "; median pro = ", median(abundance$pro), " ; median pro increase = ", (median(abundance$pro)*(1+fc_interact))))
      print(paste("i = ", i))
      i# add diff cells only to disease == 1

      len = length(unique(dummy_data[dummy_data$disease == "1",]$subject_id)) * diff
      #print(paste("unique(dummy_data[dummy_data$disease == '1',]$subject_id) : ", unique(dummy_data[dummy_data$disease == "1",]$subject_id)))
      # unique... is the number of subjects with disease == 1?
      print(paste("len = ", len))

      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "1",]$subject_id),diff)),
        celltype_df
      )
    }
    else{
      # number of disease cells in control cells
      #diff_control = ceiling(n_cells*(median(abundance$pro)*(prop_control_types))+ n_cells*median(abundance$pro))
      diff_control = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) + n_cells*median(abundance$pro))
      print(paste("Inside control cell types. diff_control = ", diff_control))
      print(paste("n_cells = ", n_cells, "; median pro = ", median(abundance$pro), " ; median pro control = ", (median(abundance$pro)*(prop_control_types))))
      len = length(unique(dummy_data[dummy_data$disease == "0",]$subject_id)) * diff_control
      print(paste("len = ", len))
      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "0",]$subject_id),diff_control)),
        celltype_df
      )

    }
  }
  dummy_data = merge(dummy_data,celltype_df,by="subject_id")
  # initial attempt to increase number of 0 (control) cells in non disease cell types
  # Shuffle rows
  set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,diff_cell_types))
}


#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
generate_pseudo_pcs_woInteraction = function(data,
                                             # number of principal components
                                             n_pcs = 20,
                                             # Among the principal components representing each attribute, arrange the components in descending order of variance.
                                             cluster_pcs = 1:20,
                                             disease_pcs = 0,
                                             sex_pcs = 0,
                                             age_pcs = 0,
                                             bmi_pcs = 0,
                                             batch_pcs = 0,

                                             scale_factor = 2,
                                             # Define the maximum percentage of each attribute
                                             cluster_ratio = 0.25,
                                             disease_ratio = 0,
                                             sex_ratio = 0,
                                             age_ratio = 0,
                                             bmi_ratio = 0,
                                             batch_ratio = 0,

                                             cluser_col = "cell_type",
                                             disease_col = "disease",
                                             sex_col = "sex",
                                             age_col = "age",
                                             bmi_col = "bmi",
                                             batch_col = "batch",
                                             seed = 1234
) {
  set.seed(seed)
  disease_pcs = sample(disease_pcs)
  set.seed(seed*2)
  sex_pcs = sample(sex_pcs)
  set.seed(seed*3)
  age_pcs = sample(age_pcs)
  set.seed(seed*4)
  bmi_pcs = sample(bmi_pcs)
  set.seed(seed*5)
  batch_pcs = sample(batch_pcs)

  sapply(1:n_pcs, function(x){

    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluser_col]))
    n_pcs <- n_pcs

    cell_clusters <- as.integer(factor(data[,cluser_col]))
    for (j in 1:3){
      cell_clusters_tmp = cell_clusters
      set.seed(x*j)
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

    variance <- 1 / (x * scale_factor)  # Set so that the variance decreases as the principal component goes from 1 to 20.
    set.seed(seed*x)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))
    set.seed(seed*x)
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters1), sd = sqrt(variance))
    # for (j in 1:3){
    #   if (j==1){
    #     pc_cluster = rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance))
    #   } else {
    #     pc_cluster = pc_cluster + rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance))
    #   }
    # }
    set.seed(seed*x)
    pc_disease = rnorm(n_cells, mean = scale(cell_diseases), sd = sqrt(variance))
    set.seed(seed*x*2)
    pc_sex = rnorm(n_cells, mean = scale(cell_sex), sd = sqrt(variance))
    set.seed(seed*x*3)
    pc_age = rnorm(n_cells, mean = scale(cell_age), sd = sqrt(variance))
    set.seed(seed*x*4)
    pc_bmi = rnorm(n_cells, mean = scale(cell_bmi), sd = sqrt(variance))
    set.seed(seed*x*5)
    pc_batch = rnorm(n_cells, mean = scale(cell_batch), sd = sqrt(variance))
    # cluster_ratio_tmp = cluster_ratio / sqrt(charmatch(x,cluster_pcs))
    cluster_ratio_tmp = cluster_ratio / charmatch(x,cluster_pcs)
    # disease_ratio_tmp = disease_ratio / sqrt(charmatch(x,disease_pcs))
    disease_ratio_tmp = disease_ratio / charmatch(x,disease_pcs)
    # sex_ratio_tmp = sex_ratio / sqrt(charmatch(x,sex_pcs))
    sex_ratio_tmp = sex_ratio / charmatch(x,sex_pcs)
    # age_ratio_tmp = age_ratio / sqrt(charmatch(x,age_pcs))
    age_ratio_tmp = age_ratio / charmatch(x,age_pcs)
    # bmi_ratio_tmp = bmi_ratio / sqrt(charmatch(x,bmi_pcs))
    bmi_ratio_tmp = bmi_ratio / charmatch(x,bmi_pcs)
    # batch_ratio_tmp = batch_ratio / sqrt(charmatch(x,batch_pcs))
    batch_ratio_tmp = batch_ratio / charmatch(x,batch_pcs)
    if (!(x %in% cluster_pcs)) cluster_ratio_tmp = 0
    if (!(x %in% disease_pcs)) disease_ratio_tmp = 0
    if (!(x %in% sex_pcs)) sex_ratio_tmp = 0
    if (!(x %in% age_pcs)) age_ratio_tmp = 0
    if (!(x %in% bmi_pcs)) bmi_ratio_tmp = 0
    if (!(x %in% batch_pcs)) batch_ratio_tmp = 0
    return(pc * (1
                 -cluster_ratio_tmp
                 -disease_ratio_tmp
                 -sex_ratio_tmp
                 -age_ratio_tmp
                 -bmi_ratio_tmp
                 -batch_ratio_tmp
    ) +
      pc_cluster * cluster_ratio_tmp +
      pc_disease * disease_ratio_tmp +
      pc_sex * sex_ratio_tmp +
      pc_age * sex_ratio_tmp +
      pc_bmi * sex_ratio_tmp +
      pc_batch * sex_ratio_tmp
    )
  })
}

#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
generate_pseudo_features = function(data,
                                    n_features = 20,
                                    cluster_features = 1:20,
                                    disease_features = 0,
                                    sex_features = 0,
                                    age_features = 0,
                                    bmi_features = 0,
                                    batch_features = 0,
                                    individual_features = 0,

                                    cluster_ratio = 0.25,
                                    disease_ratio = 0,
                                    sex_ratio = 0,
                                    age_ratio = 0,
                                    bmi_ratio = 0,
                                    batch_ratio = 0,
                                    individual_ratio = 0.1,
                                    ratio_variance = 0.5,

                                    cluster_col = "cell_type",
                                    disease_col = "disease",
                                    sex_col = "sex",
                                    age_col = "age",
                                    bmi_col = "bmi",
                                    batch_col = "batch",
                                    individual_col = "batch",
                                    seed = 1234,
                                    n_thread = 1,
                                    verbose = TRUE
) {

  # library(doParallel)
  # library(foreach)
  # library(dplyr)

  set.seed(seed)
  disease_features = sample(disease_features)
  set.seed(seed*2)
  sex_features = sample(sex_features)
  set.seed(seed*3)
  age_features = sample(age_features)
  set.seed(seed*4)
  bmi_features = sample(bmi_features)
  set.seed(seed*5)
  batch_features = sample(batch_features)
  set.seed(seed*6)
  individual_features = sample(individual_features)

  # Register parallel backend to use multiple cores
  cl <- makeCluster(n_thread) # Leave one core free for other processes
  registerDoParallel(cl)

  pcs_list = foreach(x = 1:n_features, .packages = c("dplyr")) %dopar% {

    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluster_col]))

    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    cell_individual <- as.integer(factor(data[,individual_col]))

    var_all <- c()

    # Generate dummy features reflecting cell clusters
    ## CELL CLUSTERS
    cell_clusters <- as.integer(factor(data[,cluster_col]))
    for (j in 1:2){
      cell_clusters_tmp = cell_clusters
      set.seed(x*j)
      cluster_mean = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters_means",j),cell_clusters_tmp)
      print(paste("length of cell_clusters_tmp: ", length(cell_clusters_tmp)))
    }
    # paste0("cell_clusters_means",j) include cluster-specific mean values
    ## length(cell_clusters_means1)==length(cell_clusters)
    ## length(cell_clusters_means2)==length(cell_clusters)
    # check if cluster means are same value by cluster
    ## cell_clusters_means1[cell_clusters==1]
    ## cell_clusters_means2[cell_clusters==1]
    ## cell_clusters_means1[cell_clusters==2]
    ## cell_clusters_means2[cell_clusters==2]
    variance <- 1 / cell_clusters_means2 # cell type specific variance
    set.seed(seed*x)
    ### =====
    # The following R code generates dummy PC values by cells from a normal distribution by means of random numbers.
    # Specifically, the rnorm function is used to generate n_cells of random numbers, and the mean of each PC value for each cell cluster is the (scaled) unique value (cell_clusters_means1) for each cell type. The square root of the specific variance is the standard deviation.
    #1. rnorm(n) generates n random values from a standard normal distribution (mean 0, standard deviation 1).
    #2. mean = scale(cell_clusters1) specifies the mean of the random values to be generated. Here, the scale() function standardizes the data (subtracts the mean and divides by the standard deviation), but the use of the scale function may be inappropriate in this context. Normally, the scale() function standardizes a column of vectors or data frames and converts them to a value with mean 0 and standard deviation 1. However, if the goal here is to set the mean, then the use of mean(cell_clusters1) is appropriate.
    #3. sd = sqrt(variance) specifies the standard deviation of the random value to be generated. variance is the value of its variance, and sqrt(variance) gives the standard deviation.
    ### =====
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)

    ## DISEASE
    # Similary, generate dummy PC for other covariates with changing seeds
    cell_covariate = cell_diseases
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*2)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*2)
    pc_disease = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    # check if PC values are different by disease group
    ## data.frame(cell_covariate,pc_disease) %>%
    ##   dplyr::group_by(cell_covariate) %>%
    ##   dplyr::summarize(mean_pc = mean(pc_disease),
    ##                    sd_pc = sd(pc_disease))
    var_all = c(var_all, variance)

    ## SEX
    cell_covariate = cell_sex
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*3)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*3)
    pc_sex = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)

    # Treat age as fixed effect for PC mean
    ## AGE
    cell_covariate = cell_age
    for (j in 2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*4)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*4)
    pc_age = rnorm(n_cells, mean = scale(data[,age_col]), sd = sqrt(variance))
    ## data.frame(data[,age_col],pc_age) %>%
    ##   magrittr::set_colnames(c("cell_covariate","pc_age")) %>%
    ##   dplyr::group_by(cell_covariate) %>%
    ##     dplyr::summarize(mean_pc = mean(pc_age),
    ##                      sd_pc = sd(pc_age)) %>%
    ##   as.data.frame()
    var_all = c(var_all, variance)
    ## BMI
    cell_covariate = cell_bmi
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*5)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*5)
    pc_bmi = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    ## BATCH
    cell_covariate = cell_batch
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*6)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*6)
    pc_batch = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    data.frame(cell_covariate,pc_batch) %>%
      dplyr::group_by(cell_covariate) %>%
      dplyr::summarize(mean_pc = mean(pc_batch),
                       sd_pc = sd(pc_batch))
    # SUBJECT ID?
    cell_covariate = cell_individual
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*7)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*7)
    pc_individual = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    data.frame(cell_covariate,pc_individual) %>%
      dplyr::group_by(cell_covariate) %>%
      dplyr::summarize(mean_pc = mean(pc_individual),
                       sd_pc = sd(pc_individual))

    # Generate dummy PC regardless cell types or other covariates (noise term)
    set.seed(seed*x*100)
    variance = sample(var_all,n_cells)
    set.seed(seed*x*100)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))

    cluster_ratio_tmp = cluster_ratio + runif(1, min = -cluster_ratio*ratio_variance, max = cluster_ratio*ratio_variance)
    disease_ratio_tmp = disease_ratio + runif(1, min = -disease_ratio*ratio_variance, max = disease_ratio*ratio_variance)
    sex_ratio_tmp = sex_ratio + runif(1, min = -sex_ratio*ratio_variance, max = sex_ratio*ratio_variance)
    age_ratio_tmp = age_ratio + runif(1, min = -age_ratio*ratio_variance, max = age_ratio*ratio_variance)
    bmi_ratio_tmp = bmi_ratio + runif(1, min = -bmi_ratio*ratio_variance, max = bmi_ratio*ratio_variance)
    batch_ratio_tmp = batch_ratio  + runif(1, min = -batch_ratio*ratio_variance, max = batch_ratio*ratio_variance)
    individual_ratio_tmp = individual_ratio  + runif(1, min = -individual_ratio*ratio_variance, max = individual_ratio*ratio_variance)
    if (!(x %in% cluster_features)) cluster_ratio_tmp = 0
    if (!(x %in% disease_features)) disease_ratio_tmp = 0
    if (!(x %in% sex_features)) sex_ratio_tmp = 0
    if (!(x %in% age_features)) age_ratio_tmp = 0
    if (!(x %in% bmi_features)) bmi_ratio_tmp = 0
    if (!(x %in% batch_features)) batch_ratio_tmp = 0
    if (!(x %in% individual_features)) individual_ratio_tmp = 0

    noise_ratio_tmp = (1
                       -cluster_ratio_tmp
                       -disease_ratio_tmp
                       -sex_ratio_tmp
                       -age_ratio_tmp
                       -bmi_ratio_tmp
                       -batch_ratio_tmp
                       -individual_ratio_tmp
    )
    if (noise_ratio_tmp < 0) noise_ratio_tmp = 0

    if(verbose){
      message(paste0("Feature",x,";\ncluster ratio = ",cluster_ratio_tmp,
                     "\ndisease ratio = ", disease_ratio_tmp,
                     "\nsex ratio = ", sex_ratio_tmp,
                     "\nage ratio = ", age_ratio_tmp,
                     "\nbmi ratio = ", bmi_ratio_tmp,
                     "\nbatch ratio = ", batch_ratio_tmp,
                     "\nindividual ratio = ", individual_ratio_tmp,
                     "\nnoise ratio = ", noise_ratio_tmp))
    }
    return(scale(pc) * noise_ratio_tmp +
             scale(pc_cluster) * cluster_ratio_tmp  +
             scale(pc_disease) * disease_ratio_tmp  +
             scale(pc_sex) * sex_ratio_tmp  +
             scale(pc_age) * age_ratio_tmp  +
             scale(pc_bmi) * bmi_ratio_tmp  +
             scale(pc_batch) * batch_ratio_tmp +
             scale(pc_individual) * individual_ratio_tmp
    )
  }

  stopCluster(cl)

  pcs <- do.call(cbind, pcs_list)
  return(pcs)
}

#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
conditional_permutation <- function(B, Y, num) {
  purrr::map(seq_len(num), function(i) {
    split(seq_len(length(Y)), B) %>% purrr::map(function(idx) {
      data.frame(idx, val=sample(Y[idx]))
    }) %>% dplyr::bind_rows() %>%
      dplyr::arrange(idx) %>%
      with(val)
  }) %>%
    purrr::reduce(Matrix::cbind2)
}

#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
empirical_fdrs <- function(z, znull, thresholds) {
  n <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:n, ])
  ranks <- t(tail_counts(thresholds, z)[1:n, ])

  # compute FDPs
  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)

  return(fdr)
}

#' Title (optional, can be brief)
#' @noRd
#' @keywords internal
tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
  })
}



