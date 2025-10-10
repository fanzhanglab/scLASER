


suppressWarnings({
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(ggalluvial)
  library(ggrepel)
  library(MASS)
  library(caret)

  library(Seurat)
  library(ggplot2)
  library(glue)

  library(stevemisc)
  library(stevedata)
  library(lme4)
  library(broom.mixed)


  library(doParallel)
  library(pbapply)

  library(variancePartition)
  library(pheatmap)
  library(Seurat)

  library(data.table)
  library(presto)
  library(edgeR)
  library(dplyr)
  library(MatrixEQTL)
  library(harmony)

  library(MASS)
  library(foreach)
  library(doParallel)

  library(MOFA2)
  library(shapviz)
  library(xgboost)
})



generate_dummy_metadata <- function(n_cells = 3000,
                                              sd_celltypes = 0.1,
                                              n_major_cell_types = 7,
                                              n_minor_cell_types = 3,
                                              relative_abundance = 0.1,
                                              n_major_diff_celltypes = 1,
                                              n_minor_diff_celltypes = 1,
                                              n_individuals = 30,
                                              n_batchs = 4,
                                              prop_sex = 0.5,
                                              prop_disease = 0.5,
                                              fc_increase = 0.1,
                                              seed = 1234
) {
  n_cell_types = n_major_cell_types + n_minor_cell_types


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



  major_cell_types <- n_major_cell_types

  rare_cell_types <- n_minor_cell_types

  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "subject_id")


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











  diff_clusters = c(1:n_major_diff_celltypes, (n_cell_types-n_minor_diff_celltypes+1):n_cell_types)
  diff_cell_types = LETTERS[diff_clusters]

  for (i in LETTERS[1:n_cell_types]) {
    if(i %in% diff_cell_types){

      abundance = dplyr::left_join(celltype_df,
                                   dummy_data,
                                   by="subject_id") %>%
        dplyr::filter(cell_type == i) %>%
        dplyr::group_by(subject_id) %>%
        dplyr::summarise(pro = dplyr::n()/n_cells,
                         count = dplyr::n())
      diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_increase)) - n_cells*median(abundance$pro))

      i

      len = length(unique(dummy_data[dummy_data$disease == "1",]$subject_id)) * diff

      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "1",]$subject_id),diff)),
        celltype_df
      )
    }
  }
  dummy_data = merge(dummy_data,celltype_df,by="subject_id")

  set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,diff_cell_types))
}


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

  library(doParallel)
  library(foreach)
  library(dplyr)

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


  cl <- makeCluster(n_thread)
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
    }








    variance <- 1 / cell_clusters_means2
    set.seed(seed*x)







    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)


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
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*2)
    pc_disease = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))





    var_all = c(var_all, variance)

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
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*3)
    pc_sex = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)


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
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*4)
    pc_age = rnorm(n_cells, mean = scale(data[,age_col]), sd = sqrt(variance))






    var_all = c(var_all, variance)

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
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*5)
    pc_bmi = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)

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
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*6)
    pc_batch = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    data.frame(cell_covariate,pc_batch) %>%
      dplyr::group_by(cell_covariate) %>%
      dplyr::summarize(mean_pc = mean(pc_batch),
                       sd_pc = sd(pc_batch))

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


FUN <- function(x, n, frac) {
  min_byFrac = as.integer(length(x)*frac)
  if (as.integer(length(x)) <= n) {
    return(x)
  } else if (as.integer(length(x)) > n & min_byFrac <= n) {
    return(x[x %in% sample(x, n)])
  } else if (as.integer(length(x)) > n & min_byFrac > n) {
    return(x[x %in% sample(x, min_byFrac)])
  }
}
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}




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


create_NAM = function(metadata,
                      pcs,
                      samplem_key = NULL,
                      graph_use = 'RNA_snn',
                      batches = NULL,
                      covs = NULL ,
                      nsteps = NULL,
                      verbose=TRUE ,
                      assay=NULL,
                      key='NAMPC_'
){

  meta = metadata
  rownames(meta) = 1:nrow(meta)

  m <- as(t(pcs), "dgTMatrix")
  colnames(m) = 1:ncol(m)

  obj <- Seurat::CreateSeuratObject(
    counts = m,
    meta.data = meta,
    assay = 'RNA',
    names.field = 1
  )

  harmony_embeddings_all = pcs
  rownames(harmony_embeddings_all) = 1:nrow(harmony_embeddings_all)

  obj@reductions$harmony = Seurat::CreateDimReducObject(
    embeddings = harmony_embeddings_all,
    stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
    assay = "RNA",
    key = Seurat::Key("HARMONY", quiet = TRUE)
  )

  obj <- obj %>%
    Seurat::FindNeighbors(verbose = TRUE, reduction = 'harmony', dims = 1:20, k.param = 30, nn.eps = 0)

  seurat_object = obj


  covs_keep <- NULL
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)

  if (length(names(seurat_object@graphs)) == 0) {
    stop('Must precompute graph in Seurat with FindNeighbors()')
  }
  if (is.null(graph_use)) {
    graph_use <- names(seurat_object@graphs)[[1]]
    message('Graph not specified. Using graph {graph_use}')
  } else {
    if (!graph_use %in% names(seurat_object@graphs)) {
      stop('Graph {graph_use} not in seurat object')
    }
  }
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, one_of(covs_keep))))
  obs_df <- tibble::rownames_to_column(seurat_object@meta.data, 'CellID')
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop(
      'Sample-level metadata is same length as cell-level metadata.
             Please check that samplem_vars are sample-level covariates.'
    )
  }

  rcna_data <- list(
    samplem = samplem_df,
    obs = obs_df,
    connectivities = seurat_object@graphs[[graph_use]],
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )
  data = rcna_data
  suffix=''



  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(as.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }
  max_frac_pcs=0.15
  res <- list()
  covs_mat <- as.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))

  f <- as.formula(as.character(glue('~0+{data$samplem_key}')))
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]
  s <- s[, data$samplem[[data$samplem_key]]]






  prevmedkurt <- Inf
  maxnsteps=15L









  diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(as.matrix(res))
  }





  for (i in seq_len(maxnsteps)) {
    print(i)
    s <- diffuse_step(data, s)
    print(s)
    medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis))
    if (is.null(nsteps)) {
      prevmedkurt <- medkurt
      print(prevmedkurt)
      print(medkurt)
      if (prevmedkurt - medkurt < 3 & i > 3) {
        message(glue::glue('stopping after {i} steps'))
        break
      }
    } else if (i == nsteps) {
      break
    }
  }








  snorm <- t(prop.table(s, 2))
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]

  NAM=snorm

  N <- nrow(NAM)

  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))

  } else {
    message('filtering based on batches kurtosis')
  }



  .batch_kurtosis <- function(NAM, batches_vec) {
    purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
      if (length(i)>1) {
        Matrix::colMeans(NAM[i, ])
      } else if (length(i)==1) {
        Matrix::colMeans(t(NAM[i, ]))
      }
    }) %>%
      dplyr::bind_cols() %>%
      apply(1, moments::kurtosis)
  }



  kurtoses <- .batch_kurtosis(NAM, batches_vec)

  threshold <- max(6, 2*median(kurtoses))
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}'))
  keep <- which(kurtoses < threshold)






  .res_qc_nam = list(NAM = NAM, keep = keep)

  res <- list()
  res[[paste0('NAM.T', suffix)]] <- t(.res_qc_nam[[1]])
  res[[paste0('keptcells', suffix)]] <- .res_qc_nam[[2]]
  res[[paste0('_batches', suffix)]] <- batches_vec



  if (verbose) message('Residualize NAM')
  N <- nrow(NAM)
  NAM_ <- scale(NAM, center = TRUE, scale = FALSE)
  ncols_C <- 0
  if (!is.null(covs_mat)) {
    covs_mat <- scale(covs_mat)
    ncols_C <- ncols_C + ncol(covs_mat)
  }
  if (is.null(covs_mat)) {
    M <- Matrix::Diagonal(n = N)
  } else {
    M <- Matrix::Diagonal(n = N) - covs_mat %*% solve(t(covs_mat) %*% covs_mat, t(covs_mat))
  }





  NAM_ <- M %*% NAM_












  .res_resid_nam = list(


    NAM_=scale(NAM_, center=FALSE, scale=TRUE),
    M=M,
    r=ncols_C
  )
  res[[paste0('_M', suffix)]] <- .res_resid_nam$M
  res[[paste0('_r', suffix)]] <- .res_resid_nam$r

  if (verbose) message('Decompose NAM')
  npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))

  npcs <- min(npcs, nrow(data$samplem) - 1)
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(scale(NAM_, center=FALSE, scale=TRUE))
  } else {
    svd_res <- RSpectra::svds(scale(NAM_, center=FALSE, scale=TRUE), k = npcs)
  }


  dim(svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v))



  U_df <- svd_res$u[, seq_len(npcs)]
  colnames(U_df) <- paste0('PC', seq_len(npcs))
  rownames(U_df) <- rownames(NAM_)
  V_df <- svd_res$v[, seq_len(npcs)]
  colnames(V_df) <- paste0('PC', seq_len(npcs))
  rownames(V_df) <- colnames(NAM_)
  .res_svd_nam <- list(U=U_df, svs=svd_res$d^2, V=V_df)

  res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
  res[[paste0('NAM_svs', suffix)]] <- .res_svd_nam$svs
  res[[paste0('NAM_varexp', suffix)]] <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
  res[[paste0('NAM_nbhdXpc', suffix)]] <- .res_svd_nam$V

  nam_res = res






  NAMsvd=list(
    nam_res$NAM_sampleXpc,
    nam_res$NAM_svs,
    nam_res$NAM_nbhdXpc,
    nam_res$NAM_varexp
  )
  names(NAMsvd) = c("sampleXpc","svs","nbhdXpc","varexp")
  return(NAMsvd)
}


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

empirical_fdrs <- function(z, znull, thresholds) {
  n <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:n, ])
  ranks <- t(tail_counts(thresholds, z)[1:n, ])


  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)

  return(fdr)
}

tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
  })
}

association_nam = function(seurat_object=NULL,
                           metadata=NULL,
                           pcs=NULL,
                           test_var=NULL,
                           samplem_key = NULL,
                           graph_use = 'RNA_snn',
                           batches = NULL,
                           covs = NULL ,
                           nsteps = NULL,
                           verbose=TRUE ,
                           assay=NULL,
                           key='NAMPC_',
                           maxnsteps=15L,
                           max_frac_pcs=0.15,
                           ks=NULL,
                           Nnull=1000,
                           force_permute_all=FALSE,
                           local_test=TRUE,
                           seed=1234,
                           return_nam=TRUE
){
  if (!is.null(seurat_object) && is.null(metadata) && is.null(pcs)) {
    message(paste0("will use Seurat object following analysis..."))
  } else if (is.null(seurat_object) && !is.null(metadata) && !is.null(pcs)){

    meta = metadata
    rownames(meta) = 1:nrow(meta)

    m <- as(t(pcs), "dgTMatrix")
    colnames(m) = 1:ncol(m)

    obj <- Seurat::CreateSeuratObject(
      counts = m,
      meta.data = meta,
      assay = 'RNA',
      names.field = 1
    )

    harmony_embeddings_all = pcs
    rownames(harmony_embeddings_all) = 1:nrow(harmony_embeddings_all)

    obj@reductions$harmony = Seurat::CreateDimReducObject(
      embeddings = harmony_embeddings_all,
      stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
      assay = "RNA",
      key = Seurat::Key("HARMONY", quiet = TRUE)
    )

    obj <- obj %>%
      Seurat::FindNeighbors(verbose = TRUE, reduction = 'harmony', dims = 1:20, k.param = 30, nn.eps = 0)

    seurat_object = obj
  } else if ((is.null(seurat_object) && is.null(metadata) && !is.null(pcs))|
             (is.null(seurat_object) && !is.null(metadata) && is.null(pcs))){
    stop('Must provide both metadata and precomuputed PCs')
  }


  covs_keep <- test_var
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)

  if (length(names(seurat_object@graphs)) == 0) {
    stop('Must precompute graph in Seurat with FindNeighbors()')
  }
  if (is.null(graph_use)) {
    graph_use <- names(seurat_object@graphs)[[1]]
    message('Graph not specified. Using graph {graph_use}')
  } else {
    if (!graph_use %in% names(seurat_object@graphs)) {
      stop('Graph {graph_use} not in seurat object')
    }
  }
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, one_of(covs_keep))))
  obs_df <- tibble::rownames_to_column(seurat_object@meta.data, 'CellID')
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop(
      'Sample-level metadata is same length as cell-level metadata.
             Please check that samplem_vars are sample-level covariates.'
    )
  }

  rcna_data <- list(
    samplem = samplem_df,
    obs = obs_df,
    connectivities = seurat_object@graphs[[graph_use]],
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )
  data = rcna_data
  suffix=''



  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(data.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }

  res <- list()
  covs_mat <- data.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))

  f <- as.formula(as.character(glue('~0+{data$samplem_key}')))
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]
  s <- s[, data$samplem[[data$samplem_key]]]






  prevmedkurt <- Inf









  diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    return(as.matrix(res))
  }





  for (i in seq_len(maxnsteps)) {
    s <- diffuse_step(data, s)
    medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis))
    if (is.null(nsteps)) {
      prevmedkurt <- medkurt
      if (prevmedkurt - medkurt < 3 & i > 3) {
        message(glue::glue('stopping after {i} steps'))
        break
      }
    } else if (i == nsteps) {
      break
    }
  }








  snorm <- t(prop.table(s, 2))
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]

  NAM=snorm

  N <- nrow(NAM)

  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))
  } else {
    message('filtering based on batches kurtosis')
  }



  .batch_kurtosis <- function(NAM, batches_vec) {
    purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
      if (length(i)>1) {
        Matrix::colMeans(NAM[i, ])
      } else if (length(i)==1) {
        Matrix::colMeans(t(NAM[i, ]))
      }
    }) %>%
      dplyr::bind_cols() %>%
      apply(1, moments::kurtosis)
  }



  kurtoses <- .batch_kurtosis(NAM, batches_vec)

  threshold <- max(6, 2*median(kurtoses))
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}'))
  keep <- which(kurtoses < threshold)






  .res_qc_nam = list(NAM = NAM, keep = keep)

  res <- list()
  res[[paste0('NAM.T', suffix)]] <- t(.res_qc_nam[[1]])
  res[[paste0('keptcells', suffix)]] <- .res_qc_nam[[2]]
  res[[paste0('_batches', suffix)]] <- batches_vec



  if (verbose) message('Residualize NAM')
  N <- nrow(NAM)
  NAM_ <- scale(NAM, center = TRUE, scale = FALSE)
  ncols_C <- 0
  if (!is.null(covs_mat)) {
    covs_mat <- scale(covs_mat)
    ncols_C <- ncols_C + ncol(covs_mat)
  }
  if (is.null(covs_mat)) {
    M <- Matrix::Diagonal(n = N)
  } else {
    M <- Matrix::Diagonal(n = N) - covs_mat %*% solve(t(covs_mat) %*% covs_mat, t(covs_mat))
  }





  NAM_ <- M %*% NAM_












  .res_resid_nam = list(


    NAM_=scale(NAM_, center=FALSE, scale=TRUE),
    M=M,
    r=ncols_C
  )
  res[[paste0('_M', suffix)]] <- .res_resid_nam$M
  res[[paste0('_r', suffix)]] <- .res_resid_nam$r

  if (verbose) message('Decompose NAM')
  npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))

  npcs <- min(npcs, nrow(data$samplem) - 1)
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(scale(NAM_, center=FALSE, scale=TRUE))
  } else {
    svd_res <- RSpectra::svds(scale(NAM_, center=FALSE, scale=TRUE), k = npcs)
  }


  dim(svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v))



  U_df <- svd_res$u[, seq_len(npcs)]
  colnames(U_df) <- paste0('PC', seq_len(npcs))
  rownames(U_df) <- rownames(NAM_)
  V_df <- svd_res$v[, seq_len(npcs)]
  colnames(V_df) <- paste0('PC', seq_len(npcs))
  rownames(V_df) <- colnames(NAM_)
  .res_svd_nam <- list(U=U_df, svs=svd_res$d^2, V=V_df)

  res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
  res[[paste0('NAM_svs', suffix)]] <- .res_svd_nam$svs
  res[[paste0('NAM_varexp', suffix)]] <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
  res[[paste0('NAM_nbhdXpc', suffix)]] <- .res_svd_nam$V






  nam_res = res






  NAMsvd=list(
    nam_res$NAM_sampleXpc,
    nam_res$NAM_svs,
    nam_res$NAM_nbhdXpc,
    nam_res$NAM_varexp
  )

  names(NAMsvd) = c("sampleXpc","svs","nbhdXpc","varexp")

  M=res[[paste0('_M', suffix)]]
  r=res[[paste0('_r', suffix)]]

  yvals <- rcna_data$samplem[[test_var]]
  if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer') ) {
    stop('test_var is of class {class(yvals)}. It must be numeric variable for association testing.')
  }
  y = yvals

  if (is.null(seed)) {
    set.seed(sample(1e6, 1))
  }
  if (force_permute_all) {
    batches_vec <- rep(1L, length(y))
  }


  U <- NAMsvd[[1]]
  sv <- NAMsvd[[2]]
  V <- NAMsvd[[3]]
  y <- scale(y)
  n <- length(y)

  if (is.null(ks)) {
    incr <- max(round(0.02*n), 1)
    maxnpcs <- min(4*incr, round(n/5))
    ks <- seq(incr, maxnpcs+1, incr)
  }


  .reg <- function(q, k) {
    Xpc <- U[, 1:k]
    beta <- t(Xpc) %*% q
    qhat <- Xpc %*% beta
    return(list(qhat = qhat, beta = beta))
  }

  .stats <- function(yhat, ycond, k) {
    ssefull <- crossprod(yhat - ycond)
    ssered <- crossprod(ycond)
    deltasse <-  ssered - ssefull
    f <- (deltasse / k) / (ssefull/n)
    p <- -pf(f, k, n-(1+r+k), log.p = TRUE)
    r2 <- 1 - ssefull/ssered
    return(list(p=p, r2=r2))
  }

  .minp_stats <- function(z) {
    zcond <- scale(M %*% z, center = FALSE, scale = TRUE)
    qhats <- purrr::map(ks, function(k) .reg(zcond, k)$qhat)
    .tmp <- purrr::map2(qhats, ks, function(qhat, k) .stats(qhat, zcond, k))
    ps <- purrr::map_dbl(.tmp, 'p')
    r2s <- purrr::map_dbl(.tmp, 'r2')
    k_ <- which.min(ps)
    return(list(k=ks[k_], p=ps[k_], r2=r2s[k_]))
  }



  .tmp <- .minp_stats(y)
  k <- .tmp$k
  p <- .tmp$p
  r2 <- .tmp$r2
  if (k == max(ks)) {
    warning(glue::glue('data supported use of {k} NAM PCs, which is the maximum considered. Consider allowing more PCs by using the "ks" argument.'))
  }

  ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
  .tmp <- .reg(ycond, k)
  yhat <- .tmp$qhat
  beta <- .tmp$beta
  r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))**2


  ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta/n)
  rownames(ncorrs) <- rownames(V)

  set.seed(seed)
  y_ <- conditional_permutation(batches_vec, y, Nnull)
  .tmp <- apply(y_, 2, .minp_stats)
  nullminps <- purrr::map_dbl(.tmp, 'p')
  nullr2s <- purrr::map_dbl(.tmp, 'r2')

  pfinal <- (sum(nullminps <= p+1e-8) + 1)/(Nnull + 1)
  if (sum(nullminps <= p+1e-8) == 0) {
    warning('global association p-value attained minimal possible value. Consider increasing Nnull')
  }


  fdrs <- NULL
  fdr_5p_t <- NULL
  fdr_10p_t <- NULL
  fdr_20p_t <- NULL
  fdr_30p_t <- NULL
  fdr_40p_t <- NULL
  fdr_50p_t <- NULL

  if (local_test) {
    message('computing neighborhood-level FDRs')
    Nnull <- min(1000, Nnull)
    y_ <- y_[, 1:Nnull]
    ycond_ <- scale(M %*% y_, center = FALSE, scale = TRUE)
    gamma_ <- crossprod(U[, 1:k], ycond_)
    nullncorrs <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_ / n)))

    maxcorr <- max(abs(ncorrs))
    fdr_thresholds <- seq(maxcorr/4, maxcorr, maxcorr/400)
    fdr_vals <- empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
    fdrs = data.frame(

      threshold = head(fdr_thresholds, -1),
      fdr = fdr_vals,
      num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t))
    )

    if (min(fdrs$fdr) > 0.05) {
      fdr_5p_t <- NULL
    } else {
      fdr_5p_t <- min(subset(fdrs, fdr < 0.05)$threshold)
    }
    if (min(fdrs$fdr) > 0.05) {
      fdr_10p_t <- NULL
    } else {
      fdr_10p_t <- min(subset(fdrs, fdr < 0.1)$threshold)
    }
    if (min(fdrs$fdr) > 0.05) {
      fdr_20p_t <- NULL
    } else {
      fdr_20p_t <- min(subset(fdrs, fdr < 0.2)$threshold)
    }
    if (min(fdrs$fdr) > 0.05) {
      fdr_30p_t <- NULL
    } else {
      fdr_30p_t <- min(subset(fdrs, fdr < 0.3)$threshold)
    }
    if (min(fdrs$fdr) > 0.05) {
      fdr_40p_t <- NULL
    } else {
      fdr_40p_t <- min(subset(fdrs, fdr < 0.4)$threshold)
    }
    if (min(fdrs$fdr) > 0.05) {
      fdr_50p_t <- NULL
    } else {
      fdr_50p_t <- min(subset(fdrs, fdr < 0.5)$threshold)
    }
  }

  res <- list(
    p = pfinal,
    nullminps=nullminps,
    k=k,
    ncorrs=ncorrs,
    fdrs=fdrs,
    fdr_5p_t=fdr_5p_t,
    fdr_10p_t=fdr_10p_t,
    yhat=yhat,
    ycond=ycond,
    ks=ks,
    beta=beta,
    r2=r2,
    r2_perpc=r2_perpc,
    nullr2_mean=mean(nullr2s),
    nullr2_std=sd(nullr2s)
  )

  if (return_nam) {
    res[['NAM_embeddings']] <- nam_res$NAM_nbhdXpc
    res[['NAM_loadings']] <- nam_res$NAM_sampleXpc
    res[['NAM_svs']] <- nam_res$NAM_svs
    res[['NAM_varexp']] <- nam_res$NAM_varexp
  }

  seurat_object[['cna']] <- Seurat::CreateDimReducObject(
    embeddings = res$NAM_embeddings,
    loadings = res$NAM_loadings,
    stdev = res$NAM_svs,
    assay = assay,
    key = key,
    misc = res
  )

  seurat_object@meta.data$cna_ncorrs <- ncorrs[colnames(seurat_object), , drop=TRUE]

  seurat_object@meta.data$cna_ncorrs_fdr05 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_5p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_5p_t)
    seurat_object@meta.data$cna_ncorrs_fdr05[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }

  seurat_object@meta.data$cna_ncorrs_fdr10 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_10p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_10p_t)
    seurat_object@meta.data$cna_ncorrs_fdr10[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }

  seurat_object@meta.data$cna_ncorrs_fdr20 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_20p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_20p_t)
    seurat_object@meta.data$cna_ncorrs_fdr20[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }

  seurat_object@meta.data$cna_ncorrs_fdr30 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_30p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_30p_t)
    seurat_object@meta.data$cna_ncorrs_fdr30[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }

  seurat_object@meta.data$cna_ncorrs_fdr40 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_40p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_40p_t)
    seurat_object@meta.data$cna_ncorrs_fdr40[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }

  seurat_object@meta.data$cna_ncorrs_fdr50 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_50p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$cna_ncorrs) >= fdr_50p_t)
    seurat_object@meta.data$cna_ncorrs_fdr50[idx_passed] <- seurat_object@meta.data$cna_ncorrs[idx_passed]
  }
  return(seurat_object)
}

NAM_NMF_GP = function(seurat_object=NULL,
                      metadata=NULL,
                      pcs=NULL,
                      test_var=NULL,
                      interaction_feature=NULL,
                      samplem_key = NULL,
                      graph_use = 'RNA_snn',
                      batches = NULL,
                      covs = NULL ,
                      nsteps = NULL,
                      verbose=TRUE ,
                      assay=NULL,
                      key='NAMPC_',
                      maxnsteps=15L,
                      max_frac_pcs=0.3,
                      model_p=TRUE,
                      Nnull=500,
                      force_permute_all=FALSE,
                      local_test=TRUE,
                      seed=1234,
                      return_nam=TRUE
){





  if (!is.null(seurat_object) && is.null(metadata) && is.null(pcs)) {
    message(paste0("will use Seurat object following analysis..."))
  } else if (is.null(seurat_object) && !is.null(metadata) && !is.null(pcs)){
    message(paste0("Computing nearest neighbor graph..."))
    meta = metadata
    rownames(meta) = 1:nrow(meta)

    m <- as(t(pcs), "dgTMatrix")
    colnames(m) = 1:ncol(m)

    obj <- Seurat::CreateSeuratObject(
      counts = m,
      meta.data = meta,
      assay = 'RNA',
      names.field = 1
    )

    harmony_embeddings_all = pcs
    rownames(harmony_embeddings_all) = 1:nrow(harmony_embeddings_all)

    obj@reductions$harmony = Seurat::CreateDimReducObject(
      embeddings = harmony_embeddings_all,
      stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
      assay = "RNA",
      key = Seurat::Key("HARMONY", quiet = TRUE)
    )

    obj <- obj %>%
      Seurat::FindNeighbors(verbose = FALSE, reduction = 'harmony', dims = 1:20, k.param = 30, nn.eps = 0)

    seurat_object = obj
  } else if ((is.null(seurat_object) && is.null(metadata) && !is.null(pcs))|
             (is.null(seurat_object) && !is.null(metadata) && is.null(pcs))){
    stop('Must provide both metadata and precomuputed PCs')
  }


  covs_keep <- test_var
  covs_keep <- c(covs_keep, interaction_feature)
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)

  if (length(names(seurat_object@graphs)) == 0) {
    stop('Must precompute graph in Seurat with FindNeighbors()')
  }
  if (is.null(graph_use)) {
    graph_use <- names(seurat_object@graphs)[[1]]
    message('Graph not specified. Using graph {graph_use}')
  } else {
    if (!graph_use %in% names(seurat_object@graphs)) {
      stop('Graph {graph_use} not in seurat object')
    }
  }
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, one_of(covs_keep))))
  obs_df <- tibble::rownames_to_column(seurat_object@meta.data, 'CellID')
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop(
      'Sample-level metadata is same length as cell-level metadata.
             Please check that samplem_vars are sample-level covariates.'
    )
  }

  rcna_data <- list(
    samplem = samplem_df,
    obs = obs_df,
    connectivities = seurat_object@graphs[[graph_use]],
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )
  data = rcna_data
  data$interaction_feature_1 = interaction_feature_1



  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(data.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }

  f <- as.formula(as.character(glue('~0+{data$samplem_key}')))
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]
  s <- s[, data$samplem[[data$samplem_key]]]

  diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees

    res <- (a %*% s_norm) + s_norm
    return(as.matrix(res))
  }

  prevmedkurt <- Inf

  for (i in seq_len(maxnsteps)) {
    s <- diffuse_step(data, s)
    medkurt <- median(apply(prop.table(s, 2), 1, moments::kurtosis))
    if (is.null(nsteps)) {
      prevmedkurt <- medkurt
      if (prevmedkurt - medkurt < 3 & i > 3) {
        message(glue::glue('Diffusion process stoped after {i} steps'))
        break
      }
    } else if (i == nsteps) {
      break
    }
  }

  NAM <- t(prop.table(s, 2))
  rownames(NAM) <- data$samplem[[data$samplem_key]]
  colnames(NAM) <- data$obs[[data$obs_key]]









  scaling = function (x, range = NULL)
  {
    (x - min(c(x, range), na.rm = TRUE))/(max(c(x, range), na.rm = TRUE) -
                                            min(c(x, range), na.rm = TRUE))
  }
  NAM_scaled = apply(NAM, 2, scaling)



  npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))
  npcs <- min(npcs, nrow(data$samplem) - 1)
  k=npcs

  set.seed(seed)
  out_NMF <- nnTensor::NMF(NAM_scaled,

                           J=k)

  U <- out_NMF$U[, seq_len(k)]
  colnames(U) <- paste0('Comp', seq_len(k))
  rownames(U) <- rownames(NAM_scaled)
  V <- out_NMF$V[, seq_len(k)]
  colnames(V) <- paste0('Comp', seq_len(k))
  rownames(V) <- colnames(NAM_scaled)







  yvals <- rcna_data$samplem[[test_var]]
  if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer') ) {
    stop('test_var is of class {class(yvals)}. It must be numeric variable for association testing.')
  }
  y = yvals

  y_cell <- dplyr::left_join(data.frame(subject_id = rcna_data$obs[,c(samplem_key)]),
                             data.frame(y = y,
                                        subject_id = rcna_data$samplem[[samplem_key]]),
                             by="subject_id") %>%
    .$y
  n <- length(y)

  if (is.factor(samplem_df[,interaction_feature])){
    form_df = data.frame(y = y,
                         U,
                         inter_var = samplem_df[,interaction_feature],
                         samplem_df[,c(covs,batch_col,samplem_key)]
    )
  } else if (is.numeric(samplem_df[,interaction_feature])){
    inter_bin = cut(rcna_data$samplem[,interaction_feature], breaks=ceiling(length(unique(rcna_data$samplem[,interaction_feature]))/2))
    form_df = data.frame(y = y,
                         U,
                         inter_var = inter_bin ,
                         samplem_df[,c(covs,batch_col,samplem_key)]
    )
  }

  for(i in 1:length(covs)){
    if (is.numeric(samplem_df[,covs[i]])) {
      form_df[,1+ncol(U)+1+i] = as.numeric(data.matrix(form_df[,covs[i]]))
    } else if (is.factor(samplem_df[,covs[i]])) {
      form_df[,1+ncol(U)+1+i] = as.factor(data.matrix(form_df[,covs[i]]))
    }
  }
  form_df[,batch_col] = as.factor(form_df[,batch_col])
  form_df[,samplem_key] = as.factor(form_df[,samplem_key])

  frml = paste0("y ~ ",
                paste0("gp(Comp", 1:k, ")", collapse = " + "),
                " + ",
                paste0('gp(Comp', 1:k, ")", "*zs(inter_var)", collapse = " + "),
                " + ",
                "zs(", samplem_key, ") + zs(", batch_col, ")")

  for(i in 1:length(covs)){
    if (is.numeric(samplem_df[,covs[i]])) {
      covs_frml = paste0(
        " + ",
        paste0("gp(",covs[i],")", collapse = " + "))
    } else if (is.factor(samplem_df[,covs[i]])) {
      covs_frml = paste0(
        " + ",
        paste0("zs(",covs[i],")", collapse = " + "))
    }
  }
  message(paste0("Formula of AGP model: gp(NAM-Comp1:", k, ") + gp(NAM-Comp1:", k, ")*zs(",
                 if(is.numeric(samplem_df[,interaction_feature])) {
                   paste0(interaction_feature,"[bin]")
                 } else if (is.factor(samplem_df[,interaction_feature])) {
                   interaction_feature
                 },
                 ")", covs_frml, " + zs(", samplem_key, ") + zs(", batch_col, ")"))
  frml = paste0(frml, covs_frml)
  frml = as.formula(frml)

  prior <- list(
    alpha = lgpr::normal(mu = 0, sigma = 1),
    ell = lgpr::igam(shape = 5, scale = 5),
    wrp = lgpr::log_normal(0, 1)






  )

  message("Fitting AGP model...")
  set.seed(seed)
  fit <- lgp(frml,
             data = form_df,
             prior = prior,
             chains = 4,
             cores = thread,
             control = list(adapt_delta = 0.95),

             seed = seed,
             refresh  = 0,
             quiet = TRUE,
             show_messages = FALSE)



  message("Calculating estimated Y using fitted AGP model...")
  p <- pred(fit, form_df,
            reduce = mean,
            verbose = FALSE)


  ssefull <- crossprod(as.numeric(p@y_mean) - y)
  ssefull = as.numeric(ssefull)








  sf <- fit@stan_fit
  sampler_params <- rstan::get_sampler_params(sf, inc_warmup = FALSE)







  if (is.factor(samplem_df[,interaction_feature])){
    x_pred <- data.frame(y = y_cell,
                         V[,1:k],
                         inter_var = rcna_data$obs[,interaction_feature],
                         rcna_data$obs[,c(covs,batch_col,samplem_key)],
                         CellID = rcna_data$obs[["CellID"]])
  } else if (is.numeric(samplem_df[,interaction_feature])){
    inter_bin_cell = cut(rcna_data$obs[,interaction_feature], breaks=ceiling(length(unique(rcna_data$obs[,interaction_feature]))/2))
    x_pred <- data.frame(y = y_cell,
                         V[,1:k],
                         inter_var = inter_bin_cell,
                         rcna_data$obs[,c(covs,batch_col,samplem_key)],
                         CellID = rcna_data$obs[["CellID"]])
  }

  for(i in 1:length(covs)){
    if (is.numeric(samplem_df[,covs[i]])) {
      x_pred[,1+ncol(V[,1:k])+1+i] = as.numeric(data.matrix(x_pred[,covs])[,i])
    } else if (is.factor(samplem_df[,covs[i]])) {
      x_pred[,1+ncol(V[,1:k])+1+i] = as.factor(data.matrix(x_pred[,covs])[,i])
    }
  }
  x_pred[,batch_col] = as.factor(x_pred[,batch_col])
  x_pred[,samplem_key] = as.factor(x_pred[,samplem_key])


  if(model_p){
    message("Calculating empirical p-value by comparing full and null model...")


    huber_loss <- function(y_true, y_pred, delta = 1.0) {

      error <- y_true - y_pred


      abs_error <- abs(error)
      quadratic <- pmin(abs_error, delta)
      linear <- abs_error - quadratic
      loss <- 0.5 * quadratic^2 + delta * linear

      return(mean(loss))
    }

    frml_null = paste0("y ~ ",
                       paste0("gp(Comp", 1:k, ")", collapse = " + "),
                       " + ",


                       "zs(", samplem_key, ") + zs(", batch_col, ")")
    frml_null = paste0(frml_null, covs_frml)
    frml_null = as.formula(frml_null)
    message(paste0("Formula of null model: gp(NAM-Comp1:", k, ") + zs(",
                   if(is.numeric(samplem_df[,interaction_feature])) {
                     paste0(interaction_feature,"[bin]")
                   } else if (is.factor(samplem_df[,interaction_feature])) {
                     interaction_feature
                   },
                   ")", covs_frml, " + zs(", samplem_key, ") + zs(", batch_col, ")"))


    .perm_stats = function(i){
      require(lgpr)

      set.seed(seed+i)
      fit_full <- lgp(frml,
                      data = form_df,
                      prior = prior,
                      chains = 4,
                      cores = thread,
                      control = list(adapt_delta = 0.95),
                      iter = 100,
                      seed = seed+i,
                      refresh  = 0,
                      quiet = TRUE,
                      show_messages = FALSE)
      p_cell_full <- pred(fit_full,
                          x_pred,
                          reduce = mean,
                          verbose = FALSE,
                          force = TRUE)
      hl_full = huber_loss(as.numeric(p_cell_full@y_mean), x_pred$y)

      set.seed(seed+i)
      fit_null <- lgp(frml_null,
                      data = form_df,
                      prior = prior,
                      chains = 4,
                      cores = thread,
                      control = list(adapt_delta = 0.95),
                      iter = 100,
                      seed = seed+i,
                      refresh  = 0,
                      quiet = TRUE,
                      show_messages = FALSE)
      p_cell_null <- pred(fit_null,
                          x_pred,
                          reduce = mean,
                          verbose = FALSE,
                          force = TRUE)
      hl_null = huber_loss(as.numeric(p_cell_null@y_mean), x_pred$y)

      return(list(hl_full = hl_full, hl_null = hl_null))
    }


    cl = parallel::makeCluster(thread)
    parallel::clusterExport(cl,
                            c("x_pred","frml","frml_null","prior","form_df","huber_loss","seed","thread"),
                            envir=environment()
    )
    perm_list <- pbapply::pblapply(1:Nnull, .perm_stats, cl = cl)
    stopCluster(cl = cl)
    fullhls = unlist(lapply(perm_list, function(x) x$hl_full))
    nullhls = unlist(lapply(perm_list, function(x) x$hl_null))
    pfinal <- (sum(nullhls <= fullhls) + 1)/(Nnull + 1)
    if (sum(nullhls <= fullhls) == 0) {
      warning('global association p-value attained minimal possible value. Consider increasing Nnull')
    }
    message(paste0("model p-value = ", pfinal))
  }

  message("Projecting AGP model to cell neiggborhood...")
  p_cell <- pred(fit,
                 x_pred,
                 reduce = mean,
                 verbose = FALSE,
                 force = TRUE)

  scale_by_group <- function(value, group) {
    scaled_value = data.frame(value = value,
                              group = group) %>%
      group_by(group) %>%
      mutate(scaled_value = scale(value


      )) %>%
      .$scaled_value %>%
      as.numeric()
    return(scaled_value)
  }



  if(is.factor(rcna_data$obs[,interaction_feature])){
    ncorrs = scale_by_group(value = as.numeric(p_cell@y_mean),
                            group = rcna_data$obs[,interaction_feature])
  } else if (is.numeric(rcna_data$obs[,interaction_feature])) {

    ncorrs = scale_by_group(value = as.numeric(p_cell@y_mean),
                            group = inter_bin_cell)
  }
  seurat_object@meta.data$NMFGP_ncorrs <- ncorrs


  conditional_permutation <- function(B, Y, num) {
    purrr::map(seq_len(num), function(i) {
      split(seq_len(length(Y)), B) %>%
        purrr::map(function(idx) {
          data.frame(idx, val=sample(Y[idx], replace = TRUE))
        }) %>% dplyr::bind_rows() %>%
        dplyr::arrange(idx) %>%
        with(val)
    }) %>%
      purrr::reduce(Matrix::cbind2)
  }





  empirical_fdrs <- function(z, znull, thresholds) {
    N <- length(thresholds) - 1
    tails <- t(tail_counts(thresholds, znull)[1:N, ])
    ranks <- t(tail_counts(thresholds, z)[1:N, ])


    fdp <- sweep(tails, 2, ranks, '/')
    fdr <- Matrix::colMeans(fdp)

    return(fdr)
  }

  tail_counts <- function(z, znull) {
    apply(znull, 2, function(znulli) {
      as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))
    })
  }

  fdrs <- NULL
  fdr_5p_t <- NULL
  fdr_10p_t <- NULL
  fdr_20p_t <- NULL
  fdr_30p_t <- NULL
  fdr_40p_t <- NULL
  fdr_50p_t <- NULL

  if (local_test) {
    message('Computing neighborhood-level FDRs')
    set.seed(seed)
    y_ <- conditional_permutation(batches_vec, y, Nnull)
    Nnull <- min(1000, Nnull)
    y_ <- y_[, 1:Nnull]
    gamma_ <- crossprod(U[, 1:k], y_)
    m_ = (n-k-length(c(interaction_feature,covs,samplem_key,batch_col)))
    if (m_ > 0) {
      nullncorrs <- abs(V[, 1:k] %*% (gamma_/m_))
    } else {
      nullncorrs <- abs(V[, 1:k] %*% (gamma_/n))
    }
    nullncorrs = abs(scale(nullncorrs))














    maxcorr <- max(abs(ncorrs))
    fdr_thresholds <- seq(maxcorr/4, maxcorr, maxcorr/400)
    fdr_vals <- empirical_fdrs(matrix(ncorrs,ncol=1), nullncorrs, fdr_thresholds)
    fdrs = data.frame(

      threshold = head(fdr_thresholds, -1),
      fdr = fdr_vals,
      num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t))
    )

    if (min(fdrs$fdr) > 0.05) {
      fdr_5p_t <- NULL
    } else {
      fdr_5p_t <- min(subset(fdrs, fdr < 0.05)$threshold)
    }
    if (min(fdrs$fdr) > 0.10) {
      fdr_10p_t <- NULL
    } else {
      fdr_10p_t <- min(subset(fdrs, fdr < 0.1)$threshold)
    }
    if (min(fdrs$fdr) > 0.20) {
      fdr_20p_t <- NULL
    } else {
      fdr_20p_t <- min(subset(fdrs, fdr < 0.2)$threshold)
    }
    if (min(fdrs$fdr) > 0.30) {
      fdr_30p_t <- NULL
    } else {
      fdr_30p_t <- min(subset(fdrs, fdr < 0.3)$threshold)
    }
    if (min(fdrs$fdr) > 0.40) {
      fdr_40p_t <- NULL
    } else {
      fdr_40p_t <- min(subset(fdrs, fdr < 0.4)$threshold)
    }
    if (min(fdrs$fdr) > 0.50) {
      fdr_50p_t <- NULL
    } else {
      fdr_50p_t <- min(subset(fdrs, fdr < 0.5)$threshold)
    }
  }


  seurat_object@meta.data$NMFGP_ncorrs_fdr05 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_5p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_5p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr05[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  seurat_object@meta.data$NMFGP_ncorrs_fdr10 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_10p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_10p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr10[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  seurat_object@meta.data$NMFGP_ncorrs_fdr20 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_20p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_20p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr20[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  seurat_object@meta.data$NMFGP_ncorrs_fdr30 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_30p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_30p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr30[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  seurat_object@meta.data$NMFGP_ncorrs_fdr40 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_40p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_40p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr40[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  seurat_object@meta.data$NMFGP_ncorrs_fdr50 <- rep(0, nrow(seurat_object@meta.data))
  if (!is.null(fdr_50p_t)) {
    idx_passed <- which(abs(seurat_object@meta.data$NMFGP_ncorrs) >= fdr_50p_t)
    seurat_object@meta.data$NMFGP_ncorrs_fdr50[idx_passed] <- seurat_object@meta.data$NMFGP_ncorrs[idx_passed]
  }

  res <- list(prediction = p_cell,
              model = fit,
              AGPM_params = sampler_params,
              p = pfinal,




              k=k,
              ncorrs=ncorrs,
              fdrs=fdrs,
              fdr_5p_t=fdr_5p_t,
              fdr_10p_t=fdr_10p_t
  )

  U <- out_NMF$U[, seq_len(k)]
  colnames(U) <- paste0('NAM_Comp', seq_len(k))
  rownames(U) <- rownames(NAM_scaled)
  V <- out_NMF$V[, seq_len(k)]
  colnames(V) <- paste0('NAM_Comp', seq_len(k))
  rownames(V) <- colnames(NAM_scaled)
  seurat_object[['NMFGP']] <- Seurat::CreateDimReducObject(
    embeddings = V,
    loadings = U,
    assay = assay,
    key = key,
    misc = res
  )

  return(seurat_object)
}



NormalizeDataSeurat <- function(A, scaling_factor = 1e4, do_ftt = FALSE) {
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  A@x <- scaling_factor * A@x
  if (do_ftt) {
    A@x <- sqrt(A@x) + sqrt(1 + A@x)
  } else {
    A@x <- log(1 + A@x)
  }
  return(A)
}
cosine_normalize = function(x, dim){ apply(x, MARGIN = dim, function(a) a / sqrt(sum(a^2))) }


unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


apply_model_to_column <- function(column_data) {

  form_test <- as.formula(paste0("column_data ~ (1|cell_type) + (1|sex) + (1|disease) + (1|batch) + (1|subject_id) + age + bmi"))


  fit <- lme4::lmer(form_test, data = dummy_data, REML = FALSE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


  var_stats <- variancePartition::calcVarPart(fit)

  return(var_stats)
}

apply_countmodel_to_column <- function(column_data) {

  form_test <- as.formula(paste0("column_data ~ (1|cell_type) + (1|sex) + (1|disease) + (1|batch) + (1|subject_id) + age + bmi"))


  fit <- lme4::glmer(form_test, family = "poisson", nAGQ = 0, data = dummy_data, control = glmerControl(optimizer = "nloptwrap"))


  var_stats <- variancePartition::calcVarPart(fit)

  return(var_stats)
}



compute_lisi_parallel <- function(
    X, meta_data, label_colnames, perplexity = 30, nn_eps = 0, n_thread = 1
) {
  library(lisi)
  library(parallel)
  library(RANN)
  N <- nrow(meta_data)
  dknn <- nn2(X, k = perplexity * 3, eps = nn_eps)
  lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
  lisi_df <- Reduce(cbind, mclapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message('Cannot compute LISI on missing values')
      return(rep(NA, N))
    } else {

      dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
      dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
      labels <- as.integer(factor(labels)) - 1
      n_batches <- length(unique(labels))
      simpson <- compute_simpson_index(
        t(dknn$nn.dists),
        t(dknn$nn.idx) - 1,
        labels,
        n_batches,
        perplexity
      )
      return(1 / simpson)
    }
  }, mc.cores = n_thread))
  lisi_df <- as.data.frame(lisi_df)
  colnames(lisi_df) <- label_colnames
  row.names(lisi_df) <- row.names(meta_data)
  return(lisi_df)
}


