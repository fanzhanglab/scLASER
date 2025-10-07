association_nam <- function(seurat_object=NULL,
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
