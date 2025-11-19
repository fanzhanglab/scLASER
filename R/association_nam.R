#' Association testing between NAM embeddings and a sample-level variable
#'
#' Computes Neighborhood Aggregation Matrix (NAM) embeddings from a Seurat graph,
#' residualizes them with respect to optional covariates, performs an SVD to obtain
#' NAM PCs, and tests association between the NAM subspace and a numeric
#' sample-level variable via a min-\emph{p} procedure with conditional permutations.
#' Optionally writes the results back into the Seurat object as a reduction (`cna`)
#' and per-cell neighborhood correlations (with FDR-thresholded variants).
#'
#' You can either (a) pass a Seurat object that already has a nearest-neighbor graph
#' (e.g. created by `Seurat::FindNeighbors()`), or (b) pass `metadata` and precomputed
#' low-dimensional coordinates (`pcs`), in which case a minimal Seurat object and graph
#' are built internally.
#'
#' @param seurat_object A `Seurat` object. If provided, it must contain at least one
#'   cell-cell graph in `object@graphs` (e.g. after `FindNeighbors()`).
#'   If `NULL`, you must provide both `metadata` and `pcs`.
#' @param metadata A data.frame of sample-level metadata (rows = samples, columns = covariates).
#'   Required only when `seurat_object` is `NULL`. Must include `test_var`, any `batches`/`covs`,
#'   and the `samplem_key` column that links samples to cells.
#' @param pcs A numeric matrix of precomputed embeddings (rows = samples, cols = PCs)
#'   used to build a graph when `seurat_object` is `NULL`. Ignored otherwise.
#' @param test_var Character. Name of the numeric sample-level variable to test
#'   (must be a column of `metadata` or `seurat_object@meta.data` after aggregation).
#' @param samplem_key Character. Name of the sample identifier column used to map
#'   cells to samples (present in cell-level `meta.data`) and to index `metadata`.
#' @param graph_use Character. Name of the graph in `seurat_object@graphs` to use
#'   (default `'RNA_snn'`). If `NULL`, the first available graph is used.
#' @param batches Optional character vector of sample-level batch columns used only
#'   for conditional permutation strata and batch-kurtosis QC.
#' @param covs Optional character vector of sample-level covariate columns to residualize out
#'   prior to association testing.
#' @param nsteps Integer or `NULL`. Number of graph diffusion steps to form NAM.
#'   If `NULL`, an adaptive stopping rule based on median kurtosis is used.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @param assay Character or `NULL`. Assay name for the created reduction. Default `NULL`.
#' @param key Character. Reduction key prefix for the created reduction. Default `'NAMPC_'`.
#' @param maxnsteps Integer. Maximum diffusion steps when using adaptive stopping. Default `15L`.
#' @param max_frac_pcs Numeric in (0,1]. Maximum fraction of samples to use as the
#'   upper bound on SVD PCs (also bounded by sample size). Default `0.15`.
#' @param ks Optional integer vector of candidate numbers of NAM PCs to consider
#'   in the min-\emph{p} procedure. If `NULL`, a data-driven sequence is used.
#' @param Nnull Integer. Number of conditional permutations for the global test. Default `1000`.
#' @param force_permute_all Logical. If `TRUE`, ignore `batches` and permute across all samples. Default `FALSE`.
#' @param local_test Logical. If `TRUE`, compute neighborhood-level empirical FDR curves and thresholds. Default `TRUE`.
#' @param seed Integer or `NULL`. RNG seed used for permutations. Default `1234`.
#' @param return_nam Logical. If `TRUE`, store NAM embeddings/loadings/singular values/variance explained
#'   in the returned reduction. Default `TRUE`.
#'
#' @return A `Seurat` object where:
#' \itemize{
#'   \item `reductions$cna` contains:
#'     \itemize{
#'       \item `embeddings`: neighborhood-by-PC NAM embeddings (nbhd x PC)
#'       \item `loadings`: sample-by-PC NAM loadings
#'       \item `stdev`: singular values
#'       \item `misc`: a list with fields: `p` (global p-value), `nullminps`, `k` (selected PCs),
#'         `ncorrs` (per-cell neighborhood correlations), `fdrs` (data.frame of FDR vs threshold),
#'         `fdr_5p_t`, `fdr_10p_t`, `yhat`, `ycond`, `ks`, `beta`, `r2`,
#'         `r2_perpc`, `nullr2_mean`, `nullr2_std`, projection matrix `M` and rank `r`.
#'     }
#'   \item `meta.data$cna_ncorrs`: per-cell (neighborhood) correlation score.
#'   \item `meta.data$cna_ncorrs_fdr{05,10,20,30,40,50}`: thresholded versions at target FDRs.
#' }
#'
#' @details
#' **Pipeline summary:**\
#' (1) Build NAM by diffusing sample indicators over the cell graph;\
#' (2) Batch-kurtosis QC to drop unstable neighborhoods; \
#' (3) Residualize NAM against covariates; \
#' (4) SVD to obtain NAM PCs; \
#' (5) Choose \eqn{k} via a min-\emph{p} scan over candidate `ks`; \
#' (6) Obtain a global p-value using conditional permutations (stratified by `batches` if provided); \
#' (7) Optionally compute neighborhood-level empirical FDR thresholds.
#'
#' @section Requirements:
#' Relies on \pkg{Seurat}, \pkg{Matrix}, \pkg{RSpectra}, \pkg{moments}, \pkg{dplyr},
#' \pkg{tibble}, \pkg{purrr}, and \pkg{glue}. Ensure a neighbor graph exists if supplying
#' a `Seurat` object (see `Seurat::FindNeighbors()`).
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}}, \code{\link[Seurat]{CreateSeuratObject}},
#'   \code{\link[Seurat]{CreateDimReducObject}}
#'
#' @examples
#' \dontrun{
#' # Case A: start from an existing Seurat object with a graph
#' obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:20)
#' obj <- association_nam(
#'   seurat_object = obj,
#'   test_var = "numeric_sample_trait",
#'   samplem_key = "sample_id",
#'   covs = c("age","sex"),
#'   batches = "batch",
#'   graph_use = "RNA_snn",
#'   local_test = TRUE
#' )
#'
#' # Case B: build from metadata + precomputed PCs
#' obj2 <- association_nam(
#'   seurat_object = NULL,
#'   metadata = meta_df,             # rows = samples
#'   pcs = pcs_mat,                  # rows = samples, cols = PCs
#'   test_var = "numeric_sample_trait",
#'   samplem_key = "sample_id",
#'   local_test = FALSE
#' )
#' }
#'
#' @export

association_nam_LV <- function(object = NULL,
                               seurat_object = NULL,
                               metadata = NULL,
                               pcs = NULL,
                               test_var = NULL,
                               samplem_key = NULL,
                               graph_use = 'RNA_snn',
                               batches = NULL,
                               covs = NULL,
                               nsteps = NULL,
                               verbose = TRUE,
                               assay = NULL,
                               key = 'NAMPC_',
                               maxnsteps = 15L,
                               max_frac_pcs = 0.15,   # kept for backward compatibility
                               n_pcs = NULL,          # NEW: explicitly request number of PCs
                               ks = NULL,
                               Nnull = 1000,
                               force_permute_all = FALSE,
                               local_test = TRUE,
                               seed = 1234,
                               return_nam = TRUE) {

  ## ---- NEW: pull metadata/pcs from scLASER object if provided ----
  if (!is.null(object)) {
    metadata <- object@metadata
    pcs      <- object@harmony
  }

  if (!is.null(seurat_object) && is.null(metadata) && is.null(pcs)) {
    message("will use Seurat object following analysis...")
  } else if (is.null(seurat_object) && !is.null(metadata) && !is.null(pcs)) {
    # build a minimal Seurat object from metadata + precomputed PCs (stored as 'harmony' reduction)
    meta <- metadata
    rownames(meta) <- 1:nrow(meta)

    m <- as(t(pcs), "dgTMatrix")
    colnames(m) <- 1:ncol(m)

    obj <- Seurat::CreateSeuratObject(
      counts = m,
      meta.data = meta,
      assay = 'RNA',
      names.field = 1
    )

    harmony_embeddings_all <- pcs
    rownames(harmony_embeddings_all) <- 1:nrow(harmony_embeddings_all)

    obj@reductions$harmony <- Seurat::CreateDimReducObject(
      embeddings = harmony_embeddings_all,
      stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
      assay = "RNA",
      key = Seurat::Key("HARMONY", quiet = TRUE)
    )

    # Use requested number of PCs for neighbor graph if provided
    dims_use <- 1:min(ifelse(is.null(n_pcs), 20L, as.integer(n_pcs)),
                      ncol(harmony_embeddings_all))
    obj <- obj %>%
      Seurat::FindNeighbors(verbose = TRUE, reduction = 'harmony',
                            dims = dims_use, k.param = 30, nn.eps = 0)

    seurat_object <- obj
  } else if ((is.null(seurat_object) && is.null(metadata) && !is.null(pcs)) ||
             (is.null(seurat_object) && !is.null(metadata) && is.null(pcs))) {
    stop('Must provide both metadata and precomputed PCs')
  }

  ## (1) format data
  covs_keep <- test_var
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs)) covs_keep <- c(covs_keep, covs)

  if (length(names(seurat_object@graphs)) == 0) {
    stop('Must precompute graph in Seurat with FindNeighbors()')
  }
  if (is.null(graph_use)) {
    graph_use <- names(seurat_object@graphs)[[1]]
    message(glue::glue('Graph not specified. Using graph {graph_use}'))
  } else {
    if (!graph_use %in% names(seurat_object@graphs)) {
      stop(glue::glue('Graph {graph_use} not in seurat object'))
    }
  }
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(seurat_object@meta.data, dplyr::one_of(covs_keep))))
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
  data <- rcna_data
  suffix <- ''

  # batches numeric vector
  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(data.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }

  res <- list()
  covs_mat <- data.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))

  f <- as.formula(as.character(glue::glue('~0+{data$samplem_key}')))
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue::glue('^{data$samplem_key}(.*)')), '\\\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]
  s <- s[, data$samplem[[data$samplem_key]]]  # cells x samples indicator

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
      if (prevmedkurt - medkurt < 3 && i > 3) {
        message(glue::glue('stopping after {i} steps'))
        break
      }
    } else if (i == nsteps) {
      break
    }
  }

  snorm <- t(prop.table(s, 2))  # samples x cells
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]
  NAM <- snorm

  N <- nrow(NAM)
  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))
  } else {
    message('filtering based on batches kurtosis')
  }

  .batch_kurtosis <- function(NAM, batches_vec) {
    purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
      if (length(i) > 1) {
        Matrix::colMeans(NAM[i, ])
      } else if (length(i) == 1) {
        Matrix::colMeans(t(NAM[i, ]))
      }
    }) |>
      dplyr::bind_cols() |>
      apply(1, moments::kurtosis)
  }

  kurtoses <- .batch_kurtosis(NAM, batches_vec)
  threshold <- max(6, 2 * median(kurtoses))
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}'))
  keep <- which(kurtoses < threshold)

  .res_qc_nam <- list(NAM = NAM, keep = keep)

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

  .res_resid_nam <- list(
    NAM_ = scale(NAM_, center = FALSE, scale = TRUE),
    M = M,
    r = ncols_C
  )
  res[[paste0('_M', suffix)]] <- .res_resid_nam$M
  res[[paste0('_r', suffix)]] <- .res_resid_nam$r

  if (verbose) message('Decompose NAM')

  # Decide number of PCs
  n_samples <- nrow(data$samplem)
  max_rank <- min(nrow(NAM_), ncol(NAM_)) - 1L
  if (!is.null(n_pcs)) {
    npcs <- as.integer(n_pcs)
  } else {
    npcs <- max(10L, round(max_frac_pcs * n_samples))
  }
  npcs <- max(1L, min(npcs, max_rank))

  # Full vs truncated SVD
  scaled_NAM_ <- scale(NAM_, center = FALSE, scale = TRUE)
  if (npcs > 0.5 * min(dim(NAM_))) {
    svd_full <- svd(scaled_NAM_)
    svd_res <- list(
      u = svd_full$u[, seq_len(npcs), drop = FALSE],
      v = svd_full$v[, seq_len(npcs), drop = FALSE],
      d = svd_full$d[seq_len(npcs)]
    )
  } else {
    svd_res <- RSpectra::svds(scaled_NAM_, k = npcs)
  }

  # A = U D V^T
  U_df <- svd_res$u[, seq_len(npcs), drop = FALSE]
  colnames(U_df) <- paste0('PC', seq_len(npcs))
  rownames(U_df) <- rownames(NAM_)
  V_df <- svd_res$v[, seq_len(npcs), drop = FALSE]
  colnames(V_df) <- paste0('PC', seq_len(npcs))
  rownames(V_df) <- colnames(NAM_)
  .res_svd_nam <- list(U = U_df, svs = svd_res$d^2, V = V_df)

  res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
  res[[paste0('NAM_svs', suffix)]] <- .res_svd_nam$svs
  res[[paste0('NAM_varexp', suffix)]] <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
  res[[paste0('NAM_nbhdXpc', suffix)]] <- .res_svd_nam$V

  nam_res <- res

  NAMsvd <- list(
    nam_res$NAM_sampleXpc,
    nam_res$NAM_svs,
    nam_res$NAM_nbhdXpc,
    nam_res$NAM_varexp
  )
  names(NAMsvd) <- c("sampleXpc", "svs", "nbhdXpc", "varexp")

  M <- res[[paste0('_M', suffix)]]
  r <- res[[paste0('_r', suffix)]]

  yvals <- rcna_data$samplem[[test_var]]
  if (is(yvals, 'character') | is(yvals, 'factor') | is(yvals, 'integer')) {
    stop(glue::glue('test_var is of class {class(yvals)}. It must be numeric variable for association testing.'))
  }
  y <- yvals

  if (is.null(seed)) {
    set.seed(sample(1e6, 1))
  }
  if (force_permute_all) {
    batches_vec <- rep(1L, length(y))
  }

  # prep data
  U <- NAMsvd[[1]]
  sv <- NAMsvd[[2]]
  V <- NAMsvd[[3]]
  y <- scale(y)
  n <- length(y)

  if (is.null(ks)) {
    incr <- max(round(0.02 * n), 1)
    maxnpcs_default <- min(4 * incr, round(n / 5))
    ks <- seq(incr, min(maxnpcs_default, ncol(U)), incr)
    ks <- unique(sort(pmax(1L, ks)))
  }

  .reg <- function(q, k) {
    Xpc <- U[, 1:k, drop = FALSE]
    beta <- t(Xpc) %*% q
    qhat <- Xpc %*% beta
    return(list(qhat = qhat, beta = beta))
  }

  .stats <- function(yhat, ycond, k) {
    ssefull <- as.numeric(crossprod(yhat - ycond))
    ssered <- as.numeric(crossprod(ycond))
    deltasse <- ssered - ssefull
    f <- (deltasse / k) / (ssefull / n)
    p <- -pf(f, k, n - (1 + r + k), log.p = TRUE)
    r2 <- 1 - ssefull / ssered
    return(list(p = p, r2 = r2))
  }

  .minp_stats <- function(z) {
    zcond <- scale(M %*% z, center = FALSE, scale = TRUE)
    qhats <- purrr::map(ks, function(k) .reg(zcond, k)$qhat)
    .tmp <- purrr::map2(qhats, ks, function(qhat, k) .stats(qhat, zcond, k))
    ps <- purrr::map_dbl(.tmp, 'p')
    r2s <- purrr::map_dbl(.tmp, 'r2')
    k_ <- which.min(ps)
    return(list(k = ks[k_], p = ps[k_], r2 = r2s[k_]))
  }

  # get non-null f-test p-value
  .tmp <- .minp_stats(y)
  k <- .tmp$k
  p <- .tmp$p
  r2 <- .tmp$r2
  if (k == max(ks)) {
    warning(glue::glue('data supported use of {k} NAM PCs, which is the maximum considered. Consider allowing more PCs by using the "ks" argument.'))
  }

  # compute coefficients and r2 with chosen model
  ycond <- scale(M %*% y, center = FALSE, scale = TRUE)
  .tmp <- .reg(ycond, k)
  yhat <- .tmp$qhat
  beta <- .tmp$beta
  r2_perpc <- (beta / as.numeric(sqrt(crossprod(ycond))))**2

  ncorrs <- V[, 1:k, drop = FALSE] %*% (sqrt(sv[1:k]) * beta / n)
  rownames(ncorrs) <- rownames(V)

  set.seed(seed)
  y_ <- conditional_permutation(batches_vec, y, Nnull)
  .tmp <- apply(y_, 2, .minp_stats)
  nullminps <- purrr::map_dbl(.tmp, 'p')
  nullr2s <- purrr::map_dbl(.tmp, 'r2')

  pfinal <- (sum(nullminps <= p + 1e-8) + 1) / (Nnull + 1)
  if (sum(nullminps <= p + 1e-8) == 0) {
    warning('global association p-value attained minimal possible value. Consider increasing Nnull')
  }

  # get neighborhood fdrs if requested
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
    y_ <- y_[, 1:Nnull, drop = FALSE]
    ycond_ <- scale(M %*% y_, center = FALSE, scale = TRUE)
    gamma_ <- crossprod(U[, 1:k, drop = FALSE], ycond_)
    nullncorrs <- abs(V[, 1:k, drop = FALSE] %*% (sqrt(sv[1:k]) * (gamma_ / n)))

    maxcorr <- max(abs(ncorrs))
    fdr_thresholds <- seq(maxcorr / 4, maxcorr, maxcorr / 400)
    fdr_vals <- empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
    fdrs <- data.frame(
      threshold = head(fdr_thresholds, -1),
      fdr = fdr_vals,
      num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t))
    )
    # find minimal thresholds giving desired FDRs
    if (min(fdrs$fdr) <= 0.05) fdr_5p_t  <- min(subset(fdrs, fdr < 0.05)$threshold)
    if (min(fdrs$fdr) <= 0.10) fdr_10p_t <- min(subset(fdrs, fdr < 0.10)$threshold)
    if (min(fdrs$fdr) <= 0.20) fdr_20p_t <- min(subset(fdrs, fdr < 0.20)$threshold)
    if (min(fdrs$fdr) <= 0.30) fdr_30p_t <- min(subset(fdrs, fdr < 0.30)$threshold)
    if (min(fdrs$fdr) <= 0.40) fdr_40p_t <- min(subset(fdrs, fdr < 0.40)$threshold)
    if (min(fdrs$fdr) <= 0.50) fdr_50p_t <- min(subset(fdrs, fdr < 0.50)$threshold)
  }

  res <- list(
    p = pfinal,
    nullminps = nullminps,
    k = k,
    ncorrs = ncorrs,
    fdrs = fdrs,
    fdr_5p_t = fdr_5p_t,
    fdr_10p_t = fdr_10p_t,
    yhat = yhat,
    ycond = ycond,
    ks = ks,
    beta = beta,
    r2 = r2,
    r2_perpc = r2_perpc,
    nullr2_mean = mean(nullr2s),
    nullr2_std = sd(nullr2s)
  )

  if (return_nam) {
    res[['NAM_embeddings']] <- nam_res$NAM_nbhdXpc
    res[['NAM_loadings']] <- nam_res$NAM_sampleXpc
    res[['NAM_svs']] <- nam_res$NAM_svs
    res[['NAM_varexp']] <- nam_res$NAM_varexp
    res[['NAM']] <- t(NAM)
  }

  seurat_object[['cna']] <- Seurat::CreateDimReducObject(
    embeddings = res$NAM_embeddings,
    loadings = res$NAM_loadings,
    stdev = res$NAM_svs,
    assay = assay,
    key = key,
    misc = res
  )

  seurat_object@meta.data$cna_ncorrs <- ncorrs[colnames(seurat_object), , drop = TRUE]

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

  ## ---- NEW: if scLASER object passed, save and return it instead ----
  if (!is.null(object) && return_nam) {
    object@NAM_matrix <- res[['NAM']]
    object@nam_pcs    <- res[['NAM_sampleXpc']]
    return(object)
  }

  return(seurat_object)
}
