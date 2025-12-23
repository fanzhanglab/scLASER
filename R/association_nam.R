#' Association testing between NAM embeddings and a sample-level variable
#'
#' Computes Neighborhood Aggregation Matrix (NAM) embeddings from a Seurat
#' nearest-neighbor graph, residualizes NAM with respect to optional covariates,
#' performs SVD to obtain NAM PCs, and tests association between the NAM subspace
#' and a numeric sample-level variable using a min-p scan over candidate numbers
#' of PCs with conditional permutations.
#'
#' This function can either (A) use a provided \pkg{Seurat} object that already
#' contains a neighbor graph, or (B) construct a minimal Seurat object internally
#' from \code{obj@metadata} and \code{obj@harmony} and then build the graph.
#'
#' @param obj A \code{\linkS4class{scLASER}} object. If \code{seurat_object} is
#'   \code{NULL}, \code{obj} must have both \code{@metadata} and \code{@harmony}
#'   populated so a minimal Seurat object can be constructed.
#'
#' @param seurat_object Optional. A \code{Seurat} object with a precomputed
#'   cell-cell neighbor graph (e.g., after \code{Seurat::FindNeighbors()}).
#'   If \code{NULL}, a minimal Seurat object is built from \code{obj}.
#'
#' @param test_var Character scalar. Column name in sample-level metadata
#'   (derived from \code{seurat_object@meta.data}) giving a numeric variable
#'   to test for association with NAM embeddings.
#'
#' @param samplem_key Character scalar. Column name identifying samples in the
#'   cell-level metadata (and used to derive sample-level metadata).
#'
#' @param graph_use Character scalar. Name of the graph in
#'   \code{seurat_object@graphs} to use. Default is \code{"RNA_snn"}.
#'
#' @param batches Optional character vector of column names (sample-level) used
#'   to define permutation blocks for conditional permutations. If \code{NULL},
#'   all samples are permuted together.
#'
#' @param covs Optional character vector of sample-level covariate column names
#'   to regress out (residualize) before association testing.
#'
#' @param nsteps Integer or \code{NULL}. Number of diffusion steps used to build
#'   NAM. If \code{NULL}, an adaptive stopping rule is used.
#'
#' @param verbose Logical. Print progress messages.
#'
#' @param assay Character or \code{NULL}. Assay name used when creating the
#'   NAM reduction object internally.
#'
#' @param key Character. Key prefix for NAM PC names in the created reduction.
#'
#' @param maxnsteps Integer. Maximum diffusion steps when using adaptive
#'   stopping.
#'
#' @param max_frac_pcs Numeric in (0, 1]. If \code{n_pcs} is \code{NULL}, the
#'   default number of PCs is computed as \code{max(10, round(max_frac_pcs * N))},
#'   where \code{N} is the number of samples.
#'
#' @param n_pcs Integer or \code{NULL}. Explicit number of NAM PCs to compute
#'   (and also used to build neighbors in the internally constructed Seurat
#'   object). If \code{NULL}, uses \code{max_frac_pcs} heuristic.
#'
#' @param ks Optional integer vector. Candidate numbers of NAM PCs to scan in the
#'   min-p procedure. If \code{NULL}, a default grid is constructed from sample
#'   size and available PCs.
#'
#' @param Nnull Integer. Number of conditional permutation replicates.
#'
#' @param force_permute_all Logical. If \code{TRUE}, ignores \code{batches} and
#'   permutes all samples together.
#'
#' @param seed Integer or \code{NULL}. Random seed used for permutations. If
#'   \code{NULL}, a random seed is chosen.
#'
#' @param return_nam Logical. If \code{TRUE}, stores NAM-related matrices
#'   (embeddings/loadings/variance explained and NAM itself) into the returned
#'   object.
#'
#' @return A \code{\linkS4class{scLASER}} object with updated slots:
#' \itemize{
#'   \item \code{@nam_pcs}: NAM neighborhood-by-PC matrix.
#'   \item \code{@NAM_matrix}: NAM sample-by-neighborhood matrix (transposed as stored).
#' }
#' (Additional association results may be stored in the Seurat reduction
#' internally, depending on implementation.)
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}},
#'   \code{\link[Seurat]{CreateSeuratObject}},
#'   \code{\link[Seurat]{CreateDimReducObject}}
#'
#' @examples
#' \donttest{
#' ## Example workflow (requires PCs + Harmony already computed in the scLASER object)
#' ## 1) Load a prepared scLASER object that already contains metadata and PCs.
#' ##    (You typically create this in your own pipeline and save it as an .rds.)
#' # obj <- readRDS("path/to/precomputed_scLASER_object.rds")
#'
#' ## 2) (Optional) If Harmony embeddings are not present, compute them.
#' ##    This assumes obj@pcs is populated and obj@metadata contains `sample_id`
#' ##    (or another batch column).
#' # obj <- run_harmony(obj, batch_col = "sample_id")
#'
#' ## 3) Run association testing.
#' ##    If you do not pass a Seurat object, this function will internally
#' ##    construct a minimal Seurat object from obj@metadata + obj@harmony.
#' # obj <- association_nam_scLASER(
#' #   obj          = obj,
#' #   seurat_object = NULL,
#' #   test_var     = "age",
#' #   samplem_key  = "sample_id",
#' #   batches      = c("sample_id"),
#' #   covs         = c("sex"),
#' #   Nnull        = 200,
#' #   seed         = 1
#' # )
#'
#' ## 4) Inspect outputs
#' # dim(obj@nam_pcs)
#' # dim(obj@NAM_matrix)
#' }

#' @export


association_nam_scLASER <- function(obj,seurat_object = NULL,
                                    #  metadata = NULL,
                                    # pcs = NULL,
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
                                    seed = 1234,
                                    return_nam = TRUE) {

  if (!is.null(seurat_object)) {
    message("Using provided Seurat object for analysis...")

  } else if (!is.null(obj@metadata) && !is.null(obj@harmony)) {
    # Build a minimal Seurat object from metadata + precomputed Harmony embeddings
    meta <- obj@metadata
    rownames(meta) <- seq_len(nrow(meta))

    m <- as(t(obj@harmony), "dgTMatrix")
    colnames(m) <- seq_len(ncol(m))

    obj1 <- Seurat::CreateSeuratObject(
      counts = m,
      meta.data = meta,
      assay = "RNA",
      names.field = 1
    )

    harmony_embeddings_all <- obj@harmony
    rownames(harmony_embeddings_all) <- seq_len(nrow(harmony_embeddings_all))

    obj1@reductions$harmony <- Seurat::CreateDimReducObject(
      embeddings = harmony_embeddings_all,
      stdev = as.numeric(apply(harmony_embeddings_all, 2, stats::sd)),
      assay = "RNA",
      key = Seurat::Key("HARMONY", quiet = TRUE)
    )

    # Use requested number of PCs for neighbor graph if provided
    dims_use <- 1:min(
      ifelse(is.null(n_pcs), 20L, as.integer(n_pcs)),
      ncol(harmony_embeddings_all)
    )

    obj1 <- Seurat::FindNeighbors(
      object = obj1,
      verbose = TRUE,
      reduction = "harmony",
      dims = dims_use,
      k.param = 30,
      nn.eps = 0
    )

    seurat_object <- obj1

  } else {
    stop(
      "Must provide either:\n",
      "  (1) a non-NULL `seurat_object`, or\n",
      "  (2) an scLASER object with both `@metadata` and `@harmony` populated."
    )
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
  colnames(s) <- gsub(as.character(glue::glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
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
  y_ <- .conditional_permutation(batches_vec, y, Nnull)
  .tmp <- apply(y_, 2, .minp_stats)
  nullminps <- purrr::map_dbl(.tmp, 'p')
  nullr2s <- purrr::map_dbl(.tmp, 'r2')

  pfinal <- (sum(nullminps <= p + 1e-8) + 1) / (Nnull + 1)
  if (sum(nullminps <= p + 1e-8) == 0) {
    warning('global association p-value attained minimal possible value. Consider increasing Nnull')
  }


  res <- list(
    p = pfinal,
    nullminps = nullminps,
    k = k,
    ncorrs = ncorrs,
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



  obj@nam_pcs =nam_res$NAM_nbhdXpc
  obj@NAM_matrix = t(NAM)


  return(obj)
}
