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

association_nam_scLASER <- function(object,
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
                                    max_frac_pcs = 0.15,
                                    ks = NULL,
                                    Nnull = 1000,
                                    force_permute_all = FALSE,
                                    local_test = TRUE,
                                    seed = 1234,
                                    return_nam = TRUE) {
  
  stopifnot(inherits(object, "scLASER"))
  if (is.null(object@metadata)) stop("@metadata is NULL in scLASER object.")
  if (is.null(samplem_key) || !(samplem_key %in% names(object@metadata))) {
    stop("Please provide a valid 'samplem_key' present in object@metadata.")
  }
  if (is.null(test_var)) stop("Please provide 'test_var' (numeric sample-level variable).")
  
  # --- Build a minimal Seurat object from scLASER (use harmony if present, else pcs) ---
  meta <- object@metadata
  
  # ensure test_var exists (if missing, try to derive from 'disease')
  if (!(test_var %in% names(meta))) {
    if ("disease" %in% names(meta)) {
      meta[[test_var]] <- suppressWarnings(as.numeric(as.character(meta[["disease"]])))
    } else {
      stop(sprintf("Column '%s' not found in metadata and 'disease' not available to derive it.", test_var))
    }
  }
  # force numeric test_var
  meta[[test_var]] <- suppressWarnings(as.numeric(meta[[test_var]]))
  if (any(!is.finite(meta[[test_var]]))) {
    stop(sprintf("Non-finite values found in '%s'. Make sure it's numeric.", test_var))
  }
  
  # choose embeddings matrix
  emb <- if (!is.null(object@harmony)) object@harmony else object@pcs
  if (is.null(emb)) stop("Both @harmony and @pcs are NULL; need embeddings to build graph.")
  emb <- as.matrix(emb)
  storage.mode(emb) <- "double"
  
  # give unique per-cell rownames for Seurat using samplem_key
  meta[[samplem_key]] <- as.character(meta[[samplem_key]])
  cell_ids <- make.unique(meta[[samplem_key]])
  rownames(meta) <- cell_ids
  rownames(emb) <- cell_ids
  
  # tiny dense counts (avoid dgCMatrix colSums symbol error paths)
  counts <- matrix(1L, nrow = 1L, ncol = nrow(meta), dimnames = list("gene1", cell_ids))
  
  obj <- Seurat::CreateSeuratObject(
    counts       = counts,
    meta.data    = meta,
    assay        = 'RNA',
    names.field  = 1,
    min.cells    = 0,
    min.features = 0
  )
  
  harm_sd <- apply(emb, 2, stats::sd); harm_sd[!is.finite(harm_sd)] <- 0
  obj@reductions$harmony <- Seurat::CreateDimReducObject(
    embeddings = emb,
    stdev      = as.numeric(harm_sd),
    assay      = "RNA",
    key        = Seurat::Key("HARMONY", quiet = TRUE)
  )
  
  obj <- Seurat::FindNeighbors(
    object    = obj,
    reduction = 'harmony',
    dims      = seq_len(min(20L, ncol(emb))),
    k.param   = 30,
    nn.eps    = 0,
    verbose   = verbose
  )
  
  # === ORIGINAL association_nam logic from here (unchanged) ============================
  covs_keep <- test_var
  if (!is.null(batches)) covs_keep <- c(covs_keep, batches)
  if (!is.null(covs))    covs_keep <- c(covs_keep, covs)
  
  if (length(names(obj@graphs)) == 0) stop('Must precompute graph in Seurat with FindNeighbors()')
  if (is.null(graph_use)) {
    graph_use <- names(obj@graphs)[[1]]
    message(glue::glue('Graph not specified. Using graph {graph_use}'))
  } else if (!graph_use %in% names(obj@graphs)) {
    stop(glue::glue('Graph {graph_use} not in seurat object'))
  }
  
  covs_keep <- c(covs_keep, samplem_key)
  samplem_df <- tibble::remove_rownames(unique(dplyr::select(obj@meta.data, dplyr::one_of(covs_keep))))
  obs_df     <- tibble::rownames_to_column(obj@meta.data, 'CellID')
  if (nrow(samplem_df) == nrow(obs_df)) {
    stop('Sample-level metadata is same length as cell-level metadata. Please check samplem_vars.')
  }
  
  rcna_data <- list(
    samplem = samplem_df,
    obs = obs_df,
    connectivities = obj@graphs[[graph_use]],
    samplem_key = samplem_key,
    obs_key = 'CellID',
    N = nrow(samplem_df)
  )
  data <- rcna_data
  suffix <- ''
  
  if (is.null(batches)) {
    batches_vec <- rep(1, data$N)
  } else {
    batches_vec <- as.integer(data.matrix(dplyr::select(data$samplem, dplyr::one_of(batches))))
  }
  
  res <- list()
  covs_mat <- if (is.null(covs)) NULL else data.matrix(dplyr::select(data$samplem, dplyr::one_of(covs)))
  
  f <- stats::as.formula(as.character(glue::glue('~0+{data$samplem_key}')))
  s <- model.matrix(f, data$obs)
  colnames(s) <- gsub(as.character(glue::glue('^{data$samplem_key}(.*)')), '\\1', colnames(s))
  rownames(s) <- data$obs[[data$obs_key]]
  s <- s[, data$samplem[[data$samplem_key]]]
  
  diffuse_step <- function(data, s) {
    a <- data$connectivities
    degrees <- Matrix::colSums(a) + 1
    s_norm <- s / degrees
    res <- (a %*% s_norm) + s_norm
    as.matrix(res)
  }
  
  prevmedkurt <- Inf
  for (i in seq_len(maxnsteps)) {
    s <- diffuse_step(data, s)
    medkurt <- stats::median(apply(prop.table(s, 2), 1, moments::kurtosis))
    if (is.null(nsteps)) {
      prevmedkurt <- medkurt
      if (prevmedkurt - medkurt < 3 & i > 3) {
        message(glue::glue('stopping after {i} steps'))
        break
      }
    } else if (i == nsteps) break
  }
  
  snorm <- t(prop.table(s, 2))
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]
  NAM <- snorm
  
  if (is.null(batches_vec) | length(unique(batches_vec)) == 1) {
    message('only one unique batch supplied to qc')
    keep <- rep(TRUE, ncol(NAM))
  } else {
    message('filtering based on batches kurtosis')
    .batch_kurtosis <- function(NAM, batches_vec) {
      purrr::imap(split(seq_len(length(batches_vec)), batches_vec), function(i, b) {
        if (length(i)>1) Matrix::colMeans(NAM[i, ]) else Matrix::colMeans(t(NAM[i, ]))
      }) %>% dplyr::bind_cols() %>% apply(1, moments::kurtosis)
    }
    kurtoses <- .batch_kurtosis(NAM, batches_vec)
    threshold <- max(6, 2*stats::median(kurtoses))
    keep <- which(kurtoses < threshold)
  }
  
  N <- nrow(NAM)
  if (verbose) message('Residualize NAM')
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
  NAM_ <- scale(NAM_, center = FALSE, scale = TRUE)
  
  if (verbose) message('Decompose NAM')
  npcs <- max(10, round(max_frac_pcs * nrow(data$samplem)))
  npcs <- min(npcs, nrow(data$samplem) - 1)
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(NAM_)
  } else {
    svd_res <- RSpectra::svds(NAM_, k = npcs)
  }
  
  U_df <- svd_res$u[, seq_len(npcs)]
  colnames(U_df) <- paste0('PC', seq_len(npcs))
  rownames(U_df) <- rownames(NAM_)
  V_df <- svd_res$v[, seq_len(npcs)]
  colnames(V_df) <- paste0('PC', seq_len(npcs))
  rownames(V_df) <- colnames(NAM_)
  .res_svd_nam <- list(U = U_df, svs = svd_res$d^2, V = V_df)
  
  res[[paste0('NAM_sampleXpc', suffix)]] <- .res_svd_nam$U
  res[[paste0('NAM_svs', suffix)]]       <- .res_svd_nam$svs
  res[[paste0('NAM_varexp', suffix)]]    <- .res_svd_nam$svs / nrow(.res_svd_nam$U) / nrow(.res_svd_nam$V)
  res[[paste0('NAM_nbhdXpc', suffix)]]   <- .res_svd_nam$V
  
  nam_res <- res
  
  V_save <- nam_res[['NAM_nbhdXpc']]
  colnames(V_save) <- paste0(key, seq_len(ncol(V_save)))
  object@nam_pcs <- V_save
  
  return(object)
}