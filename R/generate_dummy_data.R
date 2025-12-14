#' Title
#'
#' @param n_cells
#' @param sd_celltypes
#' @param n_major_cell_types
#' @param n_minor_cell_types
#' @param relative_abundance
#' @param n_major_interact_celltypes
#' @param n_minor_interact_celltypes
#' @param n_individuals
#' @param n_batchs
#' @param interaction_feature
#' @param time_points
#' @param test_var
#' @param prop_disease
#' @param fc_interact
#' @param interaction_type
#' @param seed
#' @param visit_effects_progressor
#' @param visit_effects_control
#' @param direction_by_cluster
#'
#' @return
#' @export
#'
#' @examples
#' Simulate longitudinal single-cell data with interacting cell types
#'
#' @param n_cells Baseline cells per major cell type per sample.
#' @param sd_celltypes Relative SD for counts.
#' @param n_major_cell_types Number of major cell types.
#' @param n_minor_cell_types Number of minor cell types.
#' @param relative_abundance Baseline minor:major ratio.
#' @param n_major_interact_celltypes How many major types are "interacting".
#' @param n_minor_interact_celltypes How many minor types are "interacting".
#' @param n_individuals Number of subjects.
#' @param n_batchs Number of batches.
#' @param interaction_feature Name of the time/visit variable (e.g. "visit").
#' @param time_points Number of longitudinal time points (>= 2).
#' @param test_var Name of the disease/exposure variable (e.g. "disease").
#' @param prop_disease Proportion of diseased subjects.
#' @param fc_interact Effect magnitude used by the default effect patterns.
#' @param interaction_type One of "specific", "differential", or "opposite"
#'   for the default longitudinal pattern (when design matrices are not given).
#' @param seed Random seed.
#' @param visit_effects_progressor Optional vector of length `time_points`
#'   for default progressor visit effects (ignored if design matrices given).
#' @param visit_effects_control Optional vector of length `time_points`
#'   for default control visit effects (ignored if design matrices given).
#' @param direction_by_cluster Optional vector of +1 / -1 per interacting
#'   cell type (ignored if design matrices given).
#' @param effect_mat_progressor Optional numeric matrix of dimension
#'   `n_cell_types x time_points` giving effect sizes for progressors.
#'   Rows correspond to cell types (by name or order), columns to visits
#'   (0..time_points-1). If provided (along with `effect_mat_control`),
#'   overrides `visit_effects_*` and `direction_by_cluster`.
#' @param effect_mat_control Optional numeric matrix of same dimension as
#'   `effect_mat_progressor` giving effects for controls.
#'
#' @return A data.frame of simulated single-cell metadata.
#' @export
generate_dummy_data <- function(
    n_cells = 3000,                  # baseline cells per major cell type per sample
    sd_celltypes = 0.10,             # relative sd for counts
    relative_abundance = 0.10,       # minor vs major baseline ratio
    n_individuals = 30,
    n_batchs = 4,

    interaction_feature = "visit",   # kept for labeling
    time_points = 4,                 # >= 2; works for 3+
    test_var = "disease",
    prop_disease = 0.50,

    fc_interact = 0.10,              # effect magnitude used by defaults below
    interaction_type = c("specific","differential","opposite"),
    seed = 1234,


    effect_mat_progressor    = NULL, # NEW: cell_type x time_points matrix for progressors
    effect_mat_control       = NULL  # NEW: same for controls
) {
  set.seed(seed)
  interaction_type <- match.arg(interaction_type)

  # --- Basic setup ----
  n_cell_types <- n_major_cell_types + n_minor_cell_types
  stopifnot(time_points >= 2, n_cell_types >= 1)

  cell_types <- LETTERS[seq_len(n_cell_types)]
  major_idx  <- seq_len(n_major_cell_types)
  minor_idx  <- if (n_minor_cell_types > 0) (n_major_cell_types + seq_len(n_minor_cell_types)) else integer(0)

  # which cell types are "interacting" (first some majors, last some minors)
  interact_idx <- c(
    head(major_idx, n_major_interact_celltypes),
    tail(seq_len(n_cell_types), n_minor_interact_celltypes)
  )
  interact_idx <- intersect(interact_idx, seq_len(n_cell_types))  # guard rails
  interact_cell_types <- cell_types[interact_idx]

  # --- Subjects and visits ----
  subject_id <- paste0("SUB_", seq_len(n_individuals))
  disease_vec <- c(rep(1L, round(n_individuals * prop_disease)),
                   rep(0L, n_individuals - round(n_individuals * prop_disease)))
  disease_vec <- sample(disease_vec, n_individuals)

  sex_vec <- sample(c(0L, 1L), n_individuals, replace = TRUE)          # 0/1
  age_vec <- sample(18:60, n_individuals, replace = TRUE)
  bmi_vec <- sample(15:35, n_individuals, replace = TRUE)
  batch_vec <- rep(seq_len(n_batchs), length.out = n_individuals)

  subjects <- data.frame(
    subject_id = subject_id,
    sex   = sex_vec,
    disease = disease_vec,           # canonical exposure storage
    age   = age_vec,
    bmi   = bmi_vec,
    batch = factor(batch_vec),
    stringsAsFactors = FALSE
  )

  visits <- data.frame(
    subject_id = rep(subject_id, each = time_points),
    visit      = rep(0:(time_points - 1), times = n_individuals),
    stringsAsFactors = FALSE
  )
  visits$sample_id <- paste0(visits$subject_id, "_V", visits$visit)

  meta <- merge(visits, subjects, by = "subject_id", sort = FALSE)

  # Make sure a column named `test_var` exists (even if test_var != "disease")
  if (!identical(test_var, "disease")) {
    meta[[test_var]] <- meta[["disease"]]
  }

  # Label and INTERACTION TERM (persist to output)
  meta$interaction    <- paste0(interaction_feature, ":", test_var)
  meta$interact_term  <- as.integer(meta[[interaction_feature]]) * as.integer(meta[[test_var]])

  # --- Baseline counts per (sample, cell_type) ----
  one_sample_counts <- function() {
    major_counts <- round(runif(n_major_cell_types,
                                min = n_cells * (1 - sd_celltypes),
                                max = n_cells * (1 + sd_celltypes)))
    minor_counts <- if (n_minor_cell_types > 0) {
      round(runif(n_minor_cell_types,
                  min = n_cells * relative_abundance * (1 - sd_celltypes),
                  max = n_cells * relative_abundance * (1 + sd_celltypes)))
    } else integer(0)
    c(major_counts, minor_counts)
  }

  counts_list <- replicate(nrow(meta), one_sample_counts(), simplify = FALSE)
  counts_df <- do.call(rbind, counts_list)
  colnames(counts_df) <- cell_types

  counts_long <- reshape(
    data.frame(sample_id = meta$sample_id, counts_df, check.names = FALSE),
    varying = cell_types, v.names = "count", timevar = "cell_type",
    times = cell_types, direction = "long"
  )
  rownames(counts_long) <- NULL

  # Merge disease/visit so we can apply effects
  counts_long <- merge(counts_long,
                       meta[, c("sample_id", "visit", "disease")],
                       by = "sample_id", sort = FALSE)

  # --- Build effect structure ----
  eff_vec <- numeric(nrow(counts_long))

  if (!is.null(effect_mat_progressor) || !is.null(effect_mat_control)) {
    ## ---- NEW: use design matrices if provided ----
    if (is.null(effect_mat_progressor) || is.null(effect_mat_control)) {
      stop("Both effect_mat_progressor and effect_mat_control must be provided if using design matrices.")
    }

    effect_mat_progressor <- as.matrix(effect_mat_progressor)
    effect_mat_control    <- as.matrix(effect_mat_control)

    if (ncol(effect_mat_progressor) != time_points ||
        ncol(effect_mat_control)    != time_points) {
      stop("effect_mat_* must have ncol = time_points.")
    }
    if (nrow(effect_mat_progressor) != n_cell_types ||
        nrow(effect_mat_control)    != n_cell_types) {
      stop("effect_mat_* must have nrow = n_major_cell_types + n_minor_cell_types.")
    }

    # Align rows to cell_types via rownames if present; otherwise assume in order
    if (!is.null(rownames(effect_mat_progressor))) {
      if (!all(sort(rownames(effect_mat_progressor)) == sort(cell_types))) {
        stop("Row names of effect_mat_progressor must match cell types: ",
             paste(cell_types, collapse = ", "))
      }
      effect_mat_progressor <- effect_mat_progressor[cell_types, , drop = FALSE]
    } else {
      rownames(effect_mat_progressor) <- cell_types
    }

    if (!is.null(rownames(effect_mat_control))) {
      if (!all(sort(rownames(effect_mat_control)) == sort(cell_types))) {
        stop("Row names of effect_mat_control must match cell types: ",
             paste(cell_types, collapse = ", "))
      }
      effect_mat_control <- effect_mat_control[cell_types, , drop = FALSE]
    } else {
      rownames(effect_mat_control) <- cell_types
    }

    # Map each row in counts_long to (cell_type, visit) -> effect
    row_idx <- match(counts_long$cell_type, cell_types)
    col_idx <- counts_long$visit + 1L  # visits are 0..time_points-1

    is_prog <- counts_long$disease == 1L

    eff_vec[is_prog]  <- effect_mat_progressor[cbind(row_idx[is_prog],  col_idx[is_prog])]
    eff_vec[!is_prog] <- effect_mat_control[   cbind(row_idx[!is_prog], col_idx[!is_prog])]

  } else {
    ## ---- OLD BEHAVIOR: vector visit effects + direction_by_cluster ----

    # --- Default visit effects (length = time_points) ----
    if (is.null(visit_effects_progressor)) {
      if (interaction_type == "specific") {
        ve <- rep(0, time_points)
        if (time_points >= 2) ve[2] <- fc_interact  # bump V1 only
        visit_effects_progressor <- ve
      } else if (interaction_type == "differential") {
        visit_effects_progressor <- rep(fc_interact, time_points)
      } else { # "opposite": alternate +/- starting at V0
        ve <- rep(0, time_points)
        ve[seq(1, time_points, by = 2)] <- +fc_interact  # V0, V2, ...
        if (time_points >= 2) ve[seq(2, time_points, by = 2)] <- -fc_interact # V1, V3, ...
        visit_effects_progressor <- ve
      }
    }
    if (is.null(visit_effects_control)) {
      visit_effects_control <- rep(0, time_points)
    }
    stopifnot(length(visit_effects_progressor) == time_points,
              length(visit_effects_control)    == time_points)

    # Direction per interacting cluster
    if (is.null(direction_by_cluster)) {
      direction_by_cluster <- rep(1L, max(1, length(interact_cell_types)))
    }
    # recycle to interacting set length
    direction_by_cluster <- rep(direction_by_cluster, length.out = length(interact_cell_types))
    names(direction_by_cluster) <- interact_cell_types

    # Apply effects only to interacting cell types
    is_interacting <- counts_long$cell_type %in% interact_cell_types
    if (any(is_interacting)) {
      # map direction per cell type
      dir_map <- setNames(direction_by_cluster, nm = names(direction_by_cluster))
      dir_ct  <- unname(dir_map[counts_long$cell_type[is_interacting]])
      # visit effects by group
      v_idx   <- counts_long$visit[is_interacting] + 1L
      is_prog <- counts_long$disease[is_interacting] == 1L
      ve      <- ifelse(is_prog, visit_effects_progressor[v_idx], visit_effects_control[v_idx])
      eff_vec[is_interacting] <- dir_ct * ve
    }
  }

  # --- Apply effects to counts ----
  counts_long$adj_count <- pmax(0L, round(counts_long$count * (1 + eff_vec)))

  # --- Expand to per-cell rows and attach metadata (INCLUDING interact_term) ----
  rep_each <- function(x, times) if (length(x) == 0) x else rep(x, times = times)
  expanded <- data.frame(
    sample_id = rep_each(counts_long$sample_id, counts_long$adj_count),
    cell_type = rep_each(counts_long$cell_type, counts_long$adj_count),
    stringsAsFactors = FALSE
  )

  # Merge metadata; keep interaction + interact_term
  keep_cols <- c("sample_id","subject_id","visit","sex","disease","age","bmi","batch",
                 "interaction","interact_term")
  dummy_data <- merge(expanded, meta[, keep_cols], by = "sample_id", sort = FALSE)

  # Shuffle rows for realism
  if (nrow(dummy_data) > 1) {
    dummy_data <- dummy_data[sample.int(nrow(dummy_data)), , drop = FALSE]
    rownames(dummy_data) <- NULL
  }

  dummy_data
}
