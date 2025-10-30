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
generate_dummy_data <- function(
  n_cells = 3000,                  # baseline cells per major cell type per sample
  sd_celltypes = 0.10,             # relative sd for counts
  n_major_cell_types = 7,
  n_minor_cell_types = 3,
  relative_abundance = 0.10,       # minor vs major baseline ratio
  n_major_interact_celltypes = 1,  # how many majors are "interacting"
  n_minor_interact_celltypes = 1,  # how many minors are "interacting"
  n_individuals = 30,
  n_batchs = 4,

  interaction_feature = "visit",   # kept for labeling
  time_points = 4,                 # >= 2; works for 3+
  test_var = "disease",
  prop_disease = 0.50,

  fc_interact = 0.10,              # effect magnitude used by defaults below
  interaction_type = c("specific","differential","opposite"),
  seed = 1234,

  visit_effects_progressor = NULL, # multiplicative change: +0.1 means +10% vs baseline
  visit_effects_control    = NULL, # default 0s
  direction_by_cluster     = NULL  # +1/-1 for interacting clusters, recycled as needed
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
    disease = disease_vec,
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
  meta$interaction <- paste0(interaction_feature, ":", test_var)

  # --- Default visit effects (length = time_points) ----
  if (is.null(visit_effects_progressor)) {
    if (interaction_type == "specific") {
      ve <- rep(0, time_points)
      if (time_points >= 2) ve[2] <- fc_interact  # bump V1 only
      visit_effects_progressor <- ve
    } else if (interaction_type == "differential") {
      visit_effects_progressor <- rep(fc_interact, time_points)
    } else { # "opposite": alternate +/-
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

  # --- Baseline counts per (sample, cell_type) ----
  # Major types centered around n_cells; minor types around n_cells * relative_abundance
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

  # --- Apply effects only to interacting cell types ----
  # effect = direction * visit_effect (progressor vs control)
  is_interacting <- counts_long$cell_type %in% interact_cell_types
  eff_vec <- numeric(nrow(counts_long))
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

  # adjusted counts (non-negative)
  counts_long$adj_count <- pmax(0L, round(counts_long$count * (1 + eff_vec)))

  # --- Expand to per-cell rows and attach metadata ----
  # Note: Simple & clear; for large sizes, this is big by design (single-cell rows)
  rep_each <- function(x, times) if (length(x) == 0) x else rep(x, times = times)
  expanded <- data.frame(
    sample_id = rep_each(counts_long$sample_id, counts_long$adj_count),
    cell_type = rep_each(counts_long$cell_type, counts_long$adj_count),
    stringsAsFactors = FALSE
  )

  # Merge metadata for each cell
  dummy_data <- merge(expanded,
                      meta[, c("sample_id","subject_id","visit","sex","disease","age","bmi","batch","interaction")],
                      by = "sample_id", sort = FALSE)

  # Shuffle rows for realism
  if (nrow(dummy_data) > 1) {
    dummy_data <- dummy_data[sample.int(nrow(dummy_data)), , drop = FALSE]
    rownames(dummy_data) <- NULL
  }

  dummy_data
}

