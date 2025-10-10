generate_dummy_data_time_3plus <- function(n_cells = 3000, # cells of major cell types per individual
         sd_celltypes = 0.1,  # standard deviation of number of cells
         n_major_cell_types = 7,
         n_minor_cell_types = 3,
         relative_abundance = 0.1, # ratio between major and rare
         n_major_interact_celltypes = 1,
         n_minor_interact_celltypes = 1,
         n_individuals = 30, # total individuals
         n_batchs = 4,

         interaction_feature = c("visit"),
         time_points = 2,

         test_var = c("disease"),
         prop_disease = 0.5,

         fc_interact = 0.1, # additional proportion of interacted cells which are from people with interacted group compared to non-interacted cell types
         interaction_type = c("specific","differential","opposite"),
         seed = 1234,

         # NEW arguments for 3+ time points
         visit_effects_progressor = NULL,  # numeric vector length = time_points (multiplicative changes vs baseline)
         visit_effects_control    = NULL,  # optional vector for controls; if NULL, 0s are used
         direction_by_cluster     = NULL   # optional vector of +1/-1 for interacting clusters; if NULL, all +1

) {
  n_cell_types = n_major_cell_types + n_minor_cell_types

  # Generate unique subject IDs
  subject_id <- paste0("SUB_", 1:n_individuals)
  sample_id <- paste0(subject_id,rep(paste0("_V",0:(time_points-1)),each=n_individuals))
  # set.seed(seed)
  sex = sample(
    c(rep(1, round(length(unique(subject_id)) * 0.5)),
      rep(0, n_individuals - round(length(unique(subject_id)) * 0.5)))
  )
  # set.seed(seed*2)
  disease = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_disease)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_disease))
    )
  )
  if (length(disease) != length(subject_id)) {
    disease = c(disease, rep(1,length(subject_id)-length(disease)))
  }
  # set.seed(seed*3)
  age <- sample(c("male","female"), n_individuals, replace = TRUE)
  # set.seed(seed*3)
  age <- sample(18:60, n_individuals, replace = TRUE)
  # set.seed(seed*4)
  bmi <- sample(15:35, n_individuals, replace = TRUE)
  batch <- rep(rep(1:n_batchs), length(subject_id))[1:length(subject_id)]
  dummy_data = data.frame(subject_id = subject_id,
                          sample_id = sample_id,
                          visit = rep(0:(time_points-1),each=n_individuals),
                          sex = rep(as.integer(sex, levels = c(0, 1)),time_points),
                          disease = rep(factor(disease, levels = c(0, 1)),time_points),
                          age = rep(age,time_points),
                          batch = rep(factor(batch, levels = c(0:max(batch))),time_points),
                          bmi = rep(bmi,time_points),
                          interaction = paste0(interaction_feature,":",test_var)
  )
  print(head(dummy_data))
  dummy_data$interact_term = dummy_data[,interaction_feature] * (as.integer(dummy_data[,test_var]))
  print(head(dummy_data))
  # Major and rare cell type counts
  ## major_cell_types <- ceiling(n_cell_types / 2)
  major_cell_types <- n_major_cell_types
  ## rare_cell_types <- n_cell_types - major_cell_types
  rare_cell_types <- n_minor_cell_types

  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "sample_id")

  # Generate baseline of cell type data.frame
  for (id in dummy_data$sample_id){
    # set.seed(seed*5*grep(id,dummy_data$sample_id)*10)
    print(paste("sample_id:", id, sep=" "))
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
    print(paste("major_cell_counts:", major_cell_counts, sep=" "))
    # set.seed(seed*6*grep(id,dummy_data$sample_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
    print(paste("rare_cell_counts:", rare_cell_counts, sep=" "))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
    print(paste("total cell_counts:", sum(major_cell_counts, rare_cell_counts), sep=" "))
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
      print(paste("n:", n, sep=" "))
      print(dim(celltype_df))
      celltype_df = rbind(celltype_df,
                          data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], n),
                                     sample_id = id)
      )
    }
  }

  # check
  ## celltype_df %>%
  ##   dplyr::group_by(sample_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  ## celltype_df %>%
  ##   dplyr::group_by(cell_type,sample_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()

  # --- LV ADDITION: generalized, visit-specific effects --------------------------------
  # all new code within the ---- marks.  Please revert back if necessary

  interact_clusters = c(1:n_major_interact_celltypes, (n_cell_types - n_minor_interact_celltypes + 1):n_cell_types)
  interact_cell_types = LETTERS[interact_clusters]

  # defaults if user did not provide visit-specific effects
  if (is.null(visit_effects_progressor)) {
    # Keep old behavior as a fallback: map interaction_type/fc_interact to visits
    # specific:     only middle visit gets bump when time_points==2; for 3+, bump V1 only
    # differential: all visits get same bump
    # opposite:     up then down, or vice versa (for 3+, up at V0, down at V1, up at V2)
    if (interaction_type[1] == "specific") {
      ve <- rep(0, time_points); if (time_points >= 2) ve[2] <- fc_interact
      visit_effects_progressor <- ve
    } else if (interaction_type[1] == "differential") {
      visit_effects_progressor <- rep(fc_interact, time_points)
    } else { # "opposite"
      visit_effects_progressor <- rep(0, time_points)
      visit_effects_progressor[seq(1, time_points, by = 2)] <- +fc_interact
      if (time_points >= 2) visit_effects_progressor[seq(2, time_points, by = 2)] <- -fc_interact
    }
  }
  if (is.null(visit_effects_control)) visit_effects_control <- rep(0, time_points)

  # direction per interacting cluster (+1 increase, -1 decrease)
  if (is.null(direction_by_cluster)) direction_by_cluster <- rep(1, length(interact_clusters))
  names(direction_by_cluster) <- LETTERS[interact_clusters]

  # Helper: get current per-sample counts for a given cell type
  get_counts <- function(ct) {
    celltype_df |>
      dplyr::filter(cell_type == ct) |>
      dplyr::count(sample_id, name = "count")
  }

  # Compute baseline median proportion per interacting cell type
  baseline_prop <- list()
  for (ct in LETTERS[1:n_cell_types]) {
    ab <- dplyr::left_join(
      get_counts(ct),
      dplyr::select(dummy_data, sample_id),
      by = "sample_id",
      copy = TRUE   # forces a local copy if needed
    )

    # counts were generated around n_cells; use that as denominator
    p <- ab$count / n_cells
    baseline_prop[[ct]] <- stats::median(p, na.rm = TRUE)
  }

  celltype_df <- tibble::as_tibble(celltype_df)
  dummy_data  <- tibble::as_tibble(dummy_data)

  # Function to add/remove rows to match a target count for each sample
  adjust_to_target <- function(ct, sample_ids, target_counts) {
    # force simple types
    sample_ids    <- as.character(sample_ids)
    target_counts <- as.integer(target_counts)
    stopifnot(length(sample_ids) == length(target_counts))

    # rows for this cell type and these samples
    rows_ct <- which(celltype_df$cell_type == ct & celltype_df$sample_id %in% sample_ids)

    # current counts per sample (default 0 when absent)
    tab <- table(celltype_df$sample_id[rows_ct])
    cur_counts <- integer(length(sample_ids))
    names(cur_counts) <- sample_ids
    if (length(tab)) cur_counts[names(tab)] <- as.integer(tab)

    # how many to add / remove per sample (never NA)
    add_needed <- pmax(0L, target_counts - cur_counts)
    rem_needed <- pmax(0L, cur_counts - target_counts)

    # ADD rows
    add_total <- sum(add_needed)
    if (add_total > 0L) {
      celltype_df <<- rbind(
        celltype_df,
        data.frame(
          cell_type = rep(ct, add_total),
          sample_id = rep(sample_ids, add_needed),
          stringsAsFactors = FALSE
        )
      )
    }

    # REMOVE rows
    rem_total <- sum(rem_needed)
    if (rem_total > 0L) {
      celltype_df <<- transform(celltype_df, row_id = seq_len(nrow(celltype_df)))
      rem_rows <- integer(0)

      # iterate by position (not by name); avoids NA indexing
      for (i in seq_along(sample_ids)) {
        need <- rem_needed[i]
        if (need > 0L) {
          sid  <- sample_ids[i]
          cand <- celltype_df$row_id[celltype_df$cell_type == ct & celltype_df$sample_id == sid]
          if (length(cand) > 0L) {
            take_n  <- min(need, length(cand))
            rem_rows <- c(rem_rows, head(cand, take_n))
          }
        }
      }

      if (length(rem_rows) > 0L) {
        celltype_df <<- subset(celltype_df, !(row_id %in% rem_rows), select = -row_id)
      } else {
        celltype_df <<- subset(celltype_df, select = -row_id)
      }
    }
  }

  # Now enforce target proportions per visit and disease group
  for (ct in interact_cell_types) {
    dir_ct <- direction_by_cluster[ct]
    base_p <- baseline_prop[[ct]]

    for (v in 0:(time_points - 1)) {
      # group sample_ids
      s_prog <- dummy_data$sample_id[dummy_data$disease == 1 & dummy_data$visit == v]
      s_ctrl <- dummy_data$sample_id[dummy_data$disease == 0 & dummy_data$visit == v]

      # target proportions
      tp_prog <- base_p * (1 + dir_ct * visit_effects_progressor[v + 1])
      tp_ctrl <- base_p * (1 +            visit_effects_control[v + 1])

      # translate to target counts per sample (cap at [0, n_cells])
      tgt_prog <- pmax(0, pmin(n_cells, round(n_cells * tp_prog)))
      tgt_ctrl <- pmax(0, pmin(n_cells, round(n_cells * tp_ctrl)))

      # adjust
      if (length(s_prog)) adjust_to_target(ct, s_prog, rep(tgt_prog, length(s_prog)))
      if (length(s_ctrl)) adjust_to_target(ct, s_ctrl, rep(tgt_ctrl, length(s_ctrl)))
    }
  }

  dummy_data = merge(dummy_data, celltype_df, by = "sample_id")
  # ------------------------------------------------------------------------------

  # Shuffle rows
  # set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(dummy_data)
}
