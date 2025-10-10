generate_dummy_data_time <- function(
    n_cells, sd_celltypes, n_major_cell_types, n_minor_cell_types,
    relative_abundance, n_major_interact_celltypes, n_minor_interact_celltypes,
    n_individuals, n_batchs, interaction_feature, time_points, test_var,
    prop_disease, fc_interact, interaction_type = c("specific","differential","opposite"),
    seed = 1234,
    visit_effects_progressor,    # length = time_points
    visit_effects_control,       # length = time_points
    direction_by_cluster         # +1/-1 per interacting cluster
) {
  if (time_points < 3) {
    generate_dummy_data_2time(
      n_cells = n_cells,
      sd_celltypes = sd_celltypes,
      n_major_cell_types = n_major_cell_types,
      n_minor_cell_types = n_minor_cell_types,
      relative_abundance = relative_abundance,
      n_major_interact_celltypes = n_major_interact_celltypes,
      n_minor_interact_celltypes = n_minor_interact_celltypes,
      n_individuals = n_individuals,
      n_batchs = n_batchs,
      interaction_feature = interaction_feature,
      time_points = time_points,
      test_var = test_var,
      prop_disease = prop_disease,
      fc_interact = fc_interact,
      interaction_type = interaction_type,
      seed = seed
    )
  } else {
    generate_dummy_data_time_3plus(
      n_cells = n_cells,
      sd_celltypes = sd_celltypes,
      n_major_cell_types = n_major_cell_types,
      n_minor_cell_types = n_minor_cell_types,
      relative_abundance = relative_abundance,
      n_major_interact_celltypes = n_major_interact_celltypes,
      n_minor_interact_celltypes = n_minor_interact_celltypes,
      n_individuals = n_individuals,
      n_batchs = n_batchs,
      interaction_feature = interaction_feature,
      time_points = time_points,
      test_var = test_var,
      prop_disease = prop_disease,
      fc_interact = fc_interact,
      interaction_type = interaction_type,
      seed = seed,
      visit_effects_progressor = visit_effects_progressor,
      visit_effects_control    = visit_effects_control,
      direction_by_cluster     = direction_by_cluster
    )
  }
}

