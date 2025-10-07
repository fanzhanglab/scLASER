generate_dummy_data_time <- function(n_cells = 3000, # cells of major cell types per individual
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

){
  if (time_points < 3){
    generate_dummy_data_2time(n_cells = 3000, # cells of major cell types per individual
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
                              seed = 1234
    )
  } else {generate_dummy_data_time_3plus(n_cells = 3000, # cells of major cell types per individual
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

  )}
}
