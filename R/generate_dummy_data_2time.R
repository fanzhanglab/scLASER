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
#'
#' @return
#' @export
#'
#' @examples
generate_dummy_data_2time <- function(n_cells = 3000, # cells of major cell types per individual
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
   # print(paste("sample_id:", id, sep=" "))
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
  #  print(paste("major_cell_counts:", major_cell_counts, sep=" "))
    # set.seed(seed*6*grep(id,dummy_data$sample_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
   # print(paste("rare_cell_counts:", rare_cell_counts, sep=" "))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
   # print(paste("total cell_counts:", sum(major_cell_counts, rare_cell_counts), sep=" "))
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
    #  print(paste("n:", n, sep=" "))
    #  print(dim(celltype_df))
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

  interact_clusters = c(1:n_major_interact_celltypes, (n_cell_types-n_minor_interact_celltypes+1):n_cell_types)
  interact_cell_types = LETTERS[interact_clusters]

  for (i in 1:n_cell_types) {
    if(i %in% interact_clusters){
      if(i %in% interact_clusters[seq(2,length(interact_clusters),2)]){

        abundance = dplyr::left_join(celltype_df,
                                     dummy_data,
                                     by="sample_id") %>%
          dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
          dplyr::group_by(sample_id) %>%
          dplyr::summarise(pro = dplyr::n()/n_cells,
                           count = dplyr::n())
        diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) - n_cells*median(abundance$pro))

        if(interaction_type == "specific"){
          # add diff cells only to both interaction category == 1

          len = length(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id)) * diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id),diff)),
            celltype_df
          )
        } else if (interaction_type == "differential"){
          # add diff cells both to interaction category 2 == 1
          len = length(unique(dummy_data[dummy_data$interact_term == 2 | dummy_data$interact_term == 1,]$sample_id)) * diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 2 | dummy_data$interact_term == 1,]$sample_id),diff)),
            celltype_df
          )
          # add new_diff cells only to both interaction category == 1 additionally
          new_abundance = dplyr::left_join(celltype_df,
                                           dummy_data,
                                           by="sample_id") %>%
            dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
            dplyr::filter(interact_term == 2 | interact_term == 1) %>%
            dplyr::group_by(sample_id) %>%
            dplyr::summarise(pro = dplyr::n()/n_cells,
                             count = dplyr::n())

          new_diff = ceiling(n_cells*(median(new_abundance$pro)*(1+fc_interact)) - n_cells*median(new_abundance$pro))

          len = length(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id)) * new_diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id),new_diff)),
            celltype_df
          )

        } else if (interaction_type == "opposite"){
          # add diff cells only to both interaction category == 1
          len = length(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id)) * diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id),diff)),
            celltype_df
          )

          # remove diff cells only to interaction category 1 == 1
          rem_rows = celltype_df %>%
            dplyr::mutate(row_id = dplyr::row_number()) %>%
            dplyr::left_join(dummy_data,
                             by="sample_id") %>%
            dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
            dplyr::filter(interact_term == 1) %>%
            dplyr::group_by(cell_type,sample_id) %>%
            dplyr::slice(1:diff) %>%
            .$row_id %>%
            as.integer()

          celltype_df = celltype_df %>%
            dplyr::mutate(row_id = dplyr::row_number()) %>%
            dplyr::filter(!(row_id %in% rem_rows)) %>%
            dplyr::select(-row_id)
        }
      } else if (i %in% interact_clusters[seq(1,length(interact_clusters),2)]){
        abundance = dplyr::left_join(celltype_df,
                                     dummy_data,
                                     by="sample_id") %>%
          dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
          dplyr::group_by(sample_id) %>%
          dplyr::summarise(pro = dplyr::n()/n_cells,
                           count = dplyr::n())
        diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) - n_cells*median(abundance$pro))

        if(interaction_type == "specific"){
          # add diff cells only to both interaction category == 1

          #len = length(unique(dummy_data[dummy_data$interact_term == 2,]$sample_id)) * diff  # JUAN CHANGED 8/11/2025
          len = length(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id)) * diff


          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id),diff)),
            celltype_df
          )
        } else if (interaction_type == "differential"){
          # add diff cells both to interaction category 2 == 1
          len = length(unique(dummy_data[dummy_data$interact_term == 2 | dummy_data$interact_term == 1,]$sample_id)) * diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 2 | dummy_data$interact_term == 1,]$sample_id),diff)),
            celltype_df
          )
          # add new_diff cells only to both interaction category == 1 additionally
          new_abundance = dplyr::left_join(celltype_df,
                                           dummy_data,
                                           by="sample_id") %>%
            dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
            dplyr::filter(interact_term == 2 | interact_term == 1) %>%
            dplyr::group_by(sample_id) %>%
            dplyr::summarise(pro = dplyr::n()/n_cells,
                             count = dplyr::n())

          new_diff = ceiling(n_cells*(median(new_abundance$pro)*(1+fc_interact)) - n_cells*median(new_abundance$pro))

          len = length(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id)) * new_diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id),new_diff)),
            celltype_df
          )

        } else if (interaction_type == "opposite"){
          # add diff cells only to both interaction category == 1
          len = length(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id)) * diff

          celltype_df = rbind(
            data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], len),
                       sample_id = rep(unique(dummy_data[dummy_data$interact_term == 1,]$sample_id),diff)),
            celltype_df
          )

          # remove diff cells only to interaction category 1 == 1
          rem_rows = celltype_df %>%
            dplyr::mutate(row_id = dplyr::row_number()) %>%
            dplyr::left_join(dummy_data,
                             by="sample_id") %>%
            dplyr::filter(cell_type == LETTERS[seq( from = 1, to = n_cell_types )][i]) %>%
            dplyr::filter(interact_term == 2) %>%
            dplyr::group_by(cell_type,sample_id) %>%
            dplyr::slice(1:diff) %>%
            .$row_id %>%
            as.integer()

          celltype_df = celltype_df %>%
            dplyr::mutate(row_id = dplyr::row_number()) %>%
            dplyr::filter(!(row_id %in% rem_rows)) %>%
            dplyr::select(-row_id)

        }
      }


    }
  }

  dummy_data = merge(dummy_data,celltype_df,by="sample_id")

  # Shuffle rows
  # set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,interact_cell_types))
}
