# Required libraries

# remotes::install_version("Matrix", version = "1.6-3")
suppressWarnings({
  library(dplyr)
  library(tidyr)
 # library(tidyverse)
  library(MASS) # Used for sampling from normal distribution
  library(caret)
  
  library(Seurat)
  library(glue)
  
  library(stevemisc)
  library(stevedata)
  library(lme4)
  library(broom.mixed)
  
  # With parallelism
  library(doParallel)
  library(foreach)
  library(pbapply)
  
  library(variancePartition)
  
  library(Seurat)
})

# Functions
## generate meta data
generate_dummy_metadata <- function(n_cells = 3000, # cells of major cell types per individual
                                    sd_celltypes = 0.1,  # standard deviation of number of cells 
                                    n_major_cell_types = 7,
                                    n_minor_cell_types = 3,
                                    relative_abundance = 0.1, # ratio between major and rare
                                    n_major_diff_celltypes = 1,
                                    n_minor_diff_celltypes = 1,
                                    n_individuals = 30, # total individuals
                                    n_batchs = 4,
                                    prop_sex = 0.5, 
                                    prop_disease = 0.5,  
                                    fc_increase = 0.1, # additional FC of specified cell types which are from people with case groups compared to control groups (e.g., if you specify 0.5, FC comparing increased cell type between case and control groups will be 1.5 in total)
                                    seed = 1234
) {
  n_cell_types = n_major_cell_types + n_minor_cell_types
  
  # Generate unique subject IDs
  subject_id <- paste0("SUB_", 1:n_individuals)
  set.seed(seed)
  sex = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_sex)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_sex)))
  )
  if (length(sex) != length(subject_id)) {
    sex = c(sex, rep(1,length(subject_id)-length(sex)))
  }
  set.seed(seed*2)
  disease = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_disease)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_disease))
    )
  )
  if (length(disease) != length(subject_id)) {
    disease = c(disease, rep(1,length(subject_id)-length(disease)))
  }
  set.seed(seed*3)
  age <- sample(18:60, n_individuals, replace = TRUE)
  set.seed(seed*4)
  bmi <- sample(15:35, n_individuals, replace = TRUE)
  batch <- rep(rep(1:n_batchs), length(subject_id))[1:length(subject_id)]
  dummy_data = data.frame(subject_id = subject_id,
                          sex = factor(sex, levels = c(0, 1)),
                          disease = factor(disease, levels = c(0, 1)),
                          age = age,
                          batch = factor(batch, levels = c(0:max(batch))),
                          bmi = bmi
  ) 
  
  # Major and rare cell type counts
  ## major_cell_types <- ceiling(n_cell_types / 2)
  major_cell_types <- n_major_cell_types
  ## rare_cell_types <- n_cell_types - major_cell_types
  rare_cell_types <- n_minor_cell_types
  
  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "subject_id")
  
  # Generate baseline of cell type data.frame 
  for (id in dummy_data$subject_id){
    set.seed(seed*5*grep(id,dummy_data$subject_id)*10)
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
    set.seed(seed*6*grep(id,dummy_data$subject_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
      celltype_df = rbind(celltype_df,
                          data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], n),
                                     subject_id = id)
      )
    }
  }
  
  # check
  ## celltype_df %>%
  ##   dplyr::group_by(subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  ## celltype_df %>%
  ##   dplyr::group_by(cell_type,subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  
  diff_clusters = c(1:n_major_diff_celltypes, (n_cell_types-n_minor_diff_celltypes+1):n_cell_types)
  diff_cell_types = LETTERS[diff_clusters]
  
  for (i in LETTERS[1:n_cell_types]) {
    if(i %in% diff_cell_types){
      
      abundance = dplyr::left_join(celltype_df,
                                   dummy_data,
                                   by="subject_id") %>%
        dplyr::filter(cell_type == i) %>%
        dplyr::group_by(subject_id) %>%
        dplyr::summarise(pro = dplyr::n()/n_cells,
                         count = dplyr::n())
      diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_increase)) - n_cells*median(abundance$pro))
      
      i# add diff cells only to disease == 1
      
      len = length(unique(dummy_data[dummy_data$disease == "1",]$subject_id)) * diff
      
      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "1",]$subject_id),diff)),
        celltype_df
      )
    }
  }
  dummy_data = merge(dummy_data,celltype_df,by="subject_id")
  # Shuffle rows
  set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,diff_cell_types))
}
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

generate_pseudo_pcs_time = function(data, 
                                    
                                    n_pcs = 20, 
                                    
                                    cluster_pcs = 1:20,
                                    disease_pcs = 0,
                                    sex_pcs = 0,
                                    age_pcs = 0,
                                    bmi_pcs = 0,
                                    batch_pcs = 0,
                                    interaction_pcs = 0,
                                    visit_pcs = 0,
                                    subject_pcs = 0,
                                    
                                    scale_factor = 2,
                                    
                                    cluster_ratio = 0.25,
                                    disease_ratio = 0,
                                    sex_ratio = 0,
                                    age_ratio = 0,
                                    bmi_ratio = 0,
                                    batch_ratio = 0,
                                    visit_ratio = 0,
                                    interaction_ratio = 0,
                                    subject_ratio = 0,
                                    
                                    cluster_col = "cell_type",
                                    disease_col = "disease",
                                    sex_col = "sex",
                                    age_col = "age",
                                    bmi_col = "bmi",
                                    batch_col = "batch",
                                    visit_col = "visit",
                                    interact_term_col = "interact_term",
                                    seed = 1234
) {
  # set.seed(seed)
  disease_pcs = sample(disease_pcs)
  # set.seed(seed*2)
  sex_pcs = sample(sex_pcs)
  # set.seed(seed*3)
  age_pcs = sample(age_pcs)
  # set.seed(seed*4)
  bmi_pcs = sample(bmi_pcs)
  # set.seed(seed*5)
  batch_pcs = sample(batch_pcs)
  
  sapply(1:n_pcs, function(x){
    
    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluster_col]))
    n_pcs <- n_pcs
    
    cell_clusters <- as.integer(factor(data[,cluster_col])) 
    for (j in 1:3){
      cell_clusters_tmp = cell_clusters
      # set.seed(x*j)
      shufle = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(shufle[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters",j),cell_clusters_tmp)
    }
    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    cell_interaction <- as.integer(factor(data[,interact_term_col]))
    
    variance <- 1 / (x * scale_factor)  
    # set.seed(seed*x)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))
    # set.seed(seed*x)
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters1), sd = sqrt(variance))
    # for (j in 1:3){
    #   if (j==1){
    #     pc_cluster = rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance))
    #   } else {
    #     pc_cluster = pc_cluster + rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance)) 
    #   }
    # }
    # set.seed(seed*x)
    pc_disease = rnorm(n_cells, mean = scale(cell_diseases), sd = sqrt(variance))
    # set.seed(seed*x*2)
    pc_sex = rnorm(n_cells, mean = scale(cell_sex), sd = sqrt(variance))
    # set.seed(seed*x*3)
    pc_age = rnorm(n_cells, mean = scale(cell_age), sd = sqrt(variance))
    # set.seed(seed*x*4)
    pc_bmi = rnorm(n_cells, mean = scale(cell_bmi), sd = sqrt(variance))
    # set.seed(seed*x*5)
    pc_batch = rnorm(n_cells, mean = scale(cell_batch), sd = sqrt(variance))
    # set.seed(seed*x*6)
    pc_interact = rnorm(n_cells, mean = scale(cell_interaction), sd = sqrt(variance))
    # set.seed(seed*x*7)
    pc_visit = rnorm(n_cells, mean = scale(data[,visit_col]), sd = sqrt(variance))
    # cluster_ratio_tmp = cluster_ratio / sqrt(charmatch(x,cluster_pcs))
    cluster_ratio_tmp = cluster_ratio / charmatch(x,cluster_pcs)
    # disease_ratio_tmp = disease_ratio / sqrt(charmatch(x,disease_pcs))
    disease_ratio_tmp = disease_ratio / charmatch(x,disease_pcs)
    # sex_ratio_tmp = sex_ratio / sqrt(charmatch(x,sex_pcs))
    sex_ratio_tmp = sex_ratio / charmatch(x,sex_pcs)
    # age_ratio_tmp = age_ratio / sqrt(charmatch(x,age_pcs))
    age_ratio_tmp = age_ratio / charmatch(x,age_pcs)
    # bmi_ratio_tmp = bmi_ratio / sqrt(charmatch(x,bmi_pcs))
    bmi_ratio_tmp = bmi_ratio / charmatch(x,bmi_pcs)
    # batch_ratio_tmp = batch_ratio / sqrt(charmatch(x,batch_pcs))
    batch_ratio_tmp = batch_ratio / charmatch(x,batch_pcs)
    # interaction_ratio_tmp = interaction_ratio / sqrt(charmatch(x,interaction_pcs))
    interaction_ratio_tmp = interaction_ratio / charmatch(x,interaction_pcs)
    # visit_ratio_tmp = visit_ratio / sqrt(charmatch(x,visit_pcs))
    visit_ratio_tmp = visit_ratio / charmatch(x,visit_pcs)
    if (!(x %in% cluster_pcs)) cluster_ratio_tmp = 0
    if (!(x %in% disease_pcs)) disease_ratio_tmp = 0
    if (!(x %in% sex_pcs)) sex_ratio_tmp = 0
    if (!(x %in% age_pcs)) age_ratio_tmp = 0
    if (!(x %in% bmi_pcs)) bmi_ratio_tmp = 0
    if (!(x %in% batch_pcs)) batch_ratio_tmp = 0
    if (!(x %in% interaction_pcs)) interaction_ratio_tmp = 0
    if (!(x %in% visit_pcs)) visit_ratio_tmp = 0
    return(pc * (1
                 -cluster_ratio_tmp
                 -disease_ratio_tmp
                 -sex_ratio_tmp
                 -age_ratio_tmp
                 -bmi_ratio_tmp
                 -batch_ratio_tmp
                 -interaction_ratio_tmp
                 -visit_ratio_tmp
    ) + 
      pc_cluster * cluster_ratio_tmp +
      pc_disease * disease_ratio_tmp +
      pc_sex * sex_ratio_tmp +
      pc_age * sex_ratio_tmp +
      pc_bmi * sex_ratio_tmp +
      pc_batch * sex_ratio_tmp +
      pc_interact * interaction_ratio_tmp +
      pc_visit * visit_ratio_tmp
    )
  })
}

generate_dummy_data_woInteraction <- function(n_cells = 3000, # cells of major cell types per individual
                                              sd_celltypes = 0.1,  # standard deviation of number of cells 
                                              n_major_cell_types = 7,
                                              n_minor_cell_types = 3,
                                              relative_abundance = 0.1, # ratio between major and rare
                                              n_major_diff_celltypes = 1,
                                              n_minor_diff_celltypes = 1,
                                              n_individuals = 30, # total individuals
                                              n_batchs = 4,
                                              prop_sex = 0.5, 
                                              prop_disease = 0.5,  
                                              fc_interact = 0.1, # additional proportion of interacted cells which are from people with interacted group compared to non-interacted cell types
                                              seed = 1234
) {
  n_cell_types = n_major_cell_types + n_minor_cell_types
  
  # Generate unique subject IDs
  subject_id <- paste0("SUB_", 1:n_individuals)
  set.seed(seed)
  sex = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_sex)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_sex)))
  )
  if (length(sex) != length(subject_id)) {
    sex = c(sex, rep(1,length(subject_id)-length(sex)))
  }
  set.seed(seed*2)
  disease = sample(
    c(rep(1, round(length(unique(subject_id)) * prop_disease)),
      rep(0, n_individuals - round(length(unique(subject_id)) * prop_disease))
    )
  )
  if (length(disease) != length(subject_id)) {
    disease = c(disease, rep(1,length(subject_id)-length(disease)))
  }
  set.seed(seed*3)
  age <- sample(18:60, n_individuals, replace = TRUE)
  set.seed(seed*4)
  bmi <- sample(15:35, n_individuals, replace = TRUE)
  batch <- rep(rep(1:n_batchs), length(subject_id))[1:length(subject_id)]
  dummy_data = data.frame(subject_id = subject_id,
                          sex = factor(sex, levels = c(0, 1)),
                          disease = factor(disease, levels = c(0, 1)),
                          age = age,
                          batch = factor(batch, levels = c(0:max(batch))),
                          bmi = bmi
  ) 
  
  # Major and rare cell type counts
  ## major_cell_types <- ceiling(n_cell_types / 2)
  major_cell_types <- n_major_cell_types
  ## rare_cell_types <- n_cell_types - major_cell_types
  rare_cell_types <- n_minor_cell_types
  
  celltype_df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(celltype_df) <- c("cell_type", "subject_id")
  
  # Generate baseline of cell type data.frame 
  for (id in dummy_data$subject_id){
    set.seed(seed*5*grep(id,dummy_data$subject_id)*10)
    major_cell_counts <- round(runif(major_cell_types, n_cells-n_cells*sd_celltypes, n_cells+n_cells*sd_celltypes))
    set.seed(seed*6*grep(id,dummy_data$subject_id)*10)
    rare_cell_counts <- round(runif(rare_cell_types, n_cells*relative_abundance-n_cells*relative_abundance*sd_celltypes, n_cells*relative_abundance+n_cells*relative_abundance*sd_celltypes))
    cell_counts <- c(major_cell_counts, rare_cell_counts)
    for (i in 1:n_cell_types) {
      n <- cell_counts[i]
      celltype_df = rbind(celltype_df,
                          data.frame(cell_type = rep(LETTERS[seq( from = 1, to = n_cell_types )][i], n),
                                     subject_id = id)
      )
    }
  }
  
  # check
  ## celltype_df %>%
  ##   dplyr::group_by(subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  ## celltype_df %>%
  ##   dplyr::group_by(cell_type,subject_id) %>%
  ##   dplyr::summarise(count = dplyr::n()) %>%
  ##   as.data.frame()
  
  diff_clusters = c(1:n_major_diff_celltypes, (n_cell_types-n_minor_diff_celltypes+1):n_cell_types)
  diff_cell_types = LETTERS[diff_clusters]
  prop_control_types <- 0.1
  
  for (i in LETTERS[1:n_cell_types]) {
    abundance = dplyr::left_join(celltype_df,
                                 dummy_data,
                                 by="subject_id") %>%
      dplyr::filter(cell_type == i) %>%
      dplyr::group_by(subject_id) %>%
      dplyr::summarise(pro = dplyr::n()/n_cells,
                       count = dplyr::n())
    if(i %in% diff_cell_types){
      
      # abundance = dplyr::left_join(celltype_df,
      #                              dummy_data,
      #                              by="subject_id") %>%
      #   dplyr::filter(cell_type == i) %>%
      #   dplyr::group_by(subject_id) %>%
      #   dplyr::summarise(pro = dplyr::n()/n_cells,
      #                    count = dplyr::n())
      
      # prop(ortion) = cells for each subject that are this cell type
      print(head(abundance))
      # diff is the number of cells that the differential cell types should have (use fc_interact to scale/increase)
      diff = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) - n_cells*median(abundance$pro))
      print(paste("Inside diff cell types. diff = ", diff))
      print(paste("n_cells = ", n_cells, "; median pro = ", median(abundance$pro), " ; median pro increase = ", (median(abundance$pro)*(1+fc_interact))))
      print(paste("i = ", i))
      i# add diff cells only to disease == 1
      
      len = length(unique(dummy_data[dummy_data$disease == "1",]$subject_id)) * diff
      #print(paste("unique(dummy_data[dummy_data$disease == '1',]$subject_id) : ", unique(dummy_data[dummy_data$disease == "1",]$subject_id)))
      # unique... is the number of subjects with disease == 1?
      print(paste("len = ", len))
      
      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "1",]$subject_id),diff)),
        celltype_df
      )
    }
    else{
      # number of disease cells in control cells
      #diff_control = ceiling(n_cells*(median(abundance$pro)*(prop_control_types))+ n_cells*median(abundance$pro))
      diff_control = ceiling(n_cells*(median(abundance$pro)*(1+fc_interact)) + n_cells*median(abundance$pro))
      print(paste("Inside control cell types. diff_control = ", diff_control))
      print(paste("n_cells = ", n_cells, "; median pro = ", median(abundance$pro), " ; median pro control = ", (median(abundance$pro)*(prop_control_types))))
      len = length(unique(dummy_data[dummy_data$disease == "0",]$subject_id)) * diff_control
      print(paste("len = ", len))
      celltype_df = rbind(
        data.frame(cell_type = rep(i, len),
                   subject_id = rep(unique(dummy_data[dummy_data$disease == "0",]$subject_id),diff_control)),
        celltype_df
      )
      
    }
  }
  dummy_data = merge(dummy_data,celltype_df,by="subject_id")
  # initial attempt to increase number of 0 (control) cells in non disease cell types
  # Shuffle rows
  set.seed(seed*7)
  dummy_data <- dummy_data[sample(nrow(dummy_data)), ]
  return(list(dummy_data,diff_cell_types))
}

generate_pseudo_pcs_woInteraction = function(data, 
                                             # number of principal components
                                             n_pcs = 20, 
                                             # Among the principal components representing each attribute, arrange the components in descending order of variance.
                                             cluster_pcs = 1:20,
                                             disease_pcs = 0,
                                             sex_pcs = 0,
                                             age_pcs = 0,
                                             bmi_pcs = 0,
                                             batch_pcs = 0,
                                             
                                             scale_factor = 2,
                                             # Define the maximum percentage of each attribute
                                             cluster_ratio = 0.25,
                                             disease_ratio = 0,
                                             sex_ratio = 0,
                                             age_ratio = 0,
                                             bmi_ratio = 0,
                                             batch_ratio = 0,
                                             
                                             cluser_col = "cell_type",
                                             disease_col = "disease",
                                             sex_col = "sex",
                                             age_col = "age",
                                             bmi_col = "bmi",
                                             batch_col = "batch",
                                             seed = 1234
) {
  set.seed(seed)
  disease_pcs = sample(disease_pcs)
  set.seed(seed*2)
  sex_pcs = sample(sex_pcs)
  set.seed(seed*3)
  age_pcs = sample(age_pcs)
  set.seed(seed*4)
  bmi_pcs = sample(bmi_pcs)
  set.seed(seed*5)
  batch_pcs = sample(batch_pcs)
  
  sapply(1:n_pcs, function(x){
    
    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluser_col]))
    n_pcs <- n_pcs
    
    cell_clusters <- as.integer(factor(data[,cluser_col])) 
    for (j in 1:3){
      cell_clusters_tmp = cell_clusters
      set.seed(x*j)
      shufle = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(shufle[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters",j),cell_clusters_tmp)
    }
    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    
    variance <- 1 / (x * scale_factor)  # Set so that the variance decreases as the principal component goes from 1 to 20.
    set.seed(seed*x)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))
    set.seed(seed*x)
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters1), sd = sqrt(variance))
    # for (j in 1:3){
    #   if (j==1){
    #     pc_cluster = rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance))
    #   } else {
    #     pc_cluster = pc_cluster + rnorm(n_cells, mean = scale(eval(parse(text=paste0("cell_clusters",j)))), sd = sqrt(variance)) 
    #   }
    # }
    set.seed(seed*x)
    pc_disease = rnorm(n_cells, mean = scale(cell_diseases), sd = sqrt(variance))
    set.seed(seed*x*2)
    pc_sex = rnorm(n_cells, mean = scale(cell_sex), sd = sqrt(variance))
    set.seed(seed*x*3)
    pc_age = rnorm(n_cells, mean = scale(cell_age), sd = sqrt(variance))
    set.seed(seed*x*4)
    pc_bmi = rnorm(n_cells, mean = scale(cell_bmi), sd = sqrt(variance))
    set.seed(seed*x*5)
    pc_batch = rnorm(n_cells, mean = scale(cell_batch), sd = sqrt(variance))
    # cluster_ratio_tmp = cluster_ratio / sqrt(charmatch(x,cluster_pcs))
    cluster_ratio_tmp = cluster_ratio / charmatch(x,cluster_pcs)
    # disease_ratio_tmp = disease_ratio / sqrt(charmatch(x,disease_pcs))
    disease_ratio_tmp = disease_ratio / charmatch(x,disease_pcs)
    # sex_ratio_tmp = sex_ratio / sqrt(charmatch(x,sex_pcs))
    sex_ratio_tmp = sex_ratio / charmatch(x,sex_pcs)
    # age_ratio_tmp = age_ratio / sqrt(charmatch(x,age_pcs))
    age_ratio_tmp = age_ratio / charmatch(x,age_pcs)
    # bmi_ratio_tmp = bmi_ratio / sqrt(charmatch(x,bmi_pcs))
    bmi_ratio_tmp = bmi_ratio / charmatch(x,bmi_pcs)
    # batch_ratio_tmp = batch_ratio / sqrt(charmatch(x,batch_pcs))
    batch_ratio_tmp = batch_ratio / charmatch(x,batch_pcs)
    if (!(x %in% cluster_pcs)) cluster_ratio_tmp = 0
    if (!(x %in% disease_pcs)) disease_ratio_tmp = 0
    if (!(x %in% sex_pcs)) sex_ratio_tmp = 0
    if (!(x %in% age_pcs)) age_ratio_tmp = 0
    if (!(x %in% bmi_pcs)) bmi_ratio_tmp = 0
    if (!(x %in% batch_pcs)) batch_ratio_tmp = 0
    return(pc * (1
                 -cluster_ratio_tmp
                 -disease_ratio_tmp
                 -sex_ratio_tmp
                 -age_ratio_tmp
                 -bmi_ratio_tmp
                 -batch_ratio_tmp
    ) + 
      pc_cluster * cluster_ratio_tmp +
      pc_disease * disease_ratio_tmp +
      pc_sex * sex_ratio_tmp +
      pc_age * sex_ratio_tmp +
      pc_bmi * sex_ratio_tmp +
      pc_batch * sex_ratio_tmp
    )
  })
}

generate_pseudo_features = function(data,
                                    n_features = 20,
                                    cluster_features = 1:20,
                                    disease_features = 0,
                                    sex_features = 0,
                                    age_features = 0,
                                    bmi_features = 0,
                                    batch_features = 0,
                                    individual_features = 0,
                                    
                                    cluster_ratio = 0.25,
                                    disease_ratio = 0,
                                    sex_ratio = 0,
                                    age_ratio = 0,
                                    bmi_ratio = 0,
                                    batch_ratio = 0,
                                    individual_ratio = 0.1,
                                    ratio_variance = 0.5,
                                    
                                    cluster_col = "cell_type",
                                    disease_col = "disease",
                                    sex_col = "sex",
                                    age_col = "age",
                                    bmi_col = "bmi",
                                    batch_col = "batch",
                                    individual_col = "batch",
                                    seed = 1234,
                                    n_thread = 1,
                                    verbose = TRUE
) {
  
  library(doParallel)
  library(foreach)
  library(dplyr)
  
  set.seed(seed)
  disease_features = sample(disease_features)
  set.seed(seed*2)
  sex_features = sample(sex_features)
  set.seed(seed*3)
  age_features = sample(age_features)
  set.seed(seed*4)
  bmi_features = sample(bmi_features)
  set.seed(seed*5)
  batch_features = sample(batch_features)
  set.seed(seed*6)
  individual_features = sample(individual_features)
  
  # Register parallel backend to use multiple cores
  cl <- makeCluster(n_thread) # Leave one core free for other processes
  registerDoParallel(cl)
  
  pcs_list = foreach(x = 1:n_features, .packages = c("dplyr")) %dopar% {
    
    n_cells <- nrow(data)
    n_clusters <- length(unique(data[,cluster_col]))
    
    cell_sex <- as.integer(factor(data[,sex_col]))
    cell_age <- as.integer(factor(data[,age_col]))
    cell_bmi <- as.integer(factor(data[,bmi_col]))
    cell_batch <- as.integer(factor(data[,batch_col]))
    cell_diseases <- as.integer(factor(data[,disease_col]))
    cell_individual <- as.integer(factor(data[,individual_col]))
    
    var_all <- c()
    
    # Generate dummy features reflecting cell clusters
    ## CELL CLUSTERS
    cell_clusters <- as.integer(factor(data[,cluster_col]))
    for (j in 1:2){
      cell_clusters_tmp = cell_clusters
      set.seed(x*j)
      cluster_mean = sample(unique(cell_clusters_tmp))*10
      for (i in 1:length(unique(cell_clusters))){
        cell_clusters_tmp <- dplyr::case_when(
          cell_clusters_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_clusters_tmp
        )
      }
      assign(paste0("cell_clusters_means",j),cell_clusters_tmp)
      print(paste("length of cell_clusters_tmp: ", length(cell_clusters_tmp)))
    }
    # paste0("cell_clusters_means",j) include cluster-specific mean values
    ## length(cell_clusters_means1)==length(cell_clusters)
    ## length(cell_clusters_means2)==length(cell_clusters)
    # check if cluster means are same value by cluster
    ## cell_clusters_means1[cell_clusters==1]
    ## cell_clusters_means2[cell_clusters==1]
    ## cell_clusters_means1[cell_clusters==2]
    ## cell_clusters_means2[cell_clusters==2]
    variance <- 1 / cell_clusters_means2 # cell type specific variance
    set.seed(seed*x)
    ### =====
    # The following R code generates dummy PC values by cells from a normal distribution by means of random numbers.
    # Specifically, the rnorm function is used to generate n_cells of random numbers, and the mean of each PC value for each cell cluster is the (scaled) unique value (cell_clusters_means1) for each cell type. The square root of the specific variance is the standard deviation.
    #1. rnorm(n) generates n random values from a standard normal distribution (mean 0, standard deviation 1).
    #2. mean = scale(cell_clusters1) specifies the mean of the random values to be generated. Here, the scale() function standardizes the data (subtracts the mean and divides by the standard deviation), but the use of the scale function may be inappropriate in this context. Normally, the scale() function standardizes a column of vectors or data frames and converts them to a value with mean 0 and standard deviation 1. However, if the goal here is to set the mean, then the use of mean(cell_clusters1) is appropriate.
    #3. sd = sqrt(variance) specifies the standard deviation of the random value to be generated. variance is the value of its variance, and sqrt(variance) gives the standard deviation.
    ### =====
    pc_cluster = rnorm(n_cells, mean = scale(cell_clusters_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    
    ## DISEASE
    # Similary, generate dummy PC for other covariates with changing seeds
    cell_covariate = cell_diseases
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*2)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*2)
    pc_disease = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    # check if PC values are different by disease group
    ## data.frame(cell_covariate,pc_disease) %>%
    ##   dplyr::group_by(cell_covariate) %>%
    ##   dplyr::summarize(mean_pc = mean(pc_disease),
    ##                    sd_pc = sd(pc_disease))
    var_all = c(var_all, variance)
    
    ## SEX
    cell_covariate = cell_sex
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*3)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*3)
    pc_sex = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    
    # Treat age as fixed effect for PC mean
    ## AGE
    cell_covariate = cell_age
    for (j in 2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*4)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*4)
    pc_age = rnorm(n_cells, mean = scale(data[,age_col]), sd = sqrt(variance))
    ## data.frame(data[,age_col],pc_age) %>%
    ##   magrittr::set_colnames(c("cell_covariate","pc_age")) %>%
    ##   dplyr::group_by(cell_covariate) %>%
    ##     dplyr::summarize(mean_pc = mean(pc_age),
    ##                      sd_pc = sd(pc_age)) %>%
    ##   as.data.frame()
    var_all = c(var_all, variance)
    ## BMI
    cell_covariate = cell_bmi
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*5)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*5)
    pc_bmi = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    ## BATCH
    cell_covariate = cell_batch
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*6)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2 # cell type specific variance
    set.seed(seed*x*6)
    pc_batch = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    data.frame(cell_covariate,pc_batch) %>%
      dplyr::group_by(cell_covariate) %>%
      dplyr::summarize(mean_pc = mean(pc_batch),
                       sd_pc = sd(pc_batch))
    # SUBJECT ID?
    cell_covariate = cell_individual
    for (j in 1:2){
      cell_covariates_tmp = cell_covariate
      set.seed(x*j*7)
      cluster_mean = sample(unique(cell_covariates_tmp))*10
      for (i in 1:length(unique(cell_covariate))){
        cell_covariates_tmp <- dplyr::case_when(
          cell_covariates_tmp == i ~ as.integer(cluster_mean[i]),
          TRUE ~ cell_covariates_tmp
        )
      }
      assign(paste0("cell_covariates_means",j),cell_covariates_tmp)
    }
    variance <- 1 / cell_covariates_means2
    set.seed(seed*x*7)
    pc_individual = rnorm(n_cells, mean = scale(cell_covariates_means1), sd = sqrt(variance))
    var_all = c(var_all, variance)
    data.frame(cell_covariate,pc_individual) %>%
      dplyr::group_by(cell_covariate) %>%
      dplyr::summarize(mean_pc = mean(pc_individual),
                       sd_pc = sd(pc_individual))
    
    # Generate dummy PC regardless cell types or other covariates (noise term)
    set.seed(seed*x*100)
    variance = sample(var_all,n_cells)
    set.seed(seed*x*100)
    pc = rnorm(n_cells, mean = 0, sd = sqrt(variance))
    
    cluster_ratio_tmp = cluster_ratio + runif(1, min = -cluster_ratio*ratio_variance, max = cluster_ratio*ratio_variance)
    disease_ratio_tmp = disease_ratio + runif(1, min = -disease_ratio*ratio_variance, max = disease_ratio*ratio_variance)
    sex_ratio_tmp = sex_ratio + runif(1, min = -sex_ratio*ratio_variance, max = sex_ratio*ratio_variance)
    age_ratio_tmp = age_ratio + runif(1, min = -age_ratio*ratio_variance, max = age_ratio*ratio_variance)
    bmi_ratio_tmp = bmi_ratio + runif(1, min = -bmi_ratio*ratio_variance, max = bmi_ratio*ratio_variance)
    batch_ratio_tmp = batch_ratio  + runif(1, min = -batch_ratio*ratio_variance, max = batch_ratio*ratio_variance)
    individual_ratio_tmp = individual_ratio  + runif(1, min = -individual_ratio*ratio_variance, max = individual_ratio*ratio_variance)
    if (!(x %in% cluster_features)) cluster_ratio_tmp = 0
    if (!(x %in% disease_features)) disease_ratio_tmp = 0
    if (!(x %in% sex_features)) sex_ratio_tmp = 0
    if (!(x %in% age_features)) age_ratio_tmp = 0
    if (!(x %in% bmi_features)) bmi_ratio_tmp = 0
    if (!(x %in% batch_features)) batch_ratio_tmp = 0
    if (!(x %in% individual_features)) individual_ratio_tmp = 0
    
    noise_ratio_tmp = (1
                       -cluster_ratio_tmp
                       -disease_ratio_tmp
                       -sex_ratio_tmp
                       -age_ratio_tmp
                       -bmi_ratio_tmp
                       -batch_ratio_tmp
                       -individual_ratio_tmp
    )
    if (noise_ratio_tmp < 0) noise_ratio_tmp = 0
    
    if(verbose){
      message(paste0("Feature",x,";\ncluster ratio = ",cluster_ratio_tmp,
                     "\ndisease ratio = ", disease_ratio_tmp,
                     "\nsex ratio = ", sex_ratio_tmp,
                     "\nage ratio = ", age_ratio_tmp,
                     "\nbmi ratio = ", bmi_ratio_tmp,
                     "\nbatch ratio = ", batch_ratio_tmp,
                     "\nindividual ratio = ", individual_ratio_tmp,
                     "\nnoise ratio = ", noise_ratio_tmp))
    }
    return(scale(pc) * noise_ratio_tmp +
             scale(pc_cluster) * cluster_ratio_tmp  +
             scale(pc_disease) * disease_ratio_tmp  +
             scale(pc_sex) * sex_ratio_tmp  +
             scale(pc_age) * age_ratio_tmp  +
             scale(pc_bmi) * bmi_ratio_tmp  +
             scale(pc_batch) * batch_ratio_tmp +
             scale(pc_individual) * individual_ratio_tmp
    )
  }
  
  stopCluster(cl)
  
  pcs <- do.call(cbind, pcs_list)
  return(pcs)
}



conditional_permutation <- function(B, Y, num) {
  purrr::map(seq_len(num), function(i) {
    split(seq_len(length(Y)), B) %>% purrr::map(function(idx) {
      data.frame(idx, val=sample(Y[idx]))
    }) %>% dplyr::bind_rows() %>% 
      dplyr::arrange(idx) %>% 
      with(val)    
  }) %>% 
    purrr::reduce(Matrix::cbind2)
}

empirical_fdrs <- function(z, znull, thresholds) {    
  n <- length(thresholds) - 1
  tails <- t(tail_counts(thresholds, znull)[1:n, ])
  ranks <- t(tail_counts(thresholds, z)[1:n, ])
  
  # compute FDPs
  fdp <- sweep(tails, 2, ranks, '/')
  fdr <- Matrix::colMeans(fdp)
  
  return(fdr)
}

tail_counts <- function(z, znull) {
  apply(znull, 2, function(znulli) {
    as.numeric(length(znulli) - cumsum(table(cut(znulli**2, c(0, z**2)))))        
  })
}

association_nam = function(seurat_object=NULL,
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
                           max_frac_pcs=0.15, # added option to pass number of PCs, for potential speedups
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
    # message(paste0("Running graph with FindNeighbors()..."))
    meta = metadata
    rownames(meta) = 1:nrow(meta)
    
    m <- as(t(pcs), "dgTMatrix") # by default, Matrix() returns dgCMatrix
    colnames(m) = 1:ncol(m)
    
    obj <- Seurat::CreateSeuratObject(
      counts = m, ## Subset expression matrix to cells in metadata
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
  
  ## (1) format data 
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
  
  # formatting and error checking
  ## For association, batches needs to be a numeric vector
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
  s <- s[, data$samplem[[data$samplem_key]]] ## Necessary? 
  # s = Matrix indicating which sample the cells that represent a row are from
  ## row: cell
  ## column: individual
  ## 1=the cell from that sample
  ### prior knowledge of connectivity among cells
  
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
  
  snorm <- t(prop.table(s, 2)) # normalization
  rownames(snorm) <- data$samplem[[data$samplem_key]]
  colnames(snorm) <- data$obs[[data$obs_key]]
  
  NAM=snorm
  
  N <- nrow(NAM)
  ## NOTE: added NULL check 
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
  # length(kurtoses) == nrow(NAM)
  threshold <- max(6, 2*median(kurtoses)) 
  message(glue::glue('throwing out neighborhoods with batch kurtosis >= {threshold}')) 
  keep <- which(kurtoses < threshold) 
  
  # keep <- rep(TRUE, ncol(NAM)) 
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
    ## NOTE: scale function in R gives slightly different results
    ##       than does similar function in python
    NAM_=scale(NAM_, center=FALSE, scale=TRUE), 
    M=M, 
    r=ncols_C
  )
  res[[paste0('_M', suffix)]] <- .res_resid_nam$M
  res[[paste0('_r', suffix)]] <- .res_resid_nam$r
  
  if (verbose) message('Decompose NAM')
  npcs <- 20 #max(10, round(max_frac_pcs * nrow(data$samplem)))
  
  npcs <- min(npcs, nrow(data$samplem) - 1) ## make sure you don't compute all SVs
  if (missing(npcs) | npcs > .5 * min(dim(NAM_))) {
    svd_res <- svd(scale(NAM_, center=FALSE, scale=TRUE))
  } else {
    svd_res <- RSpectra::svds(scale(NAM_, center=FALSE, scale=TRUE), k = npcs)
  }
  
  # A = U D V^T
  dim(svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v))
  ## d: singular value
  ## u, v: orthogonal matrix
  
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
  
  # d.uns['NAM.T'] contains the transpose of the NAM (neighborhoods by samples)
  # d.uns['NAM_sampleXpc'] contains the sample loadings of the principal components of the NAM
  # d.uns['NAM_nbhdXpc'] contains the neighborhood loadings of the principal components of the NAM
  # d.uns['NAM_svs'] contains the squared singular values of the NAM
  
  nam_res = res
  
  # d.uns['NAM.T'] contains the transpose of the NAM (neighborhoods by samples)
  # d.uns['NAM_sampleXpc'] contains the sample loadings of the principal components of the NAM
  # d.uns['NAM_nbhdXpc'] contains the neighborhood loadings of the principal components of the NAM
  # d.uns['NAM_svs'] contains the squared singular values of the NAM
  
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
  
  # prep data
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
    ps <- purrr::map_dbl(.tmp, 'p') # 
    r2s <- purrr::map_dbl(.tmp, 'r2')
    k_ <- which.min(ps) 
    return(list(k=ks[k_], p=ps[k_], r2=r2s[k_]))
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
  
  
  ncorrs <- V[, 1:k] %*% (sqrt(sv[1:k]) * beta/n)   # svs are actually eigenvalues, not SVs. I squared them to be consistent with python code.  # .res_resid_nam scales NAM columns, so this is correlation, not covariance.  
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
    y_ <- y_[, 1:Nnull]
    ycond_ <- scale(M %*% y_, center = FALSE, scale = TRUE)
    gamma_ <- crossprod(U[, 1:k], ycond_)
    nullncorrs <- abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_ / n)))
    
    maxcorr <- max(abs(ncorrs))
    fdr_thresholds <- seq(maxcorr/4, maxcorr, maxcorr/400)
    fdr_vals <- empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
    fdrs = data.frame(
      #         threshold = fdr_thresholds
      threshold = head(fdr_thresholds, -1),
      fdr = fdr_vals, 
      num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t)) 
    )
    # find maximal FDR<5% and FDR<10% sets
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
    misc = res ## Association results 
  )
  
  seurat_object@meta.data$cna_ncorrs <- ncorrs[colnames(seurat_object), , drop=TRUE]
  ## NOTE: If threshold was NULL, then no cells passed the significance threshold 
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


