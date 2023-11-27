# caitlinch/treelikeness_metrics/code/func_empirical.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).

#### Required packages ####
library(ape)
library(seqinr)



#### Wrapper functions ####
bootstrap.wrapper <- function(bs_rep_al_paths, output_directory, 
                              splitstree_path, iqtree2_path, 
                              iqtree_model = "MFP", num_iqtree2_threads = "AUTO",
                              sequence_format = "AA", redo = FALSE,
                              number_parallel_cores = 1){
  ### Function to apply treelikeness metrics to a single replicate of a parametric bootstrap and return the output for each metric
  
  ## Copy each alignment into a separate folder in the output_directory
  # Extract the id from each filepath
  unique_ids <- unlist(lapply(1:length(bs_rep_al_paths), function(i){strsplit(basename(bs_rep_al_paths), "\\.")[[i]][1]}))
  # Create a new folder for each alignment folder and make each folder
  new_folders <- paste0(output_directory, unique_ids, "/")
  lapply(new_folders, dir.create)
  # Copy each alignment into its own folder
  file_format <- unlist(lapply(1:length(bs_rep_al_paths), function(i){tail(strsplit(basename(bs_rep_al_paths), "\\.")[[i]],1)}))
  alignment_paths <- paste0(output_directory, unique_ids, "/", unique_ids, ".", file_format)
  # Copy each alignment into its own new folder
  lapply(1:length(unique_ids), function(i){file.copy(from = bs_rep_al_paths[i], to = alignment_paths[i])})
  
  ## Run treelikeness tests
  lapply(alignment_paths,  treelikeness.metrics.empirical,
         splitstree_path = splitstree_path, 
         iqtree2_path = iqtree2_path, 
         iqtree_model = iqtree_model,
         num_iqtree2_threads = num_iqtree2_threads, 
         sequence_format = sequence_format, 
         redo = redo,
         mc.cores = number_parallel_cores)
  
  ## Create a nice dataframe of all the output values
  # Identify csv files
  all_files <- list.files(replicate_folder, recursive = TRUE)
  csv_files <- paste0(replicate_folder, grep("_treelikeness_results.csv", all_files, value = T))
  # Read in csv files
  csv_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)
  csv_df <- as.data.frame(do.call(rbind, csv_list))
  # Save csv list
  name_split <- strsplit(basename(alignment_path), "\\.")[[1]]
  op_id <- paste(name_split[1:(length(name_split) - 1)], collapse = ".")
  csv_df_file <- paste0(replicate_folder, "collated_results_", op_id, ".csv")
  write.csv(csv_df, file = csv_df_file, row.names = FALSE)
  
  ## Return output csv
  return(csv_df)
}


treelikeness.metrics.with.parametric.bootstrap <- function(i, df, tl_output_directory, 
                                                           splitstree_path, iqtree2_path, 
                                                           num_iqtree2_threads = "AUTO", sequence_format = "AA", 
                                                           redo = FALSE, number_parallel_cores = 1){
  ### Function to apply treelikeness metrics with parametric bootstrap and return the output + statistical significance for each metric
  
  ## Extract row of interest
  i_row <- df[i, ]
  i_id <- paste0(i_row$ID, ".alignment")
  i_directory <- paste0(tl_output_directory, i_row$ID, "/")
  if (dir.exists(i_directory) == FALSE){dir.create(i_directory)}
  tl_df_file <- paste0(i_directory, "collated_results_", i_id, ".csv")
  p_value_path <- paste0(i_directory, i_id, ".p_values.csv")
  
  if (file.exists(csv_df_file) == FALSE){
    print("copy alignment")
    ## Copy the alignment into the tl_output_directory (treelikeness output directory)
    i_alignment_path <- paste0(i_directory, basename(i_row$alignment_path))
    file.copy(from = i_row$alignment_path, to = i_alignment_path, overwrite = TRUE)
    
    ## Generate parametric bootstrap alignments (using Alisim in IQ-Tree2)
    #  Simulate an alignment of the same length as the original alignment, using the tree and model parameters 
    #         estimated from the original alignment, and copy the same gap positions from the original alignment
    alisim_command <- paste0(iqtree2_path, " -s ", i_alignment_path, " --alisim ", i_directory, "param_bs --num-alignments 100 --out-format fasta")
    system(alisim_command)
    # Extract best model of sequence evolution from the alisim output
    i_best_model <- extract.best.model(iqtree_file = paste0(i_alignment_path, ".iqtree"))
    # Collect all alignments
    i_files <- list.files(i_directory)
    i_all_alignments <- c(i_alignment_path,
                          paste0(i_directory, grep("param_bs", i_files, value = T)))
    
    ## Run treelikeness tests for all replicates
    tl_op <- mclapply(i_all_alignments,  treelikeness.metrics.empirical,
                      splitstree_path = splitstree_path, 
                      iqtree2_path = iqtree2_path, 
                      iqtree_model = i_best_model,
                      num_iqtree2_threads = num_iqtree2_threads, 
                      sequence_format = sequence_format, 
                      redo = redo,
                      mc.cores = number_parallel_cores)
    
    ## Create a nice dataframe of all the output values
    tl_df <- as.data.frame(do.call(rbind, tl_op))
    tl_df$dataset <- i_row$dataset
    tl_df$gene <- i_row$gene
    # Save csv list
    write.csv(tl_df, file = tl_df_file, row.names = FALSE)
  } else {
    tl_df <- read.csv(tl_df_file)
  }
  
  if (file.exists(p_value_path) == FALSE){
    print("p-value")
    ## Calculate p-values
    treelikeness_p_value <- calculate.p_value(value_vector = tl_df$tree_proportion, alignment_value = tl_df$tree_proportion[grep("gene", tl_df$unique_id)])
    scdf_mean_p_value <- calculate.p_value(value_vector = tl_df$sCF_mean, alignment_value = tl_df$sCF_mean[grep("gene", tl_df$unique_id)])
    scdf_med_p_value <- calculate.p_value(value_vector = tl_df$sCF_median, alignment_value = tl_df$sCF_median[grep("gene", tl_df$unique_id)])
    delta_plot_p_value <- calculate.p_value(value_vector = tl_df$mean_delta_plot_value, alignment_value = tl_df$mean_delta_plot_value[grep("gene", tl_df$unique_id)])
    
    ## Create p-value dataframe
    p_value_df <- data.frame("dataset" = i_row$dataset, "gene" = i_row$gene,
                             "tree_proportion" = treelikeness_p_value[["alignment_value"]], "tree_proportion_p_value" = treelikeness_p_value[["p_value_ecdf"]], 
                             "sCF_mean" = scdf_mean_p_value[["alignment_value"]], "sCF_mean_p_value" = scdf_mean_p_value[["p_value_ecdf"]], 
                             "sCF_median" = scdf_med_p_value[["alignment_value"]], "sCF_median_p_value"  = scdf_med_p_value[["p_value_ecdf"]],
                             "delta_plot_mean" = delta_plot_p_value[["alignment_value"]], "delta_plot_mean_p_value" =  delta_plot_p_value[["p_value_ecdf"]])
    write.csv(p_value_df, file = p_value_path, row.names = F)
  } else {
    p_value_df <- read.csv(p_value_path)
  }
  
  ## Return output csv
  return(p_value_df)
}


treelikeness.metrics.without.bootstrap <- function(i, df, tl_output_directory, 
                                                   splitstree_path, iqtree2_path, 
                                                   num_iqtree2_threads = "AUTO", sequence_format = "AA", 
                                                   redo = FALSE){
  ### Function to apply treelikeness metrics with parametric bootstrap and return the output + statistical significance for each metric
  
  ## Extract row of interest
  i_row <- df[i, ]
  i_id <- paste0(i_row$ID, ".alignment")
  i_directory <- paste0(tl_output_directory, i_row$ID, "/")
  i_alignment_path <- i_row$alignment_path
  if (dir.exists(i_directory) == FALSE){dir.create(i_directory)}
  tl_df_file <- paste0(i_directory, "collated_results_", i_id, ".csv")
  
  ## Estimate tree in IQ0Tree
  alisim_command <- paste0(iqtree2_path, " -s ", i_alignment_path, " -m MFP -nt ", num_iqtree2_threads)
  system(alisim_command)
  # Extract best model of sequence evolution from the alisim output
  i_best_model <- extract.best.model(iqtree_file = paste0(i_alignment_path, ".iqtree"))
  
  ## Run treelikeness tests for all replicates
  tl_op <- lapply(i_alignment_path,  treelikeness.metrics.empirical,
                  splitstree_path = splitstree_path, 
                  iqtree2_path = iqtree2_path, 
                  iqtree_model = i_best_model,
                  num_iqtree2_threads = num_iqtree2_threads, 
                  sequence_format = sequence_format, 
                  redo = redo)
  
  ## Create a nice dataframe of all the output values
  tl_df <- as.data.frame(do.call(rbind, tl_op))
  tl_df$dataset <- i_row$dataset
  tl_df$gene <- i_row$gene
  # Save csv list
  write.csv(tl_df, file = tl_df_file, row.names = FALSE)
  
  ## Return output csv
  return(tl_df)
}


extract.best.model <- function(iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .iqtree file
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    ## Check for a ModelFinder section:
    # Determine whether there is a ModelFinder section
    mf_ind <- grep("ModelFinder", iq_lines)
    # Determine whether there is a line detailing the best model
    bm_ind <- grep("Best-fit model according to", iq_lines)
    ## Check for a Substitution Process section:
    # Determine the starting line of this section
    sp_ind <- grep("SUBSTITUTION PROCESS", iq_lines)
    # Determine the line detailing the model used
    mos_ind <- grep("Model of substitution", iq_lines)
    ## Extract the best fit model from the .iqtree file:
    if ((identical(mf_ind, integer(0)) == FALSE) & (identical(bm_ind, integer(0)) == FALSE)){
      # If ModelFinder was run, extract the best model from the ModelFinder section of the .iqtree file
      # Extract the line containing the best fit model
      m_line <- iq_lines[bm_ind]
    } else if ((identical(sp_ind, integer(0)) == FALSE) & (identical(mos_ind, integer(0)) == FALSE)) {
      # If there is no ModelFinder section, extract the model used from the substitution process section
      m_line <- iq_lines[mos_ind]
    } else {
      m_line <- "NA:NA"
    }
    ## Format the model nicely for output: 
    # Split the line at the colon into two parts
    m_line_split <- strsplit(m_line, ":")[[1]]
    # If the best model is a single model, the length of m_line_split will be 2
    #     One section for the explanatory text and one for the model
    # If the best model is a partition model, it will have more than two sections when split by colons
    # Extract the second part of the line onwards (contains the best fit model)
    best_model <- m_line_split[2:length(m_line_split)]
    # If best_model is longer than 1, paste it together again using colons
    if (length(best_model) >1){
      best_model <- paste(best_model, collapse = ":")
    }
    # Remove any white space from the best model
    best_model <- gsub(" ", "", best_model)
  } else if (file.exists(iqtree_file) == FALSE){
    # If the iqtree_file doesn't exist, return NA
    best_model = NA
  } # end if (file.exists(iqtree_file) == TRUE){
  # Return the best model from the iqtree_file (if the file exists)
  return(best_model)
}


#### Apply all treelikeness metrics ####
treelikeness.metrics.empirical <- function(alignment_path, splitstree_path, iqtree2_path, 
                                           iqtree_model = "MFP", num_iqtree2_threads = "AUTO",
                                           sequence_format = "AA", redo = FALSE, best.tests.only = TRUE){
  #### Function to take one simulated alignment, apply all treelikeness metrics and return results in a dataframe
  
  ### Print alignment path
  print(alignment_path)
  
  ### Prepare variables and output file names for run
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Identify the prefix for the alignment path
  file_prefix <- tail(strsplit(basename(alignment_path), "\\.")[[1]], 1)
  # Get unique id for the alignment
  unique_id <- gsub(paste0(".", file_prefix), "", basename(alignment_path))
  # Create name for output dataframes
  df_name <- paste0(replicate_folder, unique_id, "_treelikeness_results.csv")
  
  ### Convert alignment to NEXUS
  # Open alignment file
  f <- read.FASTA(alignment_path, type = "AA")
  # Write out alignment file as nexus format
  nexus_alignment_path <- paste0(replicate_folder, unique_id, ".nex")
  write.nexus.data(f, file = nexus_alignment_path, format = "protein", datablock = FALSE, interleaved = FALSE)
  
  ### Extract number of taxa in file
  number_of_taxa <- length(f)
  
  ### Prepare results dataframe
  # Check whether dataframe .csv file already exists. If it does, import the dataframe. If it doesn't, make it by running all treelikeness metrics
  if (file.exists(df_name) == TRUE & redo == FALSE){
    ## Read in the results csv file
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE | redo == TRUE){
    ### Apply all treelikeness test statistics to generate the results csv file
    
    ## Apply Site concordance factors with likelihood (Minh et. al. 2020): --scfl (iqtree2 v2.2.2)
    # Specify model for calculating site concordance factors
    if ( ( is.na(iqtree_model) == TRUE ) | (iqtree_model == "MFP") ){
      ## If either the input model is NA or is MFP, extract the best model as identified by the likelihood mapping IQ-Tree run
      ##      and apply that model to calculate the site concordance factors
      iqtree_file <- paste0(alignment_path, ".iqtree")
      # Extract the best model from the IQ-Tree run for the likelihood mapping
      scfl_model <- identify.best.model(iqtree_file)
    } else {
      ## If the model is provided as input, use that model to calculate the site concordance factors
      scfl_model <- iqtree_model
    }
    # Calculate the site concordance factors
    scfl_output <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = 100, 
                        substitution_model = scfl_model, include.prefix = TRUE, prefix = unique_id, force.redo = TRUE)
    
    ## Apply Delta plots (Holland et. al. 2002)
    mldist_file <- paste0(alignment_path, ".mldist")
    mean_delta_plot_value <- delta.plot.empirical( dist_matrix = mldist.matrix(mldist_file) )
    
    ## Apply tree proportion (new test)
    tree_proportion <- tree.proportion.long(nexus_alignment_path, sequence_format = sequence_format, model = NA, 
                                            remove_trivial_splits = TRUE, check_iqtree_log_for_identical_sequences = FALSE, 
                                            run_splitstree = TRUE, splitstree_path = splitstree_path)
    
    ## If desired, calculate the tests that didn't work well in simulation study
    if (best.tests.only == FALSE){
      ## Apply Likelihood mapping (Strimmer and von Haeseler 1997)
      # Check iqtree_model - if NA, run ModelFinder with "MFP"
      if ( is.na(iqtree_model) == TRUE ){ 
        lm_model <- "MFP"
      } else {
        lm_model <- iqtree_model
      }
      # Apply likelihood mapping method
      lm <- likelihood.mapping.empirical(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = lm_model, 
                                         number_of_taxa = number_of_taxa, sequence_format = sequence_format)
      
      ## Apply Network Treelikeness Test (Huson and Bryant 2006)
      ntlt <- network.treelikeness.test(nexus_alignment_path, splitstree_path, sequence_format = sequence_format, nexus.file.format = TRUE)
      
      ## Apply Cunningham test (Cunningham 1975)
      # Find the .mldist and .treefile output paths from IQ-Tree
      mldist_file <- paste0(alignment_path, ".mldist")
      tree_file <- paste0(alignment_path, ".treefile")
      # Calculate Cunningham metric
      cunningham_metric <- cunningham.test.empirical(mldist_file, tree_file)
      
      ## Assemble results into a dataframe and save
      results_vec <- c(unique_id, lm, scfl_output$mean_scf, scfl_output$median_scf, min(scfl_output$all_scfs), max(scfl_output$all_scfs), 
                       ntlt, mean_delta_plot_value, cunningham_metric, tree_proportion, alignment_path)
      results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
      names_vec <- c("unique_id", "LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets",
                     "LM_total_num_quartets", "LM_proportion_resolved_quartets", "sCF_mean", "sCF_median", "sCF_min", "sCF_max",
                     "NetworkTreelikenessTest", "mean_delta_plot_value", "Cunningham_test", "tree_proportion", "input_alignment_path")
      names(results_df) <- names_vec
    } else {
      ## Assemble results into a dataframe and save
      results_vec <- c(unique_id, scfl_output$mean_scf, scfl_output$median_scf, min(scfl_output$all_scfs), max(scfl_output$all_scfs),
                       mean(scfl_output$all_sdf1), median(scfl_output$all_sdf1), min(scfl_output$all_sdf1), max(scfl_output$all_sdf1),
                       mean(scfl_output$all_sdf2), median(scfl_output$all_sdf2), min(scfl_output$all_sdf2), max(scfl_output$all_sdf2),
                       mean_delta_plot_value, tree_proportion, alignment_path)
      results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
      names_vec <- c("unique_id", "sCF_mean", "sCF_median", "sCF_min", "sCF_max", 
                     "sDF1_mean", "sDF1_median", "sDF1_min", "sDF1_max",
                     "sDF2_mean", "sDF2_median", "sDF2_min", "sDF2_max",
                     "mean_delta_plot_value", "tree_proportion", "input_alignment_path")
      names(results_df) <- names_vec
    }
    
    ## Save the csv 
    write.csv(results_df, file = df_name, row.names = FALSE)
  } # end run treelikeness tests
  # Return the dataframe
  return(results_df)
} # end function



#### Treelikeness test functions ####
likelihood.mapping.empirical <- function(alignment_path, iqtree2_path, iqtree2_number_threads = 1, substitution_model = "MFP", 
                                         number_of_taxa = NA, sequence_format = "DNA"){
  # Function to call IQ-Tree and create a likelihood map for the alignment
  
  ## Create the likelihood map
  # Check whether likelihood mapping or IQ-Tree have run before. 
  # If one or both haven't run IQ-Tree to create the likelihood map
  iq_file <- paste0(alignment_path, ".iqtree")
  map_file <- paste0(alignment_path, ".lmap.eps")
  if ((file.exists(iq_file) == FALSE) | (file.exists(map_file) == FALSE)){
    number_of_quartets <- 25 * as.numeric(number_of_taxa)
    call <- paste0(iqtree2_path, " -s ", alignment_path, " -m ", substitution_model, " -nt ", iqtree2_number_threads, " -lmap ", number_of_quartets, " -redo -safe")
    system(call)
  }
  
  if (file.exists(iq_file) == TRUE){
    # Extract results from likelihood map
    iq_log <- readLines(iq_file)
    ind <- grep("Number of fully resolved  quartets",iq_log)
    resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    ind <- grep("Number of partly resolved quartets",iq_log)
    partly_resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    ind <- grep("Number of unresolved",iq_log)
    unresolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
    total_q <- (resolved_q+partly_resolved_q+unresolved_q)
    prop_resolved <- resolved_q/total_q
    # Collate results into a vector
    lm_results <- c(resolved_q, partly_resolved_q, unresolved_q, total_q, prop_resolved)
  } else if (file.exists(iq_file) == FALSE){
    # Create a vector noting that the .iqtree file does not exist
    lm_results <- rep("no_iqtree_file", 5)
  }
  
  ## Rename vector of results
  names(lm_results) <- c("num_resolved_quartets", "num_partly_resolved_quartets", "num_unresolved_quartets",
                         "total_num_quartets", "proportion_resolved_quartets")
  
  ## Return results
  return(lm_results)
}



## Delta plots (Holland et. al. 2002)
delta.plot.empirical <- function(dist_matrix){
  # This function takes an alignment, calculates a distance matrix for the alignment, and the applies the
  # `ape` function `delta.plot`. We take the mean delta plot value as the test statistic. 
  
  ## Take the pairwise distance matrix from the input distance matrix
  pdm <- dist_matrix
  ## Call ape::delta.plot function
  # Set the number of intervals for the delta plot
  dp_intervals = 100
  # Make a delta.plot based on the pairwise distance matrix
  dp <- delta.plot(pdm, k = dp_intervals, plot = FALSE)
  ## To calculate the mean delta q from ALL quartets:
  # Create two vectors, one containing the counts and one containing the midpoint of each interval
  # To determine the midpoint of each interval, first find the intervals (e.g. for k = 2, there will be 2 intervals: 0-0.5 and 0.5-1),
  #     and the midpoint of each interval will be 0.25 and 0.75 (the mean of the start and endpoint of each interval)
  interval_midpoint = (seq(0,0.999,1/(dp_intervals)) + (0.5 * (0 + seq(0,1,1/(dp_intervals))[2])))
  interval_count = dp$counts
  # To calculate the mean delta_q, calculate a weighted mean from the 
  mean_dq <- weighted.mean(interval_midpoint, interval_count)
  ## To calculate the mean delta bar (the mean value across all taxa e.g. the mean of the mean values for each taxa):
  # Calculate the mean delta bar (delta bar = the mean delta value for each observation/taxa)
  mean_db <- mean(dp$delta.bar)
  ## Return values to outside function
  # Return the mean delta bar (the mean delta q value across all taxa)
  return(mean_db)
}



cunningham.test.empirical <- function(mldist_file, tree_file){
  ## Function to estimate what proportion of the variance in the data is represented by the tree
  
  ## Test steps:
  # 1. Calculate the observed distances from the alignment (d_ij)
  # 2. Calculate the predicted distances from the tree (p_ij)
  # 3. Calculate the total sum of squares. TSS = sum of (d_ij)^2
  # 4. Calculate the residual sum of squares. RSS = sum of (p_ij - d_ij)^2
  # 5. Calculate the R^2 = (TSS - RSS)/ RSS
  
  ## Test:
  # 1. Calculate the observed distances (d_ij)
  mldist_matrix_raw <- mldist.matrix(mldist_file) # Feed in distance matrix (mldist file from IQ-Tree run)
  al_mat <- as.dist(mldist_matrix_raw)
  d_ij <- as.vector(al_mat) # observed distances between taxa i and j
  
  # 2. Calculate the predicted distances (p_ij)
  # Open the tree
  t <- read.tree(tree_file)
  # Extract the distance matrix from the tree
  t_cophenetic_mat <- cophenetic.phylo(t) # in substitutions per site
  # Now reorder the t_mat so the taxa are in the same order
  mldist_taxa <- attr(al_mat, "Labels")
  t_ordering_mat <- as.matrix(t_cophenetic_mat)[mldist_taxa, mldist_taxa]
  t_mat <- as.dist(t_ordering_mat)
  p_ij <- as.vector(t_mat) # predicted distances between taxa i and j
  
  # 3. Calculate the TSS = sum of (d_ij)^2
  TSS = sum((d_ij)^2)
  
  # 4. Calculate the RSS = sum of (p_ij - d_ij)^2
  RSS = sum((p_ij - d_ij)^2)
  
  # 5. Calculate the R^2
  r_squared = 1 - ((RSS)/(TSS))
  
  ## Return the r_squared value
  return(r_squared)
}



tiger.empirical <- function(alignment_path, fast_TIGER_path,
                            sequence_format = "DNA"){
  ## Function to take one empirical alignment, apply fast TIGER and return results in a dataframe
  
  # Print alignment path
  print(alignment_path)
  
  ## Prepare variables and output file names for run
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Get unique id for the alignment
  unique_id <- paste(gsub("_output_alignment", "", unlist(strsplit(basename(alignment_path), "\\."))[1:(length(unlist(strsplit(basename(alignment_path), "\\."))) - 1)]), collapse = ".") 
  
  # Create name for output dataframe
  df_name <- paste0(replicate_folder, unique_id, "_tiger_results.csv")
  
  if (file.exists(df_name) == TRUE){
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE){
    # Apply TIGER (Cummins and McInerney 2011)
    mean_tiger_value <- TIGER(alignment_path, fast_TIGER_path, sequence_format = sequence_format)
    
    # Assemble results into a dataframe and save
    results_vec <- c(unique_id, mean_tiger_value)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("uid", "mean_TIGER_value")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE) 
  }
  
  # Return the tiger dataframe
  return(results_df)
} # end function



#### Utility functions
mldist.matrix <- function(mldist_file){
  ## Open and format the .mldist file from IQ-Tree as a matrix
  # Open .mldist file as a text file
  mldist_text <- readLines(mldist_file)
  # Remove the first line (only gives number of taxa in matrix)
  mldist_text <- mldist_text[2:length(mldist_text)]
  # Transform text file into the ML dist matrix
  mat_list <- lapply(mldist_text, process.distance.matrix.text.line)
  dist_mat <- as.matrix(do.call(rbind, mat_list))
  # Add species names as row and column names
  mat_species <- unlist(lapply(mldist_text, extract.distance.matrix.species))
  rownames(dist_mat) <- mat_species
  colnames(dist_mat) <- mat_species
  # Return the distance matrix
  return(dist_mat)
}

process.distance.matrix.text.line <- function(line){
  ## Function to nicely format each line of the distml file output by IQ-Tree
  # Split the line at the spaces
  split_line <- strsplit(line, split = " ")[[1]]
  # Remove the first element from the line
  split_line <- split_line[2:length(split_line)]
  # Remove any empty elements
  split_line <- split_line[which( ! (split_line == "") )]
  # Change line into numbers
  split_line <- as.numeric(split_line)
  # Return the split_line as a vector
  return(split_line)
}

extract.distance.matrix.species <- function(line){
  ## Function to nicely format each line of the distml file output by IQ-Tree
  # Split the line at the spaces
  split_line <- strsplit(line, split = " ")[[1]]
  # Extract the species from the line text
  line_species <- split_line[1]
  # Return the species name
  return(line_species)
}



#### Generate alignment ####
generate.mimic.alignment <- function(output_prefix, alignment_path, iqtree_path, 
                                     substitution_model, output_format = "fasta", sequence_type = "AA",
                                     num_output_alignments = 1){
  ## This function uses the topology-unlinked partition model in Alisim to generate a sequence alignment
  #     containing multiple concatenated genes, each with its own tree topology and branch lengths
  
  # Assemble function call 
  function_call <- paste0(iqtree_path, " --alisim ", output_prefix, " -s ", alignment_path, 
                          " -m ", substitution_model, " --seqtype ", sequence_type,
                          " --out-format ", output_format, " --num-alignments ", num_output_alignments)
  # Invoke the OS command and call IQ-Tree
  system(function_call)
  
  # Print completion statement
  print("Alisim (IQ-Tree2) run complete")
}



#### Statistics ####
calculate.p_value <- function(value_vector, alignment_value){
  ## Function that given two vectors (one of test statistic values, and one of ids), calculates the p-value for that alignment
  ##    We are interested in whether the alignment value is LOWER than the bootstrap replicate values - which would indicate
  ##    that the parametric bootstrap replicates are MORE treelike then the alignment
  
  # Check whether the alignment value is present or is NA
  if (is.na(alignment_value) == TRUE){
    # If the alignment value is NA, can't calculate a score
    # If it wasn't possible to calculate a score for the alignment, output an NA
    p_value_cdf <- NA
  } else {
    # If it's possible to calculate a p-value, calculate one using the CDF
    # Compute an empirical cumulative distribution function (ecdf)
    cdf <- ecdf(value_vector)
    # Calculate the p-value by feeding the alignment test statistic value into the ecdf
    #     The probability that a randomly selected test statistic value has a value less than or equal 
    #     to the alignment value is the CDF at the alignment value
    p_value_cdf <- cdf(alignment_value)
  }
  
  # return the p-value
  op_vector <- c(alignment_value, p_value_cdf)
  names(op_vector) <- c("alignment_value", "p_value_ecdf")
  return(op_vector)
}


calculate.all.p_values <- function(output_df, test_statistic_names){
  ## Function to apply p-value to multiple test statistics at once, given a vector of column names
  
  # Apply the test to each column one at a time
  output_vector <- c()
  for (x in test_statistic_names){
    temp_p_value <- calculate.p_value(value_vector = output_df[[x]], alignment_value = output_df[[x]][which(output_df$unique_id == "alignment")])
    output_vector <- c(output_vector, temp_p_value)
  }
  # Create a dataframe
  p_value_df <- data.frame(test_statistic = test_statistic_names,
                           test_statistic_value = output_vector[c(T,F)],
                           p_value = output_vector[c(F,T)])
  # Return the dataframe of p-values
  return(p_value_df)
}


### Partition and file management
split.partitions <- function(alignment_file, partition_file, gene_output_directory){
  ## Open the dataset, split into individual genes and save each gene
  
  # Open partition file
  lines <- readLines(partition_file)
  # Extract all lines with a charset
  charset_lines <- grep("charset", lines, ignore.case = TRUE, value = TRUE)
  # Split the charset lines at the "="
  gene_name_chunks <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[1]}))
  gene_lines <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[2]}))
  # Split the genes into chunks by breaking at the commas ","
  gene_chunks <- unlist(strsplit(gene_lines, ","))
  # Format the gene chunks nicely
  gene_chunks_nospace <- gsub(" ", "", gene_chunks)
  gene_chunks_noend <- gsub(";", "", gene_chunks_nospace)
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[2]})))
  # Format the names
  gene_names <- unlist(lapply(strsplit(gene_name_chunks, " "), function(x){x[[2]]}))
  # Make a table of gene name, start, and end position
  gene_table <- data.frame(name = gene_names, start = gene_start, end = gene_end)
  # Open alignment
  al <- read.FASTA(alignment_file, type = "AA")
  al <- as.matrix(al)
  alignment_taxa <- attr(al, "dimnames")[[1]]
  # Create object for saving file paths
  output_vector <- c()
  # Iterate through each row of the table and save the genes
  for (i in 1:nrow(gene_table)){
    # Extract row using the iterator
    gene_row <- gene_table[i,]
    # Extract that section of the alignment
    gene_al <- al[, gene_row$start:gene_row$end]
    # Remove any empty taxa
    taxa_to_remove <- c()
    chars_to_check <- c("-", "?")
    for (j in 1:length(alignment_taxa)){
      unique_chars <- unique(as.character(as.list(gene_al[j, ]))[[1]])
      # If the unique_chars is only either or both these characters: "-", "?"
      #     then don't add keep this tip in the output alignment
      #     No sites for tree estimation!
      if (length(unique_chars) == 1){
        if (unique_chars == "?" | unique_chars == "-"){
          taxa_to_remove <- c(taxa_to_remove, j)
        }
      } else if (length(unique_chars) == 2){
        if (setequal(unique_chars, c("?", "-")) == TRUE){
          taxa_to_remove <- c(taxa_to_remove, j)
        }
      }
    }
    # Identify taxa to keep
    taxa_to_keep <- setdiff(1:length(alignment_taxa), taxa_to_remove)
    # Remove the tips with only gaps/unknown characters
    gene_al_complete <- gene_al[taxa_to_keep,]
    # Save gene
    gene_al_complete <- as.list(gene_al_complete)
    gene_file_path <- paste0(gene_output_directory, gene_row$name, ".fa")
    write.fasta(sequences = gene_al_complete, names = names(gene_al_complete), file.out = gene_file_path)
    # Add path to output vector
    output_vector <- c(output_vector, gene_file_path)
  }
  # Output vector of gene file paths
  return(output_vector)
}



