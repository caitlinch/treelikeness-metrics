# caitlinch/treelikeness_metrics/code/func_empirical.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).


#### Treelikeness test functions ####
treelikeness.metrics.empirical <- function(alignment_path, 
                                           iqtree2_path, splitstree_path, 
                                           phylogemetric_path, fast_TIGER_path, 
                                           supply_number_of_taxa = FALSE, number_of_taxa = NA, 
                                           num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                                           iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
                                           num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE, 
                                           run_splitstree_for_tree_proportion = FALSE, sequence_format = "DNA", 
                                           apply.TIGER = FALSE, redo = FALSE){
  ## Function to take one empirical alignment, apply all treelikeness metrics and return results in a dataframe
  
  # Print alignment path
  print(alignment_path)
  
  ## Prepare variables and output file names for run
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Get unique id for the alignment
  unique_id <- paste(gsub("_output_alignment", "", unlist(strsplit(basename(alignment_path), "\\."))[1:(length(unlist(strsplit(basename(alignment_path), "\\."))) - 1)]), collapse = ".") 
  # Get list of files in the replicate_folder
  all_folder_files <- list.files(replicate_folder)
  aln_folder_files <- grep(unique_id, all_folder_files, value = TRUE)
  # Create name for output dataframes
  df_name <- paste0(replicate_folder, unique_id, "_treelikeness_results.csv")
  collated_df_name <- paste0(replicate_folder, unique_id, "_collated_alignment_results.csv")
  
  ## Prepare results dataframe
  # Check whether dataframe .csv file already exists. If it does, import the dataframe. If it doesn't, make it by running all treelikeness metrics
  if (file.exists(df_name) == TRUE & redo == FALSE){
    ## Read in the results csv file
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE | redo == TRUE){
    ## Apply all treelikeness test statistics to generate the results csv file
    # Determine the number of taxa (needed for number of quartets in likelihood mapping and sCFs)
    if (supply_number_of_taxa == TRUE & !is.na(number_of_taxa)){
      # If the number of taxa is supplied as an input variable, use the input value
      n_tree_tips = number_of_taxa
    } else {
      # If the number of taxa isn't supplied as an input variable, determine it by identifying the number of taxa in the fasta alignment
      n_tree_tips <- length(read.FASTA(alignment_path))
    }
    
    # Estimate the ML tree using IQ-Tree, extract details about the model of sequence evolution, and conduct a likelihood mapping analysis
    # Apply Likelihood mapping (Strimmer and von Haeseler 1997)
    lm <- likelihood.mapping(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = iqtree_substitution_model, 
                             number_of_taxa = n_tree_tips, sequence_format = sequence_format)
    
    # Extract parameters for estimating the distance matrix from the .iqtree file 
    if (sequence_format == "DNA"){
      # Possible DNA models in dist.ml function: "JC69", "F81"
      if (is.na(distance_matrix_substitution_method) == TRUE){
        # If no distance_matrix_substitution_method is provided, use the default
        distance_matrix_substitution_method = "JC69"
      } else {
        # If distance_matrix_substitution_method is provided, use it
        distance_matrix_substitution_method = distance_matrix_substitution_method
      }
      
      model_params <- get.model.parameters.from.iqtree.file(alignment_path, sequence_format = "DNA")
      state_frequencies  <- model_params$state_frequencies
      Q_vector <- model_params$Q_vector
      num_rate_categories <- model_params$num_rate_categories
      best_iqtree_model = model_params$best_fit_model
      
    } # end if sequence alignment == "DNA"
    
    # Apply Site concordance factors with likelihood (Minh et. al. 2020): --scfl (iqtree2 v2.2.2)
    scfl <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                 substitution_model = iqtree_substitution_model)
    # Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(alignment_path, splitstree_path, sequence_format = sequence_format)
    # Apply Delta plots (Holland et. al. 2002)
    mean_delta_plot_value <- mean.delta.plot.value(alignment_path, sequence_format = sequence_format, substitution_model = distance_matrix_substitution_method,
                                                   base_frequencies = state_frequencies, Q_matrix = Q_vector, number_of_rate_categories = num_rate_categories)
    # Apply Q-residuals (Gray et. al. 2010)
    q_residual_results <- q_residuals(alignment_path, phylogemetric_path, sequence_format = sequence_format, phylogemetric_number_of_threads = num_phylogemetric_threads)
    mean_q_residual <- q_residual_results$mean_q_residual
    if (apply.TIGER == TRUE){
      # Apply TIGER (Cummins and McInerney 2011)
      mean_tiger_value <- TIGER(alignment_path, fast_TIGER_path, sequence_format = sequence_format)
    } else if (apply.TIGER == FALSE){
      mean_tiger_value <- "no_TIGER_run"
    }
    # Apply Cunningham test (Cunningham 1975)
    cunningham_metric <- cunningham.test(alignment_path, sequence_format, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, iqtree_substitution_model = iqtree_substitution_model, 
                                         distance_matrix_substitution_model = distance_matrix_substitution_method,
                                         base_frequencies = state_frequencies, Q_matrix = Q_vector, number_of_rate_categories = num_rate_categories)
    # Apply tree proportion (new test)
    tree_proportion <- tree.proportion.long(alignment_path, sequence_format = sequence_format, model = distance_matrix_substitution_method, 
                                            remove_trivial_splits = tree_proportion_remove_trivial_splits, check_iqtree_log_for_identical_sequences = FALSE, 
                                            run_splitstree = run_splitstree_for_tree_proportion, splitstree_path = splitstree_path,
                                            base_frequencies = state_frequencies, Q_matrix = Q_vector, number_of_rate_categories = num_rate_categories)
    
    # Assemble results into a dataframe and save
    results_vec <- c(lm, scfl$mean_scf, scfl$median_scf, min(scfl$all_scfs), max(scfl$all_scfs), ntlt, mean_delta_plot_value, mean_q_residual, mean_tiger_value,
                     cunningham_metric, tree_proportion, alignment_path, iqtree_substitution_model, best_iqtree_model, 
                     paste(state_frequencies, collapse = ","), paste(Q_vector, collapse = ","), num_rate_categories, distance_matrix_substitution_method)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets", "LM_total_num_quartets", "LM_proportion_resolved_quartets",
                   "sCF_mean", "sCF_median", "sCF_min", "sCF_max", "NetworkTreelikenessTest", "mean_delta_plot_value", "mean_Q_residual", "mean_TIGER_value",
                   "Cunningham_test", "tree_proportion", "input_alignment_path", "IQ-Tree_input_model", "IQ-Tree_best_fit_BIC_model",
                   "IQ-Tree_base_frequencies", "IQ-Tree_Q_matrix_lower_diagonal", "IQ-Tree_number_of_rate_categories", "distance_matrix_substitution_model")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE)
    
  } # end applying treelikeness tests
  
} # end function


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