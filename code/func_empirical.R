# caitlinch/treelikeness_metrics/code/func_empirical.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).



#### Wrapper functions ####
empirical.treelikeness.test.wrapper <- function(alignment_path, 
                                                iqtree2_path, splitstree_path, 
                                                phylogemetric_path, fast_TIGER_path, 
                                                supply_number_of_taxa = FALSE, number_of_taxa = NA, 
                                                num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                                                iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
                                                num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE, 
                                                run_splitstree_for_tree_proportion = FALSE, sequence_format = "DNA", 
                                                apply.TIGER = FALSE, redo = FALSE){
  ### Function to apply treelikeness metrics to an empirical dataset, perform a parametric bootstrap with 100 replicates, 
  #       and return the output and p-values for each metric
  
  ## Copy alignment to new dataframe
  
  ## Apply treelikeness metric to alignment
  
  ## Extract parameters from the IQ-Tree run of the alignment
  
  ## Perform 100 parametric bootstrap replicates
  
  ## Calculate p-values
  
  ## Assemble output csv
}




#### Apply all treelikeness metrics ####
treelikeness.metrics.empirical <- function(alignment_path, 
                                           iqtree2_path, splitstree_path, 
                                           phylogemetric_path, fast_TIGER_path, 
                                           num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                                           iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
                                           num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE, 
                                           run_splitstree_for_tree_proportion = FALSE, sequence_format = "DNA", 
                                           apply.TIGER = FALSE, redo = FALSE){
  ## Function to take one simulated alignment, apply all treelikeness metrics and return results in a dataframe
  
  # Print alignment path
  print(alignment_path)
  
  ## Prepare variables and output file names for run
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Get unique id for the alignment
  unique_id <- basename(dirname(alignment_path))
  # Create name for output dataframes
  df_name <- paste0(replicate_folder, unique_id, "_treelikeness_results.csv")
  
  ## Convert alignment to NEXUS
  # Open alignment file
  f <- read.FASTA(alignment_path, type = "AA")
  # Write out alignment file as nexus format
  nexus_alignment_path <- paste0(replicate_folder, unique_id, ".nex")
  write.nexus.data(f, file = nexus_alignment_path, format = "protein", datablock = FALSE, interleaved = FALSE)
  
  # Extract number of taxa in file
  number_of_taxa <- length(f)
  
  
  ## Prepare results dataframe
  # Check whether dataframe .csv file already exists. If it does, import the dataframe. If it doesn't, make it by running all treelikeness metrics
  if (file.exists(df_name) == TRUE & redo == FALSE){
    ## Read in the results csv file
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE | redo == TRUE){
    ## Apply all treelikeness test statistics to generate the results csv file
    
    # Apply Likelihood mapping (Strimmer and von Haeseler 1997)
    lm <- likelihood.mapping.empirical(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = iqtree_substitution_model, 
                                       number_of_taxa = number_of_taxa, sequence_format = sequence_format)
    
    # Apply Site concordance factors with likelihood (Minh et. al. 2020): --scfl (iqtree2 v2.2.2)
    scfl <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                 substitution_model = iqtree_substitution_model)
    
    # Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(nexus_alignment_path, splitstree_path, sequence_format = sequence_format, nexus.file.format = TRUE)
    
    # Apply Delta plots (Holland et. al. 2002)
    # Open ML dist matrix from IQ-Tree run
    
    # Calculate delta plot values
    mean_delta_plot_value <- delta.plot.empirical(alignment_path, sequence_format = sequence_format, substitution_model = distance_matrix_substitution_method)
    
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
    cunningham_metric <- cunningham.test(alignment_path, sequence_format, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, 
                                         iqtree_substitution_model = iqtree_substitution_model, 
                                         distance_matrix_substitution_model = distance_matrix_substitution_method,
                                         base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA)
    
    # Apply tree proportion (new test)
    tree_proportion <- tree.proportion.long(alignment_path, sequence_format = sequence_format, model = distance_matrix_substitution_method, 
                                            remove_trivial_splits = tree_proportion_remove_trivial_splits, check_iqtree_log_for_identical_sequences = FALSE, 
                                            run_splitstree = run_splitstree_for_tree_proportion, splitstree_path = splitstree_path,
                                            base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA)
    
    # Assemble results into a dataframe and save
    results_vec <- c(lm, scfl$mean_scf, scfl$median_scf, min(scfl$all_scfs), max(scfl$all_scfs), ntlt, mean_delta_plot_value, mean_q_residual, mean_tiger_value,
                     cunningham_metric, tree_proportion, alignment_path)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets", "LM_total_num_quartets", "LM_proportion_resolved_quartets",
                   "sCF_mean", "sCF_median", "sCF_min", "sCF_max", "NetworkTreelikenessTest", "mean_delta_plot_value", "mean_Q_residual", "mean_TIGER_value",
                   "Cunningham_test", "tree_proportion", "input_alignment_path")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE)
  } # end run treelikeness tests
  
} # end function



#### Treelikeness test functions ####
## Delta plots (Holland et. al. 2002)
delta.plot.empirical <- function(alignment_path, sequence_format = "AA", substitution_model = "WAG"){
  # This function takes an alignment, calculates a distance matrix for the alignment, and the applies the
  # `ape` function `delta.plot`. We take the mean delta plot value as the test statistic. 
  
  ## Open the alignment as a phyDat object
  p <- phyDat(read.FASTA(alignment_path, type = sequence_format), type = sequence_format)
  ## Calculate a distance matrix of pairwise distances from AA sequences using a model of AA substitution
  # Default model is WAG
  pdm <- dist.ml(p, model = substitution_model)
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
  }
  else if (file.exists(iq_file) == FALSE){
    # Create a vector noting that the .iqtree file does not exist
    lm_results <- rep("no_iqtree_file", 5)
  }
  
  ## Rename vector of results
  names(lm_results) <- c("num_resolved_quartets", "num_partly_resolved_quartets", "num_unresolved_quartets",
                         "total_num_quartets", "proportion_resolved_quartets")
  
  ## Return results
  return(lm_results)
}



#### Model extraction functions ####



#### Parametric bootstrap replicate functions ####




