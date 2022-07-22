# /caitlinch/treelikeness_metrics/func_metrics.R
# Caitlin Cherryh 2022

# This file contains functions for data manipulation and data analysis

collate.treelikeness.results <- function(alignment_path, experiment_number){
  # Function to take one alignment, get the parameters csv and treelikeness results csv, glue them together and return the new csv as a data.frame
  
  # Get paths to csv files
  alignment_params_file <- gsub("_output_alignment.fa", "_parameters.csv", alignment_path)
  alignment_treelikeness_file <- gsub("_output_alignment.fa", "_treelikeness_results.csv", alignment_path)
  
  # Check if files exist
  if (file.exists(alignment_params_file) == TRUE & file.exists(alignment_treelikeness_file) == TRUE){
    # Open the two csv files
    p_df <- read.csv(alignment_params_file)
    tl_df <- read.csv(alignment_treelikeness_file)
    
    # Check which experiment is being processed
    if (experiment_number == 1){
      # Remove unneeded columns from the p_df
      p_df <- p_df[,c("row_id", "uid" , "num_taxa", "num_trees", "tree_depth",  "num_reps", "alisim_gene_models", "alisim_gene_tree_length", "total_alignment_length", "sequence_type")]
    } else if (experiment_number == 2){
      # Remove unneeded columns from the p_df, making sure to retain the columns with recombination information
      p_df <- p_df[,c("row_id", "uid" , "num_taxa", "num_trees", "tree_depth", "recombination_value", "recombination_type", "num_reps", "alisim_gene_models", "alisim_gene_tree_length", 
                      "total_alignment_length", "sequence_type")]
    }
    
    # Check the uids match
    uids_identical <- grepl(p_df$uid, tl_df$input_alignment_path)
    
    # If uids are identical, combine the two dataframes
    if (uids_identical == TRUE){
      collated_df <- cbind(p_df, tl_df)
    } else {
      collated_df <- NULL
    }
    
  } else {
    # If the two files do not exist, return NULL
    collated_df <- NULL
  }
  
  # Return the collated df
  return(collated_df)
}

reformat.network.treelikeness.test.results.exp1 <- function(id, params_df, results_df){
  # For processing experiment 1 results
  # Function to collect the proportion of treelike alignments for each set of parameters
  # Reformat network treelikeness test results for pretty plotting
  
  # Extract relevant row
  row <- params_df[id, ]
  
  # Filter results_df to only rows including the relevent params
  row_results_df <- results_df[(results_df$num_taxa == row$num_taxa & 
                                  results_df$num_trees == row$num_trees & 
                                  results_df$tree_depth == row$tree_depth), ]
  # Find number of treelike and non-treelike results
  tl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Treelike"))
  ntl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Non-treelike"))
  # Find total number of alignments with these parameters
  n_als <- nrow(row_results_df)
  # Find proportion of treelike results
  prop_tl_als <- tl_results/n_als
  # Return the proportion of treelike alignments for this set of parameter values
  return(prop_tl_als)
}

reformat.network.treelikeness.test.results.exp2 <- function(id, params_df, results_df){
  # For processing experiment 2 results
  # Function to collect the proportion of treelike alignments for each set of parameters
  # Reformat network treelikeness test results for pretty plotting
  
  # Extract relevant row
  row <- params_df[id, ]
  
  # Filter results_df to only rows including the relevent params
  row_results_df <- results_df[(results_df$num_taxa == row$num_taxa & 
                                  results_df$tree_depth == row$tree_depth &
                                  results_df$recombination_value == row$recombination_value &
                                  results_df$recombination_type == row$recombination_type), ]
  # Find number of treelike and non-treelike results
  tl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Treelike"))
  ntl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Non-treelike"))
  # Find total number of alignments with these parameters
  n_als <- nrow(row_results_df)
  # Find proportion of treelike results
  prop_tl_als <- tl_results/n_als
  # Return the proportion of treelike alignments for this set of parameter values
  return(prop_tl_als)
}

