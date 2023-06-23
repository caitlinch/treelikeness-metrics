# caitlinch/treelikeness_metrics/code/func_metrics.R
# Caitlin Cherryh 2023

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
      p_df <- p_df[,c("row_id", "uid" , "num_taxa", "num_trees", "tree_depth_coalescent", "recombination_value", "recombination_type", "num_reps", "alisim_gene_models", "alisim_gene_tree_length", 
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


collate.empirical.treelikeness.results <- function(alignment_path){
  # Function to take one alignment, get the parameters csv and treelikeness results csv, glue them together and return the new csv as a data.frame
  
  # Get paths to csv files
  alignment_params_file <- gsub("_output_alignment.fa", "_parameters.csv", alignment_path)
  alignment_treelikeness_file <- gsub("_output_alignment.fa", "_treelikeness_results.csv", alignment_path)
  
  # Check if files exist
  if (file.exists(alignment_params_file) == TRUE & file.exists(alignment_treelikeness_file) == TRUE){
    # Open the two csv files
    p_df <- read.csv(alignment_params_file)
    tl_df <- read.csv(alignment_treelikeness_file)
    
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


reformat.network.treelikeness.test.results.empirical <- function(id, params_df, results_df){
  # For processing experiment 2 results
  # Function to collect the proportion of treelike alignments for each set of parameters
  # Reformat network treelikeness test results for pretty plotting
  
  # Extract relevant row
  row <- params_df[id, ]
  
  # Filter results_df to only rows including the relevent params
  row_results_df <- results_df[(results_df$num_taxa == row$num_taxa & 
                                  results_df$codon_position == row$codon_position &
                                  results_df$DNA_type == row$DNA_type), ]
  # Find number of treelike and non-treelike results
  tl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Treelike"))
  ntl_results <- length(which(row_results_df$NetworkTreelikenessTest == "Non-treelike"))
  other_results <- length(which(!row_results_df$NetworkTreelikenessTest %in% c("Treelike", "Non-treelike")))
  # Find total number of alignments with these parameters
  n_als <- nrow(row_results_df)
  # Find proportion of treelike results
  prop_tl_als <- tl_results/n_als
  # Return the proportion of treelike alignments for this set of parameter values
  return(prop_tl_als)
}


copy.empirical.alignment <- function(row_id, data_df){
  # Small function to take a row from a dataframe, save a copy of the alignment, and save the row containing the alignment parameters as a csv
  
  # Extract the row
  row <- data_df[row_id, ]
  
  # Check whether folder for this alignment exists and create it if necessary
  al_dir <- dirname(row$output_alignment_path)
  if (dir.exists(al_dir) == FALSE){dir.create(al_dir)}
  
  # Copy the alignment 
  file.copy(from = row$old_alignment_path, to = row$output_alignment_path)
  
  # Save the row as the parameters.csv file
  write.csv(row, file = row$parameters_path, row.names = FALSE)
}


format.timers <- function(time_csv){
  # A small function to read in a timing csv and read out a nice row of how long each step took
  
  # Open the csv file
  time_df <- read.csv(time_csv, stringsAsFactors = FALSE)
  # Check whether the unique identifier column is present
  if (("unique_id" %in% names(time_df)) == TRUE){
    # Remove row numbers column
    time_df <- time_df[,c("unique_id", "time_name", "timings")]
  } else if (("unique_id" %in% names(time_df)) == FALSE){
    # Extract the unique id from the file name and add it as a column to the dataframe
    time_df$unique_id <- gsub("\\.csv", "", gsub("_test_times", "", basename(time_csv)))
    # Remove row numbers column and rearrange column order
    time_df <- time_df[,c("unique_id", "time_name", "timings")]
  }
  # Remove row numbers column
  time_df <- time_df[,c("unique_id", "time_name", "timings")]
  # Change timings column to time format
  time_df$timings <- as.POSIXct(time_df$timings)
  
  # Find time for each test statistic
  test_cols <- c("uid", "likelihood_mapping", "scfs", "ntlt", "delta_plot", "q_residuals", "fast_tiger", "Cunningham_test", "tree_proportion", "total_time", "time_units")
  timings <- c(unique(time_df$unique_id)[[1]],
               format.one.time(time_df, "Start_likelihood_mapping", "End_likelihood_mapping"),
               format.one.time(time_df, "Start_scfs", "End_scfs"),
               format.one.time(time_df, "Start_ntlt", "End_ntlt"),
               format.one.time(time_df, "Start_delta_plot", "End_delta_plot"),
               format.one.time(time_df, "Start_q_residual", "End_q_residual"),
               format.one.time(time_df, "Start_fast_tiger", "End_fast_tiger"),
               format.one.time(time_df, "Start_Cunningham_test", "End_Cunningham_test"),
               format.one.time(time_df, "Start_tree_proportion", "End_tree_proportion"),
               format.one.time(time_df, "Start_time", "Results_saved_and_done"),
               "secs")
  
  # Assemble timings into a data frame row
  time_row <- as.data.frame(matrix(data = timings, nrow = 1, ncol = length(timings), byrow = TRUE))
  # Name the columns
  names(time_row) <- test_cols
  
  # Return the formatted time differences
  return(time_row)
}


format.one.time <- function(time_df, start_time_label, end_time_label){
  # Small function to take a start and end label (corresponding to objects within a column) and find the time difference in seconds between them
  
  # Find the time difference
  time_diff <- c(time_df[time_df$time_name == end_time_label, 3] - time_df[time_df$time_name == start_time_label, 3])[[1]]
  # Find what units the time difference is in
  time_unit <- attr(c(time_df[time_df$time_name == end_time_label, 3] - time_df[time_df$time_name == start_time_label, 3]), "units")
  
  # Convert the time difference to seconds
  if (time_unit == "secs"){
    # If already in seconds, do not change
    time_diff = time_diff
  } else if (time_unit == "mins"){
    # Convert minutes to seconds
    time_diff = time_diff * 60
  } else if (time_unit == "hours"){
    # Convert hours to seconds
    time_diff = time_diff * 60 * 60
  } else if (time_unit == "days"){
    # Convert days to seconds
    time_diff = time_diff * 60 * 60 * 24
  } else if (time_unit == "weeks"){
    # Convert days to seconds
    time_diff = time_diff * 60 * 60 * 24 * 7
  }
  
  # Round to 2dp
  time_diff <- round(time_diff, digits = 2)
  
  # Return the time difference
  return(time_diff)
}

