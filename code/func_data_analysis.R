# /caitlinch/treelikeness_metrics/func_metrics.R
# Caitlin Cherryh 2022

# This file contains functions for data manipulation and data analysis

collate.treelikeness.results <- function(alignment_path){
  # Function to take one alignment, get the parameters csv and treelikeness results csv, glue them together and return the new csv as a data.frame
  
  # Get paths to csv files
  alignment_params_file <- gsub("_output_alignment.fa", "_parameters.csv", alignment_path)
  alignment_treelikeness_file <- gsub("_output_alignment.fa", "_treelikeness_results.csv", alignment_path)
  
  # Check if files exist
  if (file.exists(alignment_params_file) == TRUE & file.exists(alignment_treelikeness_file) == TRUE){
    # Open the two csv files
    p_df <- read.csv(alignment_params_file)
    tl_df <- read.csv(alignment_treelikeness_file)
    
    # Remove some columns from the p_df
    p_df <- p_df[,c("row_id", "uid" , "num_taxa", "num_trees", "tree_depth",  "num_reps", "alisim_gene_models", "alisim_gene_tree_length", "total_alignment_length", "sequence_type")]
    
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




