# caitlinch/treelikeness_metrics/code/func_empirical.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).

#### Required packages ####
library(ape)



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




#### Apply all treelikeness metrics ####
treelikeness.metrics.empirical <- function(alignment_path, splitstree_path, iqtree2_path, 
                                           iqtree_model = "MFP", num_iqtree2_threads = "AUTO",
                                           sequence_format = "AA", redo = FALSE){
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
                        substitution_model = scfl_model)
    
    ## Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(nexus_alignment_path, splitstree_path, sequence_format = sequence_format, nexus.file.format = TRUE)
    
    ## Apply Delta plots (Holland et. al. 2002)
    mldist_file <- paste0(alignment_path, ".mldist")
    mean_delta_plot_value <- delta.plot.empirical( dist_matrix = mldist.matrix(mldist_file) )
    
    ## Apply Cunningham test (Cunningham 1975)
    # Find the .mldist and .treefile output paths from IQ-Tree
    mldist_file <- paste0(alignment_path, ".mldist")
    tree_file <- paste0(alignment_path, ".treefile")
    # Calculate Cunningham metric
    cunningham_metric <- cunningham.test.empirical(mldist_file, tree_file)
    
    ## Apply tree proportion (new test)
    tree_proportion <- tree.proportion.long(nexus_alignment_path, sequence_format = sequence_format, model = NA, 
                                            remove_trivial_splits = TRUE, check_iqtree_log_for_identical_sequences = FALSE, 
                                            run_splitstree = TRUE, splitstree_path = splitstree_path)
    
    ## Assemble results into a dataframe and save
    results_vec <- c(unique_id, lm, scfl_output$mean_scf, scfl_output$median_scf, min(scfl_output$all_scfs), max(scfl_output$all_scfs), 
                     ntlt, mean_delta_plot_value, cunningham_metric, tree_proportion, alignment_path)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("unique_id", "LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets",
                   "LM_total_num_quartets", "LM_proportion_resolved_quartets", "sCF_mean", "sCF_median", "sCF_min", "sCF_max",
                   "NetworkTreelikenessTest", "mean_delta_plot_value", "Cunningham_test", "tree_proportion", "input_alignment_path")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE)
  } # end run treelikeness tests
  
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
calculate.p_value <- function(value_vector,id_vector){
  ## Function that given two vectors (one of test statistic values, and one of ids), calculates the p-value for that alignment
  # Create dataframe of id and values
  p_value_df <- data.frame(value_vector, id_vector, stringsAsFactors = FALSE)
  names(p_value_df) <- c("value","id")
  alignment_value <- p_value_df[which(p_value_df$id == "alignment"),1]
  if (is.na(alignment_value) == TRUE){
    # If the alignment value is NA, can't calculate a score
    # If it wasn't possible to calculate a score (will usually be a PHI score) for the alignment, output an NA
    p_value_2tail <- NA
  } else {
    # If it's possible to calculate a p-value, calculate one
    # Exclude NA rows
    p_value_df <-  p_value_df[!is.na(p_value_df$value),]
    # Find the number of bootstrap replicates and where the actual alignment value is located
    num_rows <- nrow(p_value_df) # number of bootstrap replicates + alignment value
    p_value_df <- p_value_df[order(p_value_df$value),] # order values from smallest to largest
    alignment_row <- which(p_value_df$id == "alignment") # find the ranking of the alignment value
    alignment_value <- p_value_df[alignment_row,1] # find the alignment's test statistic value
    # check whether there are other values that are the same as the alignment value
    identical_df <- subset(p_value_df,value == alignment_value)
    # if there are identical values, you don't know where the alignment actually falls within that list
    if (nrow(identical_df)>1){
      # get all the indexes of identical values
      identical_inds <- grep(alignment_value,p_value_df$value)
      # pick an ind at random
      random_identical_row <- sample(identical_inds,1)
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5 alignments less than or equal to the alignment value
      p_value_left <- random_identical_row/num_rows
    } else if (nrow(identical_df) == 1){
      # else, simply calculate the p value using the formula 
      # For left tail probability: want to find the number of observations less than or equal to the alignment value, then divide by the number of bootstrap observations
      # If there are 8 rows and the alignment is the 5th row, then there will be 5/8 alignments less than or equal to the alignment value
      p_value_left <- alignment_row/num_rows
    }
    # An alternative way to calculate this is to use the CDF:
    cdf <- ecdf(value_vector)
    p_value_cdf <- cdf(alignment_value)
  }
  
  # return the p-value
  op_vector <- c(p_value_left, p_value_cdf)
  names(op_vector) <- c("left_tail", "ecdf")
  return(op_vector)
}



