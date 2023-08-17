# caitlinch/treelikeness_metrics/code/func_empirical.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).



#### Wrapper functions ####
empirical.treelikeness.test.wrapper <- function(alignment_path, 
                                                iqtree2_path, splitstree_path, 
                                                sequence_format = "AA", 
                                                num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                                                iqtree_substitution_model = "LG", distance_matrix_substitution_method = "LG", 
                                                tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = FALSE, 
                                                redo = FALSE, number_parallel_cores = 1){
  ### Function to apply treelikeness metrics to an empirical dataset, perform a parametric bootstrap with 100 replicates, 
  #       and return the output and p-values for each metric
  
  ## Estimate IQ-Tree with MFP
  # Call IQ-Tree with MFP to determine the best model for this alignment
  iqtree_run <- paste0(iqtree2_path, " -s ", alignment_path, " -m MFP -nt ", num_iqtree2_threads)
  system(iqtree_run)
  # Extract the treefile path
  tree_path <- paste0(alignment_path, ".treefile")
  
  ## Extract parameters from the IQ-Tree run of the alignment
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Get IQ-Tree file names
  iqtree_file <- paste0(alignment_path, ".iqtree")
  log_file <- paste0(alignment_path, ".log")
  
  ## Apply treelikeness metric to alignment
  treelikeness.metrics.empirical(alignment_path = alignment_path, 
                                 iqtree2_path = iqtree2_path, 
                                 splitstree_path = iqtree2_path, 
                                 num_iqtree2_threads = num_iqtree2_threads, 
                                 num_iqtree2_scf_quartets = num_iqtree2_scf_quartets, 
                                 iqtree_substitution_model = iqtree_substitution_model, 
                                 distance_matrix_substitution_method = distance_matrix_substitution_method, 
                                 tree_proportion_remove_trivial_splits = tree_proportion_remove_trivial_splits, 
                                 run_splitstree_for_tree_proportion = run_splitstree_for_tree_proportion, 
                                 sequence_format = sequence_format, 
                                 redo = redo)
  
  ## Generate 100 replicate alignments
  sim_al_prefix <- "bs_rep_"
  generate.mimic.alignment(output_prefix = paste0(replicate_folder, sim_al_prefix), 
                           alignment_path = alignment_path, 
                           iqtree_path = iqtree2_path, 
                           substitution_model = "MFP", 
                           output_format = "fasta", 
                           sequence_type = "AA",
                           num_output_alignments = 100)
  
  ## Perform 100 parametric bootstrap replicates
  # Extract all simulated alignments
  all_files <- list.files(replicate_directory())
  # Run treelikeness tests
  lapply(alignment_path,  treelikeness.metrics.empirical,
         iqtree2_path = iqtree2_path, 
         splitstree_path = iqtree2_path, 
         num_iqtree2_threads = num_iqtree2_threads, 
         num_iqtree2_scf_quartets = num_iqtree2_scf_quartets, 
         iqtree_substitution_model = iqtree_substitution_model, 
         distance_matrix_substitution_method = distance_matrix_substitution_method, 
         tree_proportion_remove_trivial_splits = tree_proportion_remove_trivial_splits, 
         run_splitstree_for_tree_proportion = run_splitstree_for_tree_proportion, 
         sequence_format = sequence_format, 
         redo = redo,
         mc.cores = number_parallel_cores)
  
  ## Create a nice dataframe of all the output values
  

## Calculate p-values

## Assemble output csv
}




#### Apply all treelikeness metrics ####
treelikeness.metrics.empirical <- function(alignment_path, 
                                           iqtree2_path, splitstree_path, 
                                           sequence_format = "AA", num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                                           iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
                                           tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = FALSE, 
                                           redo = FALSE){
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
    scfl_output <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                        substitution_model = iqtree_substitution_model)
    
    # Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(nexus_alignment_path, splitstree_path, sequence_format = sequence_format, nexus.file.format = TRUE)
    
    # Apply Delta plots (Holland et. al. 2002)
    mean_delta_plot_value <- delta.plot.empirical(alignment_path, sequence_format = sequence_format, substitution_model = distance_matrix_substitution_method)
    
    # Apply Cunningham test (Cunningham 1975)
    cunningham_metric <- cunningham.test.empirical(alignment_path, sequence_format, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, 
                                                   iqtree_substitution_model = iqtree_substitution_model, 
                                                   distance_matrix_substitution_model = distance_matrix_substitution_method)
    
    # Apply tree proportion (new test)
    tree_proportion <- tree.proportion.long(nexus_alignment_path, sequence_format = sequence_format, model = distance_matrix_substitution_method, 
                                            remove_trivial_splits = tree_proportion_remove_trivial_splits, check_iqtree_log_for_identical_sequences = FALSE, 
                                            run_splitstree = run_splitstree_for_tree_proportion, splitstree_path = splitstree_path)
    
    # Assemble results into a dataframe and save
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



cunningham.test.empirical <- function(alignment_path, sequence_format = "AA", iqtree2_path, iqtree2_number_threads = "AUTO", 
                                      iqtree_substitution_model = "LG", distance_matrix_substitution_model = "LG"){
  ## Function to estimate what proportion of the variance in the data is represented by the tree
  
  ## Test steps:
  # 1. Calculate the observed distances from the alignment (d_ij)
  # 2. Calculate the predicted distances from the tree (p_ij)
  # 3. Calculate the total sum of squares. TSS = sum of (d_ij)^2
  # 4. Calculate the residual sum of squares. RSS = sum of (p_ij - d_ij)^2
  # 5. Calculate the R^2 = (TSS - RSS)/ RSS
  
  ## Test:
  # 1. Calculate the observed distances (d_ij)
  p <- phyDat(read.FASTA(alignment_path, type = sequence_format), type = sequence_format)
  dna_mat <- dist.ml(p, model = distance_matrix_substitution_model)
  d_ij <- as.vector(dna_mat) # observed distances between taxa i and j
  
  # 2. Calculate the predicted distances (p_ij)
  # Infer a tree (if one does not already exist) and open it
  if (file.exists(paste0(alignment_path, ".treefile")) == FALSE){
    call.iqtree2(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", redo_flag = FALSE, safe_flag = FALSE, bootstraps = NA, model = iqtree_substitution_model)
  }
  # Open the tree
  t <- read.tree(paste0(alignment_path, ".treefile"))
  # Extract the distance matrix from the tree
  t_cophenetic_mat <- cophenetic.phylo(t) # in substitutions per site
  # Now reorder the t_mat so the taxa are in the same order
  dna_order <- attr(dna_mat, "Labels")
  t_ordering_mat <- as.matrix(t_cophenetic_mat)[dna_order, dna_order]
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



#### Model extraction functions ####
get.simulation.parameters <- function(dotiqtree_file){
  # Given a .iqtree file, this function will extract the relevant parameters
  
  # read in the IQ-TREE file to get substitution model and parameters
  iq_file <- readLines(dotiqtree_file)
  # extract the file name
  ind      <- grep("Input file name:",iq_file)
  op1      <- substr(iq_file[[ind]],18,nchar(iq_file[[ind]]))
  # extract the number of taxa and extract the length of the alignment
  ind         <- grep("Input data:",iq_file)
  input_str   <- iq_file[[ind]] # get the line that contains this info
  input_ls    <- strsplit(input_str," ")
  op2         <- input_ls[[1]][3] # extract number of sequences (number of taxa)
  op3         <- input_ls[[1]][6] # extract number of sites 
  # Extract the model of substitution (same for amino-acid and nucleotide files)
  ind         <- grep("Model of substitution: ",iq_file)
  sub_str     <- iq_file[[ind]] # get the line that contains this info
  sub_ls      <- strsplit(sub_str," ")
  op4         <- sub_ls[[1]][4]
  # Extract information about the sequence alignment
  ind <- grep("Number of constant sites:",iq_file)
  num_lines <- iq_file[c(ind:(ind+3))] # take the four lines listing the number of different kinds of sites
  num_split <- unlist(strsplit(num_lines,":")) # Split the lines at the colon
  num_names <- num_split[c(TRUE,FALSE)] # before the colon = name
  num_vals <- num_split[c(FALSE,TRUE)] # after the colon = value
  num_vals_regx <- regmatches(num_vals, gregexpr("[[:digit:]]+", num_vals)) # extract all the numbers after the colon
  # four strings = four lists of numbers: take first value from each list to get number of sites
  num_vals <- c(num_vals_regx[[1]][1],num_vals_regx[[2]][1],num_vals_regx[[3]][1],num_vals_regx[[4]][1]) 
  
  # check the type of sites - amino acid or DNA
  # if the input is DNA (nucleotide sites), gather that information
  if (input_ls[[1]][7]=="nucleotide"){
    # Extract the rate parameters
    rate1 <- as.numeric(strsplit(iq_file[[grep("A-C",iq_file)]],":")[[1]][2]) # A-C rate (same as code above, but combined 4 lines into 1 line)
    rate2 <- as.numeric(strsplit(iq_file[[grep("A-G",iq_file)]],":")[[1]][2]) # A-G rate
    rate3 <- as.numeric(strsplit(iq_file[[grep("A-T",iq_file)]],":")[[1]][2]) # A-T rate
    rate4 <- as.numeric(strsplit(iq_file[[grep("C-G",iq_file)]],":")[[1]][2]) # C-G rate
    rate5 <- as.numeric(strsplit(iq_file[[grep("C-T",iq_file)]],":")[[1]][2]) # C-T rate
    rate6 <- as.numeric(strsplit(iq_file[[grep("G-T",iq_file)]],":")[[1]][2]) # G-T rate
    
    # Extract the state frequencies
    state_freq_line <- iq_file[[grep("State frequencies",iq_file)]]
    if (state_freq_line == "State frequencies: (equal frequencies)"){
      # If the state frequencies are all equal, assign them all to 0.25 (1/4)
      sf1 <- 0.25 # pi(A) - A freq.
      sf2 <- 0.25 # pi(C) - C freq.
      sf3 <- 0.25 # pi(G) - G freq.
      sf4 <- 0.25 # pi(T) - T freq.
    } else {
      # If the state frequencies are not all equal, extract what they are
      sf1 <- as.numeric(strsplit(iq_file[[grep("pi\\(A\\)",iq_file)]],"=")[[1]][2]) # pi(A) - A freq. Remember to double backslash to escape before brackets
      sf2 <- as.numeric(strsplit(iq_file[[grep("pi\\(C\\)",iq_file)]],"=")[[1]][2]) # pi(C) - C freq.
      sf3 <- as.numeric(strsplit(iq_file[[grep("pi\\(G\\)",iq_file)]],"=")[[1]][2]) # pi(G) - G freq.
      sf4 <- as.numeric(strsplit(iq_file[[grep("pi\\(T\\)",iq_file)]],"=")[[1]][2]) # pi(T) - T freq.
    }
    
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    mrh2_name <- gsub(" ","_",mrh2_name) # change the name to be easy to parse
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model", "A-C_rate","A-G_rate","A-T_rate","C-G_rate","C-T_rate","G-T_rate","A_freq","C_freq",
               "G_freq","T_freq","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name","model_of_rate_heterogeneity_line2_value")
    # Make a list of the output rows for the output dataframe
    op <- c(op1,"DNA",op2,op3,num_vals,op4,rate1,rate2,rate3,rate4,rate5,rate6,sf1,sf2,sf3,sf4,mrh1,mrh2_name,mrh2)
    # Create the output dataframe
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    # Name the columns
    names(op_df) <- c("parameter","value")
    
    # Create the rate matrix Q
    Q_start <- grep("Rate matrix Q:",iq_file)+2
    Q_end   <- Q_start+3
    # Create the columns
    c1 <- c("A","C","G","T")
    c2 <- c()
    c3 <- c()
    c4 <- c()
    c5 <- c()
    # For each row in the iqtree file rate matrix
    for (i in Q_start:Q_end){
      # Split the row
      row <- strsplit(iq_file[[i]]," ")[[1]]
      row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
      # Add the resulting values to the relevant columns
      c2 <- c(c2,as.numeric(row[1])) # convert to numeric so can use the numbers more easily later
      c3 <- c(c3,as.numeric(row[2]))
      c4 <- c(c4,as.numeric(row[3]))
      c5 <- c(c5,as.numeric(row[4]))
    }
    # Create a dataframe of the rate matrix Q
    q_df <- data.frame(c1,c2,c3,c4,c5, stringsAsFactors = FALSE)
    #Rename the columns
    names(q_df) <- c("nucleotide","A","C","G","T")
    
    # Check if the model for rate heterogeneity is uniform
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      # If the model isn't uniform, need to create a matrix to collect and store the gamme category information
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # Start collecting info for the matrix
      g1 <- c() # initialise columns to store data in
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]]," ")[[1]] # split the rows on the (large amount of) " "'s in the middle
        row <- row[str_detect(row,"([0-9])")] # take only the numeric elements of the vector
        g1 <- c(g1,as.numeric(row[1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[2]))
        g3 <- c(g3,as.numeric(row[3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns
    }
    
    # Create a list of the three dataframes
    # This will be the output 
    params <- list(op_df,g_df,q_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories","Q_rate_matrix")
    
  } else if (input_ls[[1]][7]=="amino-acid"){ # alternatively if the data is amino acid sites
    # Extract model of rate heterogeneity
    mrh1      <- strsplit(iq_file[[grep("Model of rate heterogeneity:",iq_file)]],":")[[1]][2] # Extract model of rate heterogeneity 
    mrh2      <- as.numeric(strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][2]) # Line after the "model of rate heterogeneity" varies - extract it regardless of what it is 
    mrh2_name <- strsplit(iq_file[[(grep("Model of rate heterogeneity:",iq_file)+1)]],":")[[1]][1] # As the line varies, extract the name for the output dataframe
    # Extract state frequencies
    sf1      <- strsplit(iq_file[[grep("State frequencies:",iq_file)]],":")[[1]][2]
    
    # make a list of the rows for the output dataframe
    names <- c("file_name","sequence_type","n_taxa","n_sites",num_names,"substitution_model","model_of_rate_heterogeneity","model_of_rate_heterogeneity_line2_name",
               "model_of_rate_heterogeneity_line2_value","state_frequencies")
    # Make a list of the output rows for the first output dataframe
    op <- c(op1,"amino-acid",op2,op3,num_vals,op4,mrh1,mrh2_name,mrh2,sf1)
    op_df <- data.frame(names,op, stringsAsFactors = FALSE)
    names(op_df) <- c("parameter","value")
    
    # Check whether a gamma matrix is needed
    mrh1_check <- gsub(" ","",mrh1)
    if (mrh1_check=="Uniform"){
      # If the model for rate heterogeneity is uniform, don't need to create a matrix for discrete gamma rate categories
      g_df <- "Uniform"
    } else {
      #Create the matrix for discrete gamma categories
      g_start <- grep(" Category",iq_file)+1 # get the index for the first line of the gamma categories matrix
      empty   <- which(iq_file=="") # get indexes of all empty lines
      empty   <- empty[empty>g_start] # get empty lines above gamma categories matrix
      g_end   <- empty[1]-1 # get end index for gamma categories matrix (one less than next empty line)
      end_line <- iq_file[g_end]
      # if the end isn't an empty line, subtract one from the end count 
      # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
      # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
      check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1]," ")[[1]])
      if (check_line > 3){
        # If the check_line is longer than 3 objects, it won't be a group for the gamma categories but an instruction
        # Eg. the relative rates line has 17 objects, because each word in the sentence is an object
        # Instructions can be excluded from the gamma matrix (but categories can't)
        g_end = g_end - 1
      }
      # initialise columns to store data in
      g1 <- c()
      g2 <- c()
      g3 <- c()
      # Iterate through rows in gamma matrix
      for (i in g_start:g_end){
        row <- strsplit(iq_file[[i]],"        ") # split the rows on the long strong of 0's in the middle
        g1 <- c(g1,as.numeric(row[[1]][1])) # add the values to the columns
        g2 <- c(g2,as.numeric(row[[1]][2]))
        g3 <- c(g3,as.numeric(row[[1]][3]))
      }
      g_df <- data.frame(g1,g2,g3, stringsAsFactors = FALSE) # create a dataframe of the information
      names(g_df) <- c("category","relative_rate","proportion") # name the columns 
    }
    
    # Check whether state frequencies are needed
    sf1_squashed <- gsub(" ","",sf1)
    if (sf1_squashed == "(empiricalcountsfromalignment)"){
      # Get starting line for frequencies
      start_ind <- grep("State frequencies:",iq_file) + 2
      # Take the 20 lines containing AA frequencies
      freq_lines <- iq_file[start_ind:(start_ind+19)]
      # Split up the frequency lines into the label and the frequency
      freq_split <- unlist(strsplit(freq_lines,"="))
      # Get the frequency
      freq_nums <- freq_split[c(FALSE,TRUE)]
      # Remove any spaces (from IQTree formatting)
      freq_nums <- gsub(" ","",freq_nums)
      # Get corresponding AA letter
      freq_names <- freq_split[c(TRUE,FALSE)]
      # Remove IQTree formatting
      freq_names <- gsub("pi\\(","",freq_names)
      freq_names <- gsub("\\)","",freq_names)
      freq_names <- gsub(" ","", freq_names)
      # Create a nice dataframe
      f_df <- data.frame("amino_acid" = freq_names,
                         "frequency" = freq_nums,
                         stringsAsFactors = FALSE)
    } else {
      f_df <- "State frequencies from model"
    }
    
    # Create a list of the dataframes - this will be the output
    params <- list(op_df,g_df,f_df)
    # Name the parameters so they're easy to access once you've outputted the data
    names(params) <- c("parameters","gamma_categories", "frequency")
  }
  
  # Now the information has been collected, create an output dataframe
  # make the output dataframe
  return(params)
}


