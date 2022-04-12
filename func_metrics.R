# /caitlinch/treelikeness_metrics/func_metrics.R
# Caitlin Cherryh 2022

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).



## Load required packages
library(ape) # for general tree/alignment wrangling, and the delta.plots function
library(ips) # to determine the indices of the parsimony informative sites before running TIGER
library(phangorn) # for splits and networks, for midpoint rooting trees 



#### Treelikeness metric functions ####
tree.proportion <- function(alignment_path, sequence_format = "DNA", model = "JC69", remove_trivial_splits = TRUE){
  ## Function to calculate the tree proportion: the proportion of split weights in the phylogenetic network captured by the minimum evolution tree
  
  ## Open alignment
  # Identify file type of alignment
  suffix <- tolower(tail(strsplit(alignment_path,"\\.")[[1]],1))
  
  # Open alignment
  if (suffix == "fa" | suffix == "fasta" | suffix == "fas" | suffix == "fna" | suffix == "faa" | suffix == "frn"){
    # Open alignment
    al <- read.FASTA(alignment_path, type = sequence_format)
  } else if (suffix == "nex" | suffix == "nexus") {
    al <- as.DNAbin(read.nexus.data(alignment_path))
  }
  
  ## NeighborNet network
  ## Estimate a NeighborNet network from the distance matrix and order splits from strongest to weakest
  # Compute pairwise distances for the taxa using the specified model of sequence evolution
  mldist <- dist.ml(al, model)
  # Create a NeighbourNet network from the alignment
  nnet <- neighborNet(mldist)
  
  ## Greedy tree using an adapted version of Kruskal's algorithm
  # Extract the splits from the NeighborNet network
  unordered_nw_splits <- as.splits(nnet)
  # Rearrange the splits in order from strongest to weakest (decreasing order by weight)
  nw_splits <- unordered_nw_splits[c(order(attr(unordered_nw_splits, "weight"), decreasing = TRUE))]
  # Find the list of compatible splits in the maximum weight tree
  # Add the first edge (i.e. the edge with the largest weight) to the set of edges comprising the maximum weight spanning tree
  compatible_splits <- c(1)
  # Iterate through each split. If the split is compatible with all other splits, add it to the set of edges comprising the maximum weight spanning tree
  # If the split is not compatible, it will add reticulation and the tree will become a network. Discard any non-compatible splits.
  for (i in 2:length(nw_splits)){
    # Assign the number of i to be the current split being tested for compatibility
    current_split <- i
    # Test whether the current split is compatible with all other splits that have been added
    test_compatibility <- is.compatible.bitsplits(as.bitsplits(nw_splits[c(compatible_splits, current_split)]))
    # If the split is compatible, add it to the list of compatible splits
    if (test_compatibility == TRUE){
      compatible_splits <- c(compatible_splits, current_split)
    }
  }
  # Take the tree as the set of compatible splits
  t_splits <- nw_splits[compatible_splits]
  
  ## Manage trivial splits
  # If requested, remove all trivial splits (using phangorn::removeTrivialSplits)
  if (remove_trivial_splits == TRUE){
    t_splits <- removeTrivialSplits(t_splits)
    nw_splits <- removeTrivialSplits(nw_splits)
  }
  
  ## Calculate tree proportion
  # Calculate the proportion of split weights included in the network are present in the tree
  t_split_weight_sum <- sum(attr(t_splits, "weight"))
  nw_split_weight_sum <- sum(attr(nw_splits, "weight"))
  tree_proportion <- t_split_weight_sum/nw_split_weight_sum
  
  ## Return result
  # Return the tree proportion value
  return(tree_proportion)
}



## Cunningham test (Cunningham 1978) ##
cunningham.test <- function(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", iqtree_substitution_model = "JC", distance_matrix_substitution_model = "JC69"){
  ## Function to estimate what proportion of the variance in the data is represented by the tree
  
  ## Test steps:
  # 1. Calculate the observed distances from the alignment (d_ij)
  # 2. Calculate the predicted distances from the tree (p_ij)
  # 3. Calculate the total sum of squares. TSS = sum of (d_ij)^2
  # 4. Calculate the residual sum of squares. RSS = sum of (p_ij - d_ij)^2
  # 5. Calculate the R^2 = (TSS - RSS)/ RSS
  
  ## Test:
  # 1. Calculate the observed distances (d_ij)
  dna <- read.dna(alignment_path, format = "fasta")
  dna_mat <- dist.ml(dna, model = distance_matrix_substitution_model) 
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



## Network tree-likeness test (Huson and Bryant 2006)
network.treelikeness.test <- function(alignment_path, splitstree_path, sequence_format = "DNA"){
  ## Uses Splitstree4.17.1 software to implement the Network Treelikeness Test described in Huson and Bryant (2006)
  # Software available from:
  # https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/
  
  ## Test steps:
  #   1. Infer a split network N from some data
  #         - Can be done in Splitstree using Networks -> NeighborNet
  #   2. Construct a confidence network for N with level 1 - alpha
  #         - Construct the set of splits with CIs excluding 0
  #         - Can be done in Splitstree4 using `Analysis -> Bootstrap...` then `Analysis -> Show Confidence Network`
  #   3. If the confidence network does not contain a tree, reject the null hypothesis that the data originated in a tree
  #         - Reject the null hypothesis if and only if this set is incompatible
  
  ## Construct and bootstrap a NeighborNet network
  # Convert fasta to nexus
  nexus_alignment_path <- convert.to.nexus(alignment_path, sequence_format = "DNA", include_taxablock = TRUE)
  # Name output path
  confidence_path <- paste0(alignment_path, "_confidence.nexus")
  output_path <- paste0(alignment_path, "_Splitstree_output.nex")
  # Run Splitstree4 if the confidence_path and output_path files do not exist
  if ((file.exists(confidence_path) == FALSE) | (file.exists(output_path) == FALSE)){
    # Assemble the SplitsTree 4 command
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", nexus_alignment_path, ";",
                                 " ASSUME chartransform=Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true;", 
                                 " ASSUME disttransform=NeighborNet;",
                                 " bootstrap runs=100;",
                                 " confidence_splits level=95 file=", confidence_path, ";",
                                 " export file=", output_path, " REPLACE=yes;",
                                 " quit;'")
    # Call SplitsTree 4
    system(splitstree_command)
  }
  
  ## Construct a confidence network using the bootstrap splits
  # Read in the nexus splits from the confidence network
  splits <- suppressWarnings(read.nexus.splits(confidence_path))
  
  ## Read in the splits from the confidence network and turn text file into a dataframe
  # Find the starting line for the splits (line after "MATRIX") and read in from that line down as tsv
  splits_start <- grep("MATRIX", readLines(confidence_path))
  splits_df <- read.delim(confidence_path, header = FALSE, skip = splits_start)
  # Remove rows that do not contain splits from data frame
  splits_df <- splits_df[1:length(splits),]
  # Format data frame columns
  names(splits_df) <- c("split_and_size", "weight", "interval", "taxa")
  splits_df$taxa <- gsub(",", "", splits_df$taxa)
  ss <- gsub(" ","",gsub("size=","",gsub("\\]","",gsub("\\[","",unlist(strsplit(splits_df$split_and_size, ","))))))
  splits_df$split <- as.numeric(ss[c(TRUE,FALSE)])
  splits_df$size <- as.numeric(ss[c(FALSE, TRUE)])
  si <- gsub(" ","",gsub("\\)","",gsub("\\(","",unlist(strsplit(splits_df$interval, ",")))))
  splits_df$interval_start <- as.numeric(si[c(TRUE, FALSE)])
  splits_df$interval_end <- as.numeric(si[c(FALSE, TRUE)])
  # Reorder data frame columns
  splits_df <- splits_df[,c("split","size","interval","interval_start","interval_end","taxa")]
  
  ## Construct the set of splits with confidence intervals excluding 0 (for the network treelikeness test)
  test_df <- splits_df[which(splits_df$interval_start > 0 & splits_df$interval_end > 0), ]
  test_splits <- splits[which(splits_df$interval_start > 0 & splits_df$interval_end > 0)]
  
  # A set of splits is compatible if all pairwise comparisons between splits are compatible
  # Check using the ape function is.compatible (requires splits to be in bitsplits format, which is also in ape)
  compatibility <- is.compatible.bitsplits(as.bitsplits(test_splits))
  
  # A set of splits is compatible if all pairwise comparisons between splits are compatible
  if (compatibility == FALSE){
    # Some pairwise comparisons between splits are incompatible: therefore, the null hypothesis that data was originated in a tree is rejected
    ntlt_result <- "Non-treelike"
  } else if (compatibility == TRUE){
    # All pairwise comparisons between splits are compatible: therefore, the null hypothesis that data was originated in a tree is accepted
    ntlt_result <- "Treelike"
  }
  
  ## Create an output vector for the results
  output_vector <- c(ntlt_result)
  names(output_vector) <- c("NetworkTreelikenessTest")
  
  ## Return Network Treelikeness Test results
  return(output_vector)
}



## Q-residuals (Gray et. al. 2010)
q_residuals <- function(alignment_path, phylogemetric_path, sequence_format = "DNA", phylogemetric_number_of_threads = NA){
  ## Uses software`phylogemetric` (available at https://github.com/SimonGreenhill/phylogemetric) to calculate
  #     the Q-residuals of Gray et. al. (2010)
  # Outputs a value between 0 and 1, where 0 indicates the quartets are treelike
  
  ## Change fasta file (from alisim) into nexus file
  # Nexus format required for `phylogemetric` software
  f <- read.FASTA(file = alignment_path, type = sequence_format)
  nexus_path <- paste0(alignment_path, ".nex")
  if (sequence_format == "DNA"){
    nexus_format = "dna"
  } else if ((sequence_format == "AA") | (sequence_format == "Protein")){
    nexus_format = "protein"
  }
  write.nexus.data(f, nexus_path, format = nexus_format, interleaved = FALSE)
  # Remove the "INTERLEAVE = FALSE" section (so other programs aren't confused)
  txt <- readLines(nexus_path)
  ind <- grep("FORMAT", txt)
  txt[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=-;"
  write(txt, file = nexus_path, append = FALSE)
  
  ## Apply `phylogemetric` and save output as a vector
  if (is.na(phylogemetric_number_of_threads)){
    # Run basic command
    call <- paste0(phylogemetric_path, " qresidual ", nexus_path)
  } else {
    # Use multiple cores with the -w/--workers argument
    call <- paste0(phylogemetric_path, " -w ", phylogemetric_number_of_threads, " qresidual ", nexus_path)
  }
  program_output <- system(call, intern = TRUE)
  # Format the program output into a nice vector containing just the mean Q-residual value per taxa
  q_residuals <- as.numeric(unlist(strsplit(program_output, "\t"))[c(FALSE, TRUE)])
  
  ## Return mean Q-residual value
  mean_q_residual <- mean(q_residuals)
  return(mean_q_residual)
}


## TIGER (Cummins and McInerney 2011)
TIGER <- function(alignment_path, fast_TIGER_path, sequence_format = "DNA"){
  # Tree Independent Genertion of Evolutionary Rates
  # Function to calculate TIGER values from a multiple sequence alignment using the 
  #   pbfrandsen/fast_TIGER software (available here: https://github.com/pbfrandsen/fast_TIGER)
  
  ## Remove phylogenetically uninformative sites from the alignment
  # Uninformative sites can bias the TIGER values - see List (2021)
  alignment <- read.FASTA(alignment_path, type = sequence_format)
  # Detect informative sites (the sites to keep)
  is_inds <- pis(as.matrix(alignment), what = "index")
  # Index the alignment, to reduce it to just the informative sites
  informative_alignment <- as.matrix.DNAbin(alignment)[,c(is_inds)]
  # Write the informative alignment to disk
  informative_alignment_path <- paste0(alignment_path, "_ParsimonyInformativeSites_only.phy")
  write.phy(informative_alignment, file = informative_alignment_path)
  
  ## Run fast_TIGER
  # Change sequence format to either "dna" or "aa" (only allowable data type names)
  if (sequence_format == "DNA"){
    data_type = "dna"
  } else if ((sequence_format == "AA") | (sequence_format == "Protein")){
    data_type = "protein"
  }
  # Create system command and call fast_TIGER
  call <- paste0(fast_TIGER_path, " ", data_type, " ", informative_alignment_path)
  system(call)
  
  ## Open the results file and output mean TIGER value
  # TIGER values range from 0 to 1, with values closer to 1 indicating more stable/consistent sites
  TIGER_file <- paste0(informative_alignment_path, "_r8s.txt")
  TIGER_values <- as.numeric(readLines(TIGER_file))
  # Calculate the mean TIGER value
  mean_TIGER_value <- mean(TIGER_values)
  # Return the mean value
  return(mean_TIGER_value)
}



## Likelihood mapping (Strimmer and von Haeseler 1997)
likelihood.mapping <- function(alignment_path, iqtree2_path, iqtree2_number_threads = 1, substitution_model = "MFP", number_of_taxa){
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
  
  ## Extract results from likelihood map
  iq_log <- readLines(iq_file)
  ind <- grep("Number of fully resolved  quartets",iq_log)
  resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  ind <- grep("Number of partly resolved quartets",iq_log)
  partly_resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  ind <- grep("Number of unresolved",iq_log)
  unresolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  total_q <- (resolved_q+partly_resolved_q+unresolved_q)
  prop_resolved <- resolved_q/total_q
  
  ## Collate results into a vector
  lm_results <- c(resolved_q, partly_resolved_q, unresolved_q, total_q, prop_resolved)
  names(lm_results) <- c("num_resolved_quartets", "num_partly_resolved_quartets", "num_unresolved_quartets",
                         "total_num_quartets", "proportion_resolved_quartets")
  
  ## Return results
  return(lm_results)
}



## Site concordance factors (Minh et. al. 2020)
scf <- function(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", number_scf_quartets = 100, substitution_model = "MFP", 
                add.likelihood.map = FALSE, number_of_taxa = NA){
  # Function to calculate the site concordance factors for an alignment, given a maximum likelihood tree estimated in IQ-Tree
  # Optional: can perform likelihood map (so tree doesn't have to be estimated multiple times to do both likelihood mapping and sCF)
  
  ## Check that the treefile already exists: if it doesn't, run IQ-Tree and create it
  if (add.likelihood.map == FALSE){
    if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
      # Given an alignment, estimate the maximum likelihood tree with IQ-Tree2
      call <- paste0(iqtree2_path," -s ",alignment_path," -m ",substitution_model," -nt ", iqtree2_number_threads, " -redo -safe")
      system(call)
    }
  } else if (add.likelihood.map == TRUE){
    if (((file.exists(paste0(alignment_path, ".iqtree")) == FALSE) | (file.exists(paste0(alignment_path, ".lmap.eps")) == FALSE)) &
        (is.na(number_of_taxa) == FALSE)) {
      number_of_quartets <- 25 * as.numeric(number_of_taxa)
      call <- paste0(iqtree2_path," -s ",alignment_path," -m ",substitution_model," -nt ", iqtree2_number_threads, " -lmap ",number_of_quartets," -redo -safe")
      system(call)
    }
  }
  ## Check if the site concordance factors have already been calculated: if they have not, calculate them
  if (file.exists(paste0(alignment_path,".treefile.cf.stat")) == FALSE){
    # Create the command and call it in the system
    # for sCF: iqtree -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -nt 10
    treefile <- paste0(alignment_path,".treefile")
    call <- paste0(iqtree2_path," -t ",treefile," -s ",alignment_path," --scf ",number_scf_quartets," -nt 1 -redo -safe")
    system(call)
  }
  ## Retrieve the site concordance factors from the output table
  scf_table <- read.table(paste0(alignment_path,".treefile.cf.stat"), header = TRUE, sep = "\t")
  scf_results <- list(mean_scf = round(mean(scf_table$sCF), digits = 2), 
                      median_scf = round(median(scf_table$sCF), digits = 2), 
                      all_scfs = scf_table$sCF, 
                      branch_ids = scf_table$ID)
  ## Return the site concordance factor results
  return(scf_results)
}



## Delta plots (Holland et. al. 2002)
mean.delta.plot.value <- function(alignment_path, sequence_format = "DNA", substitution_model = "JC69"){
  # This function takes an alignment, calculates a distance matrix for the alignment, and the applies the
  # `ape` function `delta.plot`. We take the mean delta plot value as the test statistic. 
  
  ## Open the alignment as a DNAbin object
  alignment <- read.FASTA(alignment_path, type = sequence_format)
  ## Calculate a distance matrix of pairwise distances from DNA sequences using a model of DNA substitution
  # Default model of DNA substitution is JC ("JC69") - it's used to simulate the sequences for the simulations
  pdm <- dist.dna(alignment, model = substitution_model)
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



#### Tree estimation functions ####
estimate.iqtree2.gene.trees <- function(gene_folder, iqtree2_path, iqtree2_number_threads = "AUTO", redo_flag = FALSE, safe_flag = FALSE, bootstraps = NA, model = "MFP"){
  ## Function to take a folder full of genes and estimate a gene tree for each one
  # Get the list of file names
  all_gene_paths <- paste0(gene_folder, list.files(gene_folder))
  all_gene_paths <- all_gene_paths[grep("\\.fa", all_gene_paths)]
  all_gene_paths <- all_gene_paths[grep("\\.fa\\.", all_gene_paths, invert = TRUE)]
  # Run IQ-Tree2 for each of those file names
  lapply(all_gene_paths, call.iqtree2, iqtree2_path, iqtree2_number_threads, redo_flag, safe_flag, bootstraps, model)
}

call.iqtree2<- function(gene_path, iqtree2_path, iqtree2_number_threads = "AUTO", redo_flag = FALSE, safe_flag = FALSE, bootstraps = NA, model = "MFP"){
  ## Small function to call IQ-Tree2 for one alignment
  # Use _flag commands from function call to assemble IQ-Tree2 call
  if (redo_flag == TRUE){
    redo_call = "-redo"
  } else if (redo_flag == FALSE){
    redo_call = ""
  }
  if (safe_flag == TRUE){
    safe_call = "-safe"
  } else if (safe_flag == FALSE){
    safe_call = ""
  }
  if (is.na(bootstraps)){
    bootstraps_call = ""
  } else if (is.na(bootstraps) == FALSE){
    bootstraps_call = bootstraps
  }
  # Assemble call
  call <- paste0(iqtree2_path, " -s ", gene_path, " -nt ", iqtree2_number_threads, " -m ", model, " ", redo_call, " ", safe_call, " ", bootstraps_call)
  print(call)
  # Invoke OS command
  system(call)
}

estimate.ASTRAL.species.tree <- function(gene_tree_file, species_tree_file, log_file, ASTRAL_path){
  ## Function to estimate a species tree using ASTRAL
  # Assemble ASTRAL command from input file names
  astral_command <- paste0("java -jar ", ASTRAL_path, " -i ", gene_tree_file, " -o ", species_tree_file, " 2> ", log_file)
  system(astral_command)
}

estimate.ASTRAL.multilocus.bootstrapping <- function(gene_tree_file, species_tree_file, log_file, ASTRAL_path, bootstraps_file){
  ## Function to estimate a species tree using ASTRAl and perform multi-locus bootstrapping
  # Assemble ASTRAL command from input file names
  astral_command <- paste0("java -jar ", ASTRAL_path, " -i ", gene_tree_file, " -b ", bootstraps_file, " -o ", species_tree_file, " 2> ", log_file)
  system(astral_command)
}



#### Network estimation functions ####
SPECTRE.estimate.network <- function(alignment_path, netmake_path, netme_path, sequence_format = "DNA"){
  ### Function to take an alignment, estimate a NeighborNet network in Netmake and estimate a minimum evolution spanning tree in NetME
  ### Requires SPECTRE software to run (both programs exist within SPECTRE)
  
  ## Name files
  suffix <- tolower(unlist(strsplit(basename(nexus_al_path), "\\."))[length(unlist(strsplit(basename(nexus_al_path), "\\.")))])
  al_name <- paste0(unlist(strsplit(basename(nexus_al_path), "\\."))[1:(length(unlist(strsplit(basename(nexus_al_path), "\\.")))-1)])
  alignment_dir <- dirname(nexus_al_path)
  
  ## Convert alignment to nexus (if it isn't already)
  if (suffix == "fasta" |suffix == "fa" | suffix == "fna" | suffix == "ffn" | suffix == "faa" | suffix == "frn" | suffix == "fas"){
    nexus_al_path <- convert.to.nexus(alignment_path, sequence_format, include_taxablock = FALSE)
  } else if (suffix == "nexus" | suffix == "nex"){
    nexus_al_path <- alignment_path
  }
  
  ## Set directory as alignment_dir
  setwd(alignment_dir)
  
  ## Run Netmake to create a compatible split system with a circular ordering from a distance matrix
  # Assemble netmake command
  netmake_command <- paste0(netmake_path, " -o ", al_name, " ", basename(nexus_al_path))
  # Run netmake command
  system(netmake_command)
  # Assemble name of netmake output file
  netmake_op_file <- paste0(alignment_dir, "/", al_name, ".network.nex")
  
  ## Run NetME to construct a minimum evolution tree from the specified split network with an implied circular order
  netme_command <- paste0(netme_path, " -o ", al_name, " ", basename(netmake_op_file))
  # Run NetME command
  system(netme_command)
  # Assemble names of NetME output files
  netme_nw_file <- paste0(alignment_dir, "/", al_name, ".min-evo.nex")
  netme_tl_file <- paste0(alignment_dir, "/", al_name, ".treelength")
  
  ## Return file names
  output_vec <- c(netmake_op_file, netme_nw_file, netme_tl_file)
  return(output_vec)
}



#### Functions to apply multiple test statistics ####
treelikeness.metrics.simulations <- function(replicate_folder, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
                                             num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                                             delta_plot_substitution_method = "JC69", num_phylogemetric_threads = NA, tree_proportion_substitution_method = "JC69",
                                             tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA"){
  ## Function to take one alignment, apply all treelikeness metrics and return results in a dataframe
  
  # Get alignment file
  folder_files <- list.files(replicate_folder)
  alignment_path <- paste0(replicate_folder, grep("output_alignment\\.fa\\.", grep("output_alignment\\.fa", folder_files, value = TRUE), invert = TRUE, value = TRUE))
  df_name <- gsub("output_alignment.fa", "treelikeness_tests.csv", alignment_path)
  
  # Check whether dataframe .csv file already exists. If it does, import the dataframe. If it doesn't, make it.
  if (file.exists(df_name) == TRUE){
    ## Read in the results csv file
    results_df <- read.csv(df_name)
  } else if (file.exists(df_name) == FALSE){
    ## Apply all treelikeness test statistics to generate the results csv file
    
    # Determine the number of taxa (needed for number of quartets in likelihood mapping and sCFs)
    if (grepl("exp1", replicate_folder)){
      random_trees_file <- paste0(replicate_folder, grep("random_trees", folder_files, value = TRUE))
      random_trees <- read.tree(random_trees_file)
      n_tree_tips <- unique(Ntip(random_trees))[[1]]
    } else if ((grepl("exp2", replicate_folder))|(grepl("exp3", replicate_folder))){
      start_tree_file <- paste0(replicate_folder, grep("starting_tree", folder_files, value = TRUE))
      start_tree <- read.tree(start_tree_file)
      n_tree_tips <- Ntip(start_tree)
    }
    
    # Apply Likelihood mapping (Strimmer and von Haeseler 1997)
    lm <- likelihood.mapping(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = iqtree_substitution_model, 
                             number_of_taxa = n_tree_tips)
    
    # Apply Site concordance factors (Minh et. al. 2020)
    scfs <- scf(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                substitution_model = iqtree_substitution_model, add.likelihood.map = FALSE, number_of_taxa = n_tree_tips)
    
    # Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(alignment_path, splitstree_path, sequence_format)
    
    # Apply Delta plots (Holland et. al. 2002)
    mean_delta_plot_value <- mean.delta.plot.value(alignment_path, sequence_format, substitution_model = delta_plot_substitution_method)
    
    # Apply Q-residuals (Gray et. al. 2010)
    mean_q_residual <- q_residuals(alignment_path, phylogemetric_path, sequence_format, phylogemetric_number_of_threads = num_phylogemetric_threads)
    
    # Apply TIGER (Cummins and McInerney 2011)
    mean_tiger_value <- TIGER(alignment_path, fast_TIGER_path, sequence_format)
    
    # Apply Cunningham test (Cunningham 1975)
    cunningham_metric <- cunningham.test(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", iqtree_substitution_model = "JC", distance_matrix_substitution_model = "JC69")
    
    # Apply tree proportion (new test)
    tree_proportion <- tree.proportion(alignment_path, sequence_format = sequence_format, model = tree_proportion_substitution_method,
                                       remove_trivial_splits = tree_proportion_remove_trivial_splits)
    
    ## Assemble results into a dataframe and save
    results_vec <- c(lm, scfs$mean_scf, scfs$median_scf, min(scfs$all_scfs), max(scfs$all_scfs), ntlt, mean_delta_plot_value, mean_q_residual, mean_tiger_value,
                     cunningham_metric, tree_proportion)
    results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
    names_vec <- c("LM_num_resolved_quartets", "LM_num_partly_resolved_quartets", "LM_num_unresolved_quartets", "LM_total_num_quartets", "LM_proportion_resolved_quartets",
                   "sCF_mean", "sCF_median", "sCF_min", "sCF_max", "NetworkTreelikenessTest", "mean_delta_plot_value", "mean_Q_residual", "mean_TIGER_value",
                   "Cunningham_test", "tree_proportion")
    names(results_df) <- names_vec
    write.csv(results_df, file = df_name, row.names = FALSE)
  }
  
  # Return results_df so it can be collated using lapply
  return(results_df)
}



#### Utility functions ####
convert.to.nexus <- function(alignment_path, sequence_format = "DNA", include_taxablock = FALSE){
  ### Convert fasta file to nexus file (if there is no existing nexus file with the same name)
  
  ## Prepare parameters for file conversion
  # Name nexus file by simply appending ".nex" to end of existing file name
  nexus_alignment_path <- paste0(alignment_path,".nex")
  # Extract file type from alignment path
  suffix <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  # Set format for output nexus file
  if ((sequence_format == "DNA") | (sequence_format == "dna")){
    nexus_format = "dna"
  } else if ((sequence_format == "Protein") | (sequence_format == "protein") | 
             (sequence_format == "AA") | (sequence_format == "aa")){
    nexus_format = "protein"
  }
  
  # Create a variable to specify whether to include a single DATA block (datablock = TRUE) or separate TAXA and CHARACTER boxes (datablock = FALSE)
  if (include_taxablock == TRUE){
    datablock_bool = FALSE
  } else if (include_taxablock == FALSE){
    datablock_bool = TRUE
  }
  
  ## Convert to nexus using functions based on suffix
  if (suffix == "fasta" |suffix == "fa" | suffix == "fna" | suffix == "ffn" | suffix == "faa" | suffix == "frn" | suffix == "fas"){
    ## If the file is a fasta file, convert it to nexus file format (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE){
      # Read in the fasta data
      data <- read.FASTA(alignment_path, type = sequence_format)
      # Write out the nexus data
      write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = datablock_bool, interleaved = FALSE) # write the output as a nexus file)
    }
  } else if (suffix == "phy" | suffix == "phylip"){
    ## If the file is a phy file, convert it to nexus file format (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE){
      data <- read.phy(alignment_path)
      write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = datablock_bool, interleaved = FALSE) # write the output as a nexus file)
    }
  }
  
  ## Open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(nexus_alignment_path)
  ind <- grep("BEGIN CHARACTERS",nexus)+2
  if ((sequence_format == "DNA") | (sequence_format == "dna")){
    nexus[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=DNA;"
  } else if ((sequence_format == "Protein") | (sequence_format == "protein") | 
             (sequence_format == "AA") | (sequence_format == "aa")){
    nexus[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN;"
  }
  # Write the edited nexus file out 
  writeLines(nexus,nexus_alignment_path)
  
  
  ## Output file name and path for nexus file
  return(nexus_alignment_path)
}



genes.from.alignment <- function(alignment_path, partition_path, gene_folder, sequence_format = "DNA", save.gene.info = FALSE){
  ## Function to take an alignment and a partition file, and save separate alignments for each gene (and a data frame with info about the genes)
  
  # Open the alignment as a matrix
  if (sequence_format == "DNA" | sequence_format == "dna"){
    aln_mat <- as.matrix(read.FASTA(alignment_path, type = "DNA"))
  } else if (sequence_format == "AA" | sequence_format == "aa" | sequence_format == "protein" | sequence_format == "Protein"){
    aln_mat <- as.matrix(read.FASTA(alignment_path, type = "AA"))
  }
  
  # Open partition file and extract genes
  partition_file <- readLines(partition_path)
  charpartition_inds <- grep("charset", partition_file)
  charpartition_lines <- partition_file[charpartition_inds]
  
  # Save each gene individually
  gene_info_list <- lapply(1:length(charpartition_lines), save.one.gene, charpartitions = charpartition_lines, alignment_matrix = aln_mat, output_folder = gene_folder)
  
  # If desired, save information about the genes
  if (save.gene.info == TRUE){
    gene_info_df <- do.call(rbind.data.frame, gene_info_list)
    names(gene_info_df) <- c("gene_name", "gene_start_position", "gene_end_position", "gene_length", "sequence_format")
    
    write.csv(gene_info_df, file = paste0(dirname(alignment_path), "/genes_from_alignment_info.csv"))
  }
}

save.one.gene <- function(index, charpartitions, alignment_matrix, output_folder, sequence_format = "DNA"){
  ## Small function to process a single charpartition row and save the associated gene
  
  # Extract the i^th row using the index and split the line to determine the gene name and location
  c_line <- charpartitions[index]
  c_line <- gsub(";", "", c_line)
  gene_name_split <- unlist(strsplit(c_line, "="))
  gene_name <- gsub("\tcharset ", "", gene_name_split[1])
  gene_name <- gsub(" ", "", gene_name)
  gene_range <- unlist(strsplit(c_line, ","))[2]
  gene_range <- gsub(" ", "", gene_range)
  gene_range_split <- strsplit(gene_range, "-")[[1]]
  gene_start <- as.numeric(gene_range_split[[1]])
  gene_end <- as.numeric(gene_range_split[[2]])
  gene_length = length(gene_start:gene_end)
  
  # Subset the matrix to include only the sites for gene 1 for all taxa
  # matrix[1:x, 1:y] where x is the number of taxa and y is the number of sites
  gene_alignment <- alignment_matrix[, gene_start:gene_end]
  
  # Name and output gene file
  gene_filepath <- paste0(output_folder, gene_name, ".fa")
  write.FASTA(gene_alignment, file = gene_filepath)
  
  # Assemble a vector of information about the gene
  gene_info <- c(gene_name, gene_start, gene_end, gene_length, sequence_format)
  names(gene_info) <- c("gene_name", "gene_start_position", "gene_end_position", "gene_length", "sequence_format")
  return(gene_info)
}




is.split.trivial <- function(split){
  ## Small function to determine whether a split is trivial or non-trivial (trivial = terminal branch)
  # Extract the bipartition subsets from the tree
  ss1 <- split[[1]] # get the indices of all taxa in the split
  taxa <- attr(split,"labels") # get the names of all taxa in the tree
  ss1_taxa <- taxa[ss1] # use the split indices to get the taxa names from the split
  ss2_taxa <- setdiff(taxa,ss1_taxa) # take the elements from taxa not in subset 1 to get the elements in ss2
  # Test whether the split is trivial (if it is trivial, either ss1 or ss1 will be 1)
  if ((length(ss1_taxa) == 1) || (length(ss2_taxa) == 1)){
    trivial <- TRUE
  } else {
    trivial <- FALSE
  }
  return(trivial)
}


trivial.splits.wrapper <- function(index, set_of_splits){
  ## Function that wraps around is.split.trivial for better application of lapply
  split <- set_of_splits[index]
  is.trivial <- is.split.trivial(split)
  return(is.trivial)
}



are.splits.compatible <- function(split1_index, split2_index, set_of_splits){
  ## Function to compare a pair of splits and determine whether they are compatible
  
  ## Open splits
  # A split "A|B" is a bipartition of a taxon set "X" into two non-empty sets
  split1 <- set_of_splits[split1_index]
  split2 <- set_of_splits[split2_index]
  # Get indices of all taxa in splits
  si1 <- split1[[1]]
  si2 <- split2[[1]]
  
  ## Determine whether splits are compatible
  #Get list of taxa for each side of each split
  taxa <- attr(set_of_splits, "labels")
  A1 <- taxa[si1]
  B1 <- setdiff(taxa, A1)
  A2 <- taxa[si2]
  B2 <- setdiff(taxa, A2)
  # Two splits S1 = A1|B1 and S2 = A2|B2 on X are compatible if one of the four possible intersections 
  #   of their split parts is empty: A1 ∩ A2, A1 ∩ B2, B1 ∩ A2, B1 ∩ B2
  #   Otherwise, the splits are incompatible
  #   A set of splits is compatible if all pairs of splits in S are compatible
  i1 <- intersect(A1, A2)
  i2 <- intersect(A1, B2)
  i3 <- intersect(B1, A2)
  i4 <- intersect(B1, B2)
  if ((length(i1) == 0) | (length(i2) == 0) | (length(i3) == 0) | (length(i4) == 0)){
    result = "Compatible"
  } else {
    result = "Incompatible"
  }
  
  ## Return output
  output = c(split1_index, split2_index, result)
  names(output) <- c("split1", "split2", "compatibility")
  return(output)
}


pairwise.compatibility <- function(index, set_of_splits){
  ## Function to take one split and compare pairwise to all other splits within a set to determine
  #     whether they are compatible
  
  # Get indices to compare split to all splits except itself
  comparison_indices <- setdiff(1:length(set_of_splits), index)
  # Compare split to all splits except itself
  pc_df <- do.call(rbind.data.frame, lapply(comparison_indices, are.splits.compatible, index, set_of_splits))
  names(pc_df) <- c("split1", "split2", "compatibility")
  # Return pairwise compatibility data frame
  return(pc_df)
}


remove.duplicate.splits <- function(df){
  ## Small function to remove duplicate splits from the compatibility dataframe
  
  # Process each row to nicely format the two splits being compared (to make comparison of rows easier)
  df$split_ids <- unlist(lapply(1:nrow(df), process.one.compatibility.matrix.row, df))
  # Remove any duplicates split_ids
  df <- df[!duplicated(df$split_ids), ]
  # Return dataframe with each comparison of two splits included only once
  return(df)
}


process.one.compatibility.matrix.row <- function(row_id, df){
  ## Small function to process one row of the compatibility dataframe and return a nice ordered string of the splits that are involved
  
  # Identify row and extract split ids for the two splits being compared in this row
  row <- df[row_id, ]
  split_ids <- c(as.numeric(row$split1), as.numeric(row$split2))
  # Combine the two split ids (small one first, larger one second) and separate with a column
  output_string <- paste0(min(split_ids), ",", max(split_ids))
  # Return the newly formatted string
  return(output_string)
}



