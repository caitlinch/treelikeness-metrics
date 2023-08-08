# caitlinch/treelikeness_metrics/code/func_metrics.R
# Caitlin Cherryh 2023

# This file contains functions to apply tests for treelikeness to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).



## Load required packages
library(ape) # for general tree/alignment wrangling, and the delta.plots function
library(ips) # to determine the indices of the parsimony informative sites before running TIGER
library(phangorn) # for splits and networks, for midpoint rooting trees 



#### Treelikeness metric functions ####
tree.proportion.long <- function(alignment_path, sequence_format = "DNA", model = "JC69", remove_trivial_splits = TRUE, 
                                 check_iqtree_log_for_identical_sequences = FALSE, run_splitstree = FALSE, splitstree_path = NA, 
                                 base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA){
  ## Function to calculate the tree proportion: the proportion of split weights in the phylogenetic network captured by the minimum evolution tree
  
  ## Check whether multiple sequences in the alignment are identical
  if (check_iqtree_log_for_identical_sequences == TRUE) {
    # Check whether this alignment contains any identical sequences
    identical_check <- check.iqtree.log.for.identical.sequences(alignment_path, sequence_format = sequence_format)
    identical_sequences_present <- as.logical(identical_check[["identical_sequences_present"]])
    
    # If there are any identical sequences, do not run the tree proportion code
    if (identical_sequences_present == TRUE) {
      # If there are identical sequences, do not run tree proportion code: will not estimate network nicely
      run_tree_proportion = FALSE
    } else if (identical_sequences_present == FALSE) {
      # If there are no identical sequences, run tree proportion code
      run_tree_proportion = TRUE
    }
    
  } else  if (check_iqtree_log_for_identical_sequences == FALSE) {
    # If not checking for identical sequences in IQ-Tree log file, just proceed to tree proportion code
    run_tree_proportion = TRUE
    identical_sequences_present = NA
  }
  
  ## Run the tree proportion metric
  if (run_tree_proportion == TRUE) {
    if (run_splitstree == FALSE) {
      ## Calculate NeighborNet network in R
      
      # Estimate a NeighborNet network from the distance matrix and order splits from strongest to weakest
      # Compute pairwise distances for the taxa using the specified model of sequence evolution
      mldist <- calculate.dna.pairwise.distance.matrix(alignment_path, sequence_format, model, base_frequencies, Q_matrix, number_of_rate_categories)
      # Create a NeighbourNet network from the alignment
      nnet <- neighborNet(mldist)
      
      # Extract the splits from the NeighborNet network
      unordered_nw_splits <- as.splits(nnet)
    } else if (run_splitstree == TRUE & is.na(splitstree_path) == FALSE) {
      ## Calculate NeighborNet network using Splitstree 
      unordered_nw_splits <- make.splitstree.neighbornet(alignment_path, splitstree_path, return.splits = TRUE)
    }
    
    ## Identify splits in greedy tree
    if (length(unordered_nw_splits) == 1) {
      # If only one split in network, then it is compatible with itself. Therefore the set of splits in the tree is identical to the set of splits in the network
      nw_splits <- unordered_nw_splits
      t_splits <- nw_splits
    } else if (length(unordered_nw_splits) > 1) {
      ## Greedy tree using an adapted version of Kruskal's algorithm (starting with unordered splits from NeighborNet network)
      # Rearrange the splits in order from strongest to weakest (decreasing order by weight)
      nw_splits <- unordered_nw_splits[c(order(attr(unordered_nw_splits, "weight"), decreasing = TRUE))]
      # Find the list of compatible splits in the maximum weight tree
      # Add the first edge (i.e. the edge with the largest weight) to the set of edges comprising the maximum weight spanning tree
      compatible_splits <- c(1)
      # Iterate through each split. If the split is compatible with all other splits, add it to the set of edges comprising the maximum weight spanning tree
      # If the split is not compatible, it will add reticulation and the tree will become a network. Discard any non-compatible splits.
      for (i in 2:length(nw_splits)) {
        # print(i)
        # i is the current split being tested for compatibility
        # Test whether the current split is compatible with all other splits that have been added
        compatibility_matrix <- as.matrix(compatible(nw_splits))[c(compatible_splits, i), c(compatible_splits, i)]
        # print(compatibility_matrix)
        # Sum up all the values in the commpatibility matrix
        check_sum = sum(compatibility_matrix)
        # Check compatibility (by checking the sum of the matrix values)
        if (check_sum == 0){
          # If the matrix values sum to 0, then the set of splits are compatible
          check_compatibility = TRUE
        } else if (check_sum > 0){
          check_compatibility = FALSE
        }
        # print(check_compatibility)
        # If the split is compatible, add it to the list of compatible splits
        if (check_compatibility == TRUE){
          # If there are any incompatible splits in the matrix, they are represented by a 1
          # Therefore is the sum of matrix values is larger than 0, the added split is incompatible with the existing set of splits
          compatible_splits <- c(compatible_splits, i)
        }
      }
      # Take the tree as the set of compatible splits
      t_splits <- nw_splits[compatible_splits]
    }
    
    ## Remove trivial splits (if specified), then calculate tree proportion
    if (remove_trivial_splits == TRUE) {
      # Check for presence and number of trivial splits
      check_t <- trivial.splits.present(t_splits)
      check_nw <- trivial.splits.present(nw_splits)
      # If present, remove trivial splits in tree
      if (check_t$TrivialSplitsPresent == TRUE & check_t$Num_non_trivial_splits > 0) {
        # If one or more trivial splits present in tree, remove all trivial splits (using phangorn::removeTrivialSplits)
        t_splits <- removeTrivialSplits(t_splits)
      } else if (check_t$TrivialSplitsPresent == TRUE & check_t$Num_non_trivial_splits == 0) {
        # If only one split is present in tree and it is trivial, set t_splits to NA
        t_splits <- NA
      }
      # If present, remove trivial splits in network
      if (check_nw$TrivialSplitsPresent == TRUE & check_nw$Num_non_trivial_splits > 0){
        # If one or more trivial splits present in network, remove all trivial splits (using phangorn::removeTrivialSplits)
        nw_splits <- removeTrivialSplits(nw_splits)
      } else if (check_nw$TrivialSplitsPresent == TRUE & check_nw$Num_non_trivial_splits == 0) {
        # If only one split is present in network and it is trivial, set nw_splits to NA
        nw_splits <- NA
      }
      ## Calculate tree proportion
      ## Check class: can only calculate tree proportion if both t_splits and nw_splits are class "splits"
      if (class(t_splits) != "splits" & class(nw_splits) != "splits") {
        # If both tree and network are NA, that means there were no non-trivial splits in the alignment
        # Return that value
        tree_proportion <- "Only_trivial_splits_in_network_and_tree"
      } else if (class(t_splits) == "splits" & class(nw_splits) == "splits") {
        # If neither tree or network are NA, that means there were non-trivial splits in both
        ## Calculate tree proportion
        t_split_weight_sum <- sum(attr(t_splits, "weight"))
        nw_split_weight_sum <- sum(attr(nw_splits, "weight"))
        tree_proportion <- t_split_weight_sum/nw_split_weight_sum
      } else if (class(t_splits) != "splits" & class(nw_splits) == "splits") {
        # If only trivial splits in tree and other splits in network, means no splits from network are included in tree
        # Therefore tree proportion must be 0 (no trees in network also in tree: tree proportion = 0/x = 0)
        tree_proportion <- "0_no_splits_in_tree"
      } else if (class(t_splits) != "splits" & class(nw_splits) != "splits") {
        # If only trivial splits in network and not in tree, something must have gone wrong
        # Report that
        tree_proportion <- "0_no_splits_in_network"
      }
    } else if (remove_trivial_splits == FALSE) {
      ## Do NOT check for trivial splits
      ## Calculate tree proportion
      # Calculate the proportion of split weights included in the network are present in the tree
      t_split_weight_sum <- sum(attr(t_splits, "weight"))
      nw_split_weight_sum <- sum(attr(nw_splits, "weight"))
      tree_proportion <- t_split_weight_sum/nw_split_weight_sum
    }
  } else if (run_tree_proportion == FALSE) {
    if (identical_sequences_present == TRUE) {
      # If there are identical sequences, the neighborNet function will not work nicely to calculate tree proportion
      # Report that information
      tree_proportion = "identical_seqs_no_tree_proportion"
    } else {
      tree_proportion = "did_not_run"
    }
  }
  
  ## Return result
  # Return the tree proportion value
  return(tree_proportion)
}



## Cunningham test (Cunningham 1978) ##
cunningham.test <- function(alignment_path, sequence_format = "DNA", iqtree2_path, iqtree2_number_threads = "AUTO", iqtree_substitution_model = "JC", distance_matrix_substitution_model = "JC69",
                            base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA){
  ## Function to estimate what proportion of the variance in the data is represented by the tree
  
  ## Test steps:
  # 1. Calculate the observed distances from the alignment (d_ij)
  # 2. Calculate the predicted distances from the tree (p_ij)
  # 3. Calculate the total sum of squares. TSS = sum of (d_ij)^2
  # 4. Calculate the residual sum of squares. RSS = sum of (p_ij - d_ij)^2
  # 5. Calculate the R^2 = (TSS - RSS)/ RSS
  
  ## Test:
  # 1. Calculate the observed distances (d_ij)
  dna_mat <- calculate.dna.pairwise.distance.matrix(alignment_path, sequence_format, distance_matrix_substitution_model, base_frequencies, Q_matrix, number_of_rate_categories)
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
  
  ## Determine Network Treelikeness Test result
  # First, determine which splits have confidence intervals excluding 0
  which_splits <- which(splits_df$interval_start > 0 & splits_df$interval_end > 0)
  # Second, use that set to determine the Network Treelikeness Test result
  if (length(which_splits) == 0){
    ntlt_result <- "Zero_splits_where_confidence_intervals_exclude_0"
  } else if (length(which_splits) == 1){
    # If there is only one split with confidence intervals excluding 0:
    # The set is compatible by default
    ntlt_result <- "Treelike"
  } else if (length(which_splits) > 1){
    # If there is one or more split with confidence intervals excluding 0:
    # Construct the set of splits with confidence intervals excluding 0 (for the network treelikeness test)
    test_df <- splits_df[which_splits, ]
    test_splits <- splits[which_splits]
    
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
  taxa_order <- gsub(" ", "", unlist(strsplit(program_output, "\t"))[c(TRUE, FALSE)])
  names(q_residuals) <- taxa_order
  # Calculate mean value
  mean_q_residual <- mean(q_residuals)
  
  ## Return results
  output_list <- list("mean_q_residual" = mean_q_residual,
                      "all_q_residuals" = q_residuals)
  
  return(output_list)
}



## TIGER (Cummins and McInerney 2011)
TIGER <- function(alignment_path, fast_TIGER_path, sequence_format = "DNA"){
  # Tree Independent Genertion of Evolutionary Rates
  # Function to calculate TIGER values from a multiple sequence alignment using the 
  #   pbfrandsen/fast_TIGER software (available here: https://github.com/pbfrandsen/fast_TIGER)
  
  ## Remove phylogenetically uninformative sites from the alignment
  # Uninformative sites can bias the TIGER values - see List (2021)
  alignment <- read.FASTA(alignment_path, type = sequence_format)
  n_sites <- length(alignment[[1]])
  # Detect informative sites (the sites to keep)
  is_inds <- pis(as.matrix(alignment), what = "index")
  
  if (length(is_inds) > 0) {
    ## If there is one or more informative sites:
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
  } else if (length(is_inds) == 0){
    ## If there are 0 informative sites:
    # Record and return that information
    mean_TIGER_value <- "0_informative_sites"
  }
  
  # Return the output value
  return(mean_TIGER_value)
}


TIGGER <- function(alignment_path, TIGGER_path, sequence_format = "DNA"){
  # Tree Independent Genertion of Evolutionary Rates
  # Function to calculate TIGER values from a multiple sequence alignment using the 
  #   brettc/tigger software (available here: https://github.com/brettc/tigger)
  # Unimplemented - could not get TIGGER working
}



## Likelihood mapping (Strimmer and von Haeseler 1997)
likelihood.mapping <- function(alignment_path, iqtree2_path, iqtree2_number_threads = 1, substitution_model = "MFP", 
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
  
  # Check for identical taxa
  identical_check <- check.iqtree.log.for.identical.sequences(alignment_path, sequence_format = sequence_format)
  identical_sequences_present <- as.logical(identical_check[["identical_sequences_present"]])
  num_unique_taxa <- as.numeric(identical_check[["number_unique_taxa"]])
  
  ## Need four or more taxa to conduct likelihood mapping
  if (num_unique_taxa >= 4){
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
  } else {
    # Create a vector noting that there are insufficient unique taxa to create a likelihood map
    lm_results <- rep(paste0(num_unique_taxa, "_unique_taxa_no_likelihood_map"), 5)
  }
  
  ## Rename vector of results
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


scfl <- function(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", number_scf_quartets = 100, substitution_model = "MFP", force.redo = FALSE){
  # Function to calculate the site concordance factors with maximum likelihood, given a maximum likelihood tree estimated in IQ-Tree
  
  ## Check that the treefile already exists: if it doesn't, run IQ-Tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, estimate the maximum likelihood tree with IQ-Tree2
    call <- paste0(iqtree2_path," -s ",alignment_path," -m ",substitution_model," -nt ", iqtree2_number_threads, " -redo -safe")
    system(call)
  }
  
  # Call IQ-Tree2.2.2 and calculate the scfl (site concordance factors with likelihood)
  if (force.redo == TRUE){
    redo_call <- "-redo "
  } else {
    redo_call <- ""
  }
  # for sCFl: iqtree -t concat.treefile -s ALN_FILE --scfl 100 --prefix concord -nt 10
  treefile <- paste0(alignment_path, ".treefile")
  call <- paste0(iqtree2_path, " -t ", treefile, " -s ", alignment_path, " -m ", substitution_model, " --scfl ", number_scf_quartets,
                 " -nt ", iqtree2_number_threads, " ", redo_call, "-safe -pre scfl")
  system(call)
  
  ## Retrieve the site concordance factors from the output table
  treefile <- paste0(dirname(alignment_path), "/scfl.cf.tree")
  t <- read.tree(treefile)
  scfl_vec <- as.numeric(t$node.label[which(t$node.label != "")])
  scfl_branch_ids <- c(t$edge[,2])[which(t$edge[,2]>Ntip(t))] 
  scfl_results <- list(mean_scf = round(mean(scfl_vec), digits = 2), 
                       median_scf = round(median(scfl_vec), digits = 2), 
                       all_scfs = scfl_vec, 
                       branch_ids = scfl_branch_ids )
  ## Return the site concordance factor results
  return(scfl_results)
}



## Delta plots (Holland et. al. 2002)
mean.delta.plot.value <- function(alignment_path, sequence_format = "DNA", substitution_model = "JC69", 
                                  base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA){
  # This function takes an alignment, calculates a distance matrix for the alignment, and the applies the
  # `ape` function `delta.plot`. We take the mean delta plot value as the test statistic. 
  
  ## Calculate a distance matrix of pairwise distances from DNA sequences using a model of DNA substitution
  # Default model of DNA substitution is JC ("JC69") - it's used to simulate the sequences for the simulations
  pdm <- calculate.dna.pairwise.distance.matrix(alignment_path, sequence_format, substitution_model, base_frequencies, Q_matrix, number_of_rate_categories)
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
recalculate.scf <- function(alignment_path, iqtree2_path,  num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                            iqtree_substitution_model = "MFP", force.redo = FALSE){
  ## Function to take one simulated alignment, apply all treelikeness metrics and return results in a dataframe
  
  ## Prepare variables and output file names for run
  # Get directory path
  replicate_folder <- paste0(dirname(alignment_path), "/")
  # Get unique id for the alignment
  unique_id <- paste(gsub("_output_alignment", "", unlist(strsplit(basename(alignment_path), "\\."))[1:(length(unlist(strsplit(basename(alignment_path), "\\."))) - 1)]), collapse = ".") 
  # Create output folder for sclf
  df_name <- paste0(replicate_folder, unique_id, "_iqtree2.2.2_scfl.csv")
  
  ## Move to replicate folder
  setwd(replicate_folder)
  
  ## Apply Site concordance factors (Minh et. al. 2020)
  # Apply scfl- site concordance factors with likelihood (iqtree2 v2.2.2)
  scfl <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
              substitution_model = iqtree_substitution_model, force.redo = force.redo)
  
  # Assemble results into a dataframe and save
  results_vec <- c(unique_id, scfl$mean_scf, scfl$median_scf, min(scfl$all_scfs), max(scfl$all_scfs), alignment_path)
  results_df <- as.data.frame(matrix(data = results_vec, nrow = 1, ncol = length(results_vec), byrow = TRUE))
  names_vec <- c("unique_id", "sCF_mean", "sCF_median", "sCF_min", "sCF_max", "input_alignment_path")
  names(results_df) <- names_vec
  write.csv(results_df, file = df_name, row.names = FALSE)
  
  # Return the scfl results
  return(results_df)
}



treelikeness.metrics.simulations <- function(alignment_path, 
                                             iqtree2_path, splitstree_path, 
                                             phylogemetric_path, fast_TIGER_path, 
                                             supply_number_of_taxa = FALSE, number_of_taxa = NA, 
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
      # If the number of taxa isn't supplied as an input variable, determine it by finding the number of taxa from a tree in the folder files for the alignment
      if ((grepl("exp1", unique_id)) | (!identical(agrep("random_trees", aln_folder_files), integer(0)))) {
        # If either the unique id contains "exp1" OR there is a file name containing the phrase "random_trees",
        #    open the first random tree and see how many taxa are present
        random_trees_file <- paste0(replicate_folder, grep("random_trees", aln_folder_files, value = TRUE))
        random_trees <- read.tree(random_trees_file)
        n_tree_tips <- unique(Ntip(random_trees))[[1]]
      } else if ((grepl("exp3", unique_id)) | (!identical(agrep("starting_tree", aln_folder_files), integer(0)))) {
        # If either the unique id contains "exp3" OR there is a file name containing the phrase "starting_tree",
        #    open the starting tree and see how many number of taxa present
        start_tree_file <- paste0(replicate_folder, grep("starting_tree", aln_folder_files, value = TRUE))
        start_tree <- read.tree(start_tree_file)
        n_tree_tips <- Ntip(start_tree)
      } 
    }
    
    # Apply all treelikeness tests:
    # Apply Likelihood mapping (Strimmer and von Haeseler 1997)
    lm <- likelihood.mapping(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, substitution_model = iqtree_substitution_model, 
                             number_of_taxa = n_tree_tips, sequence_format = sequence_format)

    # Apply Site concordance factors with likelihood (Minh et. al. 2020): --scfl (iqtree2 v2.2.2)
    scfl <- scfl(alignment_path, iqtree2_path, iqtree2_number_threads = num_iqtree2_threads, number_scf_quartets = num_iqtree2_scf_quartets, 
                 substitution_model = iqtree_substitution_model)

    # Apply Network Treelikeness Test (Huson and Bryant 2006)
    ntlt <- network.treelikeness.test(alignment_path, splitstree_path, sequence_format = sequence_format)

    # Apply Delta plots (Holland et. al. 2002)
    mean_delta_plot_value <- mean.delta.plot.value(alignment_path, sequence_format = sequence_format, substitution_model = distance_matrix_substitution_method,
                                                   base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA)

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



#### Utility functions ####
make.splitstree.neighbornet <- function(alignment_path, splitstree_path, return.splits = TRUE){
  ## Construct a NeighborNet network using SplitsTree
  # Convert fasta to nexus (if the conversion has not already occured)
  check_nexus_path <- paste0(alignment_path,"_converted.nex")
  if (file.exists(check_nexus_path) == FALSE){
    nexus_alignment_path <- convert.to.nexus(alignment_path, sequence_format = "DNA", include_taxablock = TRUE)
  } else if (file.exists(check_nexus_path) == TRUE) {
    nexus_alignment_path <- check_nexus_path
  }
  # Name output path
  splits_output_path <- paste0(alignment_path, "_Splitstree_NeighborNet_splits.nex")
  # Run Splitstree4 if the confidence_path and output_path files do not exist
  if (file.exists(splits_output_path) == FALSE){
    # Assemble the SplitsTree 4 command
    splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", nexus_alignment_path,
                                 "; ASSUME chartransform =Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true;", 
                                 "ASSUME disttransform=NeighborNet; SAVE FILE=",
                                 splits_output_path," REPLACE=yes; QUIT'")
    # Call SplitsTree 4
    system(splitstree_command)
  }
  
  # Open the splits from SplitsTree
  splits <- read.nexus.splits(splits_output_path)
  
  # Return either the splits (when return.splits == TRUE) or the file name of the splits_output_path (when return.splits == FALSE)
  if (return.splits == TRUE){
    return(splits)
  } else if (return.splits == FALSE){
    return(splits_output_path)
  }
}


convert.to.nexus <- function(alignment_path, sequence_format = "DNA", include_taxablock = FALSE){
  ### Convert fasta file to nexus file (if there is no existing nexus file with the same name)
  
  ## Prepare parameters for file conversion
  # Name nexus file by simply appending "_converted.nex" to end of existing file name
  nexus_alignment_path <- paste0(alignment_path,"_converted.nex")
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
      if (include_taxablock == TRUE){
        # write the output as a nexus file with a taxa block (single data block = FALSE)
        write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = FALSE, interleaved = FALSE)
      } else if (include_taxablock == FALSE){
        # write the output as a nexus file without a taxa block - only a single datablock (single data block = TRUE)
        write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = TRUE, interleaved = FALSE)
      }
    }
  } else if (suffix == "phy" | suffix == "phylip"){
    ## If the file is a phy file, convert it to nexus file format (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE){
      data <- read.phy(alignment_path)
      # Write out the nexus data
      if (include_taxablock == TRUE){
        # write the output as a nexus file with a taxa block (single data block = FALSE)
        write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = FALSE, interleaved = FALSE)
      } else if (include_taxablock == FALSE){
        # write the output as a nexus file without a taxa block - only a single datablock (single data block = TRUE)
        write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, datablock = TRUE, interleaved = FALSE)
      }
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


check.iqtree.log.for.identical.sequences <- function(alignment_path, sequence_format = "DNA"){
  ## Function to check whether every sequence in an alignment is unique, using the IQ-Tree log file
  
  # Open IQ-Tree log file
  iqtree_log_file <- paste0(alignment_path, ".log")
  log_lines <- readLines(iqtree_log_file) 
  
  # Find number of sequences in original alignment
  num_taxa_line <- strsplit(log_lines[grepl("Alignment has", log_lines)], " ")[[1]]
  num_total_taxa <- num_taxa_line[3]
  
  # Check log file to see if identical sequences are present
  identical_check_1 <- grep("identical sequences \\(see below\\) will be ignored for subsequent analysis", log_lines)
  identical_check_2.1 <- grep("is identical to", log_lines)
  identical_check_2.2 <- grep("but kept for subsequent analysis", log_lines)
  identical_check_3 <- grep("WARNING: Your alignment contains too many identical sequences!", log_lines)
  identical_check_4.1 <- grep("\\(identical to", log_lines)
  identical_check_4.2 <- grep("\\) is ignored but added at the end", log_lines)
  
  # If one or more line indicating identical sequences is present, then identical sequences are present
  if ((length(identical_check_1) > 0) |
      (length(identical_check_2.1) > 0 & length(identical_check_2.2) > 0) |
      (length(identical_check_3) > 0) |
      (length(identical_check_4.1) > 0 & length(identical_check_4.2) > 0)){
    identical_sequences_present = TRUE
  } else {
    identical_sequences_present = FALSE
  }
  
  # Determine the number of unique sequences 
  # If no identical sequences are present, the number of unique sequences will be equal to the number of taxa
  if (identical_sequences_present == TRUE){
    # Find the alignment with only unique sequences
    unique_seq_check_1 <- grep("For your convenience alignment with unique sequences printed to", log_lines)[1]
    if (length(unique_seq_check_1) > 0){
      unique_seq_path_1 <- gsub(" ", "", gsub("For your convenience alignment with unique sequences printed to","",log_lines[unique_seq_check_1]))
    } else {
      unique_seq_path_1 <- NA
    }
    unique_seq_check_2 <- grep("Alignment was printed to ", log_lines)[1]
    if (length(unique_seq_check_2) > 0){
      unique_seq_path_2 <- gsub(" ", "", gsub("Alignment was printed to ","",log_lines[unique_seq_check_2]))
    } else {
      unique_seq_path_2 <- NA
    }
    
    # If one or more paths exist, open the phylip file containing the identical sequences
    if (length(unique_seq_check_1) > 0 & length(unique_seq_check_2) > 0 &
        is.na(unique_seq_check_1) == FALSE & is.na(unique_seq_check_2) == FALSE){
      if (identical(unique_seq_path_1, unique_seq_path_2) == TRUE){
        # If both paths are identical, doesn't matter which one is selected
        unique_seq_path <- unique_seq_path_1
      } else {
        # If paths are not identical, select the one that explicitly was stated to contain unique sequences
        unique_seq_path <- unique_seq_path_1
      }
      
      # Open the unique sequences file and determine how many unique taxa are present
      unique_seq_dna <- read.phyDat(unique_seq_path, format = "phylip", type = sequence_format)
      number_unique_taxa <- length(names(unique_seq_dna))
    } else if (is.na(unique_seq_check_1) == TRUE | is.na(unique_seq_check_2) == TRUE){
      # If there is no unique sequences file, identify the number of unique sequences by opening and looking at the original alignment
      number_unique_taxa <- length(unique(as.list(read.dna(alignment_path, "fasta"))))
    } else {
      # Backup condition for if something goes wrong
      number_unique_taxa = NA
    }
  } else if (identical_sequences_present == FALSE) {
    # If no paths to uniqueseq files exist, set number of identical sequences to NA
    number_unique_taxa <- num_total_taxa
  }
  
  # Prepare output
  output_vector <- c(identical_sequences_present, num_total_taxa, number_unique_taxa)
  names(output_vector) <- c("identical_sequences_present", "number_total_taxa", "number_unique_taxa")
  
  return(output_vector)
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


trivial.splits.present <- function(s){
  ## Small  function to check whether trivial splits are present and return the number of trivial splits and the number of non-trivial splits
  nSplits <- length(s)
  nTips <- length(attr(s, "labels"))
  l <- lengths(s)
  split_number_trivial <- as.logical((l == 0L) | (l == 1L) | (l == nTips) | (l == (nTips - 1L)))
  trivial_splits_present <- TRUE %in% split_number_trivial
  num_trivial_splits <- length(which((l == 0L) | (l == 1L) | (l == nTips) | (l == (nTips - 1L))))
  num_non_trivial_splits <- nSplits - num_trivial_splits
  # Return output
  op <- list("TrivialSplitsPresent" = trivial_splits_present, "Num_splits" = nSplits, 
             "Num_trivial_splits" = num_trivial_splits, "Num_non_trivial_splits" = num_non_trivial_splits)
  return(op)
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


calculate.dna.pairwise.distance.matrix <- function(alignment_path, sequence_format = "DNA", substitution_model = "JC69", base_frequencies = NA, Q_matrix = NA, number_of_rate_categories = NA){
  ## Calculate a distance matrix of pairwise distances from DNA sequences using a model of DNA substitution
  
  # Identify file type of alignment
  suffix <- tolower(tail(strsplit(alignment_path,"\\.")[[1]],1))
  # Read alignment file
  if (sequence_format == "DNA"){
    # Open DNA alignment
    if (suffix == "fa" | suffix == "fasta" | suffix == "fas" | suffix == "fna" | suffix == "faa" | suffix == "frn"){
      # Open alignment
      alignment <- read.FASTA(alignment_path, type = sequence_format)
    } else if (suffix == "nex" | suffix == "nexus") {
      alignment <- as.DNAbin(read.nexus.data(alignment_path))
    } else if (suffix == "phy"){
      alignment <- as.DNAbin(read.phy(alignment_path))
    }
  }
  
  # Default model of DNA substitution is JC ("JC69") - it's used to simulate the sequences for the simulations
  # Default is to use model of substitution only
  # Can also use number of rate categories, Q matrix, and base frequencies (or any combination of the above) by specifying in the function call
  if ((is.na(number_of_rate_categories) == TRUE) & ((NA %in% Q_matrix) == TRUE) & ((NA %in% base_frequencies) == TRUE)){
    # If none of the number_of_rate_categories, Q_matrix or base frequencies are provided, simply estimate distance matrix using model of substitution
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, " only"))
    pdm <- dist.ml(alignment, model = substitution_model)
  } else if ((is.na(number_of_rate_categories) == FALSE) & ((NA %in% Q_matrix) == TRUE) & ((NA %in% base_frequencies) == TRUE)){
    # Use number of rate categories only
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", number of rate categories"))
    pdm <- dist.ml(alignment, model = substitution_model, k = number_of_rate_categories)
  } else if ((is.na(number_of_rate_categories) == TRUE) & ((NA %in% Q_matrix) == FALSE) & ((NA %in% base_frequencies) == TRUE)){
    # Use Q matrix only
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", Q matrix"))
    pdm <- dist.ml(alignment, model = substitution_model, Q = Q_matrix)
  } else if ((is.na(number_of_rate_categories) == TRUE) & ((NA %in% Q_matrix) == TRUE) & ((NA %in% base_frequencies) == FALSE)){
    # Use base frequencies only
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", base frequencies"))
    pdm <- dist.ml(alignment, model = substitution_model, bf = base_frequencies)
  } else if ((is.na(number_of_rate_categories) == FALSE) & ((NA %in% Q_matrix) == FALSE) & ((NA %in% base_frequencies) == TRUE)){
    # Use number of rate categories and Q matrix
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", Q matrix, number of rate categories"))
    pdm <- dist.ml(alignment, model = substitution_model, Q = Q_matrix, k = number_of_rate_categories)
  } else if ((is.na(number_of_rate_categories) == FALSE) & ((NA %in% Q_matrix) == TRUE) & ((NA %in% base_frequencies) == FALSE)){
    # Use number of rate categories and base frequencies
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", base frequencies, number of rate categories"))
    pdm <- dist.ml(alignment, model = substitution_model, bf = base_frequencies, k = number_of_rate_categories)
  } else if ((is.na(number_of_rate_categories) == TRUE) & ((NA %in% Q_matrix) == FALSE) & ((NA %in% base_frequencies) == FALSE)){
    # Use Q matrix and base frequencies
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", base frequencies, Q matrix"))
    pdm <- dist.ml(alignment, model = substitution_model, bf = base_frequencies, Q = Q_matrix)
  } else if ((is.na(number_of_rate_categories) == FALSE) & ((NA %in% Q_matrix) == FALSE) & ((NA %in% base_frequencies) == FALSE)){
    # Use number of rate categories, Q matrix, and base frequencies
    print(paste0("Estimating pairwise distance matrix with model of sequence evolution ", substitution_model, ", base frequencies, Q matrix, and number of rate categories"))
    pdm <- dist.ml(alignment, model = substitution_model, bf = base_frequencies, Q = Q_matrix, k = number_of_rate_categories)
  }
  
  # Return the pairwise distance matrix
  return(pdm)
}


get.model.parameters.from.iqtree.file <- function(alignment_path, sequence_format = "DNA"){
  # Extract parameters for estimating the distance matrix from the .iqtree file 
  if (sequence_format == "DNA"){
    ## Open log file
    iqtree_file <- paste0(alignment_path, ".iqtree")
    iqtree_lines <- readLines(iqtree_file)
    
    ## Construct a folder for storing the files generated by this function
    al_directory <- paste0(dirname(alignment_path), "/")
    if (dir.exists(al_directory) == FALSE){dir.create(al_directory)}
    
    ## Base frequencies
    # Check whether state frequencies are equal
    check_line <- iqtree_lines[grep("State frequencies: ", iqtree_lines)]
    if (grepl("equal frequencies", check_line) == TRUE){
      # If the .iqtree file line reads "State frequencies: (equal frequencies)", then all frequencies are equal
      state_frequencies= c(0.25, 0.25, 0.25, 0.25)
    } else if (grepl("empirical counts from alignment", check_line) == TRUE){
      start_line_id <- grep("pi\\(A\\)", iqtree_lines)
      freqs_lines <- iqtree_lines[start_line_id:(start_line_id+3)]
      state_frequencies <- as.numeric(unlist(lapply(strsplit(freqs_lines, "="), `[[`, 2)))
    } else {
      # If neither equal frequencies or empirical counts from alignment, set state_frequencies to NA
      state_frequencies = NA
    }
    
    ## Q matrix
    # Identify the lower diagonal Q matrix and the number of rate categories (by extracting from the .iqtree file generated by IQ-Tree during tree estimation)
    # Find start and end of list of models 
    start_line_id <- grep("Rate matrix Q:", iqtree_lines) + 2
    end_line_id <- grep("Model of rate heterogeneity:", iqtree_lines) - 2
    # Get lines with Q matrix 
    Q_lines <- iqtree_lines[start_line_id:end_line_id]
    # Save Q matrix as file
    Q_lines_file <- paste0(al_directory, "Q_Matrix_results.txt")
    write(Q_lines, file = Q_lines_file)
    # Read in Q matrix file as a tsv and neaten the matrix
    Q_df <- read.table(Q_lines_file, sep = "")
    names(Q_df) <- c("base", "A", "C", "G", "T")
    rownames(Q_df) <- c("A", "C", "G", "T")
    Q_df <- Q_df[,2:5]
    # Get the vector of lower diagonal Q matrix
    # Q vector should be: CA, GA, GC, TA, TC, TG (where for ij i is the row and j is the column)
    Q_vector <- c(Q_df[2,1], Q_df[3,1], Q_df[3,2], Q_df[4,1], Q_df[4,2], Q_df[4,3])
    
    ## Rate categories
    # Check whether rate categories are uniform 
    start_line_id <- grep("Model of rate heterogeneity:", iqtree_lines)
    if (grepl("Uniform", iqtree_lines[start_line_id])){
      # If the model of rate heterogeneity is uniform, there is one rate category
      num_rate_categories = 1
    } else {
      # Find the start and end of the number of rate categories
      start_line_id <- grep("Category  Relative_rate  Proportion", iqtree_lines)
      end_line_id <- grep("LIKELIHOOD MAPPING ANALYSIS", iqtree_lines) - 2
      if (grepl("Relative rates are computed", iqtree_lines[end_line_id])){
        # If the end_line_id line is "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category.",
        #   then subtract one from the end_line_id
        # This ensures only the lines for the gamma categories are selected
        end_line_id = end_line_id - 1
      }
      if (length(start_line_id) > 0){
        # If there are rate categories, get them and use that information for the dist.ml function
        # Get lines with rate categories
        rc_lines <- iqtree_lines[start_line_id:end_line_id]
        # Save Q matrix as file
        rc_lines_file <- paste0(al_directory, "rate_categories_results.txt")
        write(rc_lines, file = rc_lines_file)
        # Read rate categories matrix in as tsv
        rc_df <- read.table(rc_lines_file, sep = "", skip = 1)
        names(rc_df) <- c("Category", "Relative rate", "Proportion")
        # Get the number of rate categories
        num_rate_categories <- nrow(rc_df)
      } else {
        # If there are no rate categories, do not make rate category dataframe and return NA
        num_rate_categories = NA
      } 
    }
    
    ## Best model according to BIC
    # Find the start and end of the number of rate categories
    start_line_id <- grep("Best-fit model according to BIC:", iqtree_lines)
    # Get line of text containing best model
    model_line <- iqtree_lines[start_line_id]
    # Extract model from line of text
    best_model <- gsub(" ", "", unlist(strsplit(model_line, ":"))[[2]])
    
  } # end if sequence_format == "DNA"
  
  # Assemble results into a vector
  output_list <- list(state_frequencies = state_frequencies, Q_vector = Q_vector, num_rate_categories = num_rate_categories, best_fit_model = best_model)
  # Return output
  return(output_list)
  
} # end function

