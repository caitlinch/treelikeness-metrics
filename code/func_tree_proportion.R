# /caitlinch/treelikeness_metrics/func_tree_proportion.R
# Caitlin Cherryh 2022

# This file contains functions to calculate the tree proportion for any alignment
# Some functions require IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, or SplitsTree (4.17.2 or above).

library(phangorn)
library(ape)

#### Tree proportion function ####
tree.proportion <- function(alignment_path, sequence_format = "DNA", remove_trivial_splits = TRUE, network.method = "SplitsTree4" , splitstree_path = NA, dist.ml.model = NA){
  ### Function to calculate the tree proportion: the proportion of split weights in the phylogenetic network captured by the minimum evolution tree
  # Network method can be "SplitsTree4" (to estimate the network in SplitsTree4) or "R" (to estimate the network in R using ape)
  # If using "SplitsTree4" as the network method, need to supply a file path to the SplitsTree4 java executeable
  # If using "R" as the network method, need to supply a model for estimating the distance matrix using ape::dist.ml. All models from ape::dist.ml are allowed.
  
  ### Check whether the alignment file exists
  # If the alignment file does not exist, raise an error and stop the function
  if (is.na(alignment_path)){
    stop('Must supply valid file path to a sequence alignment.' )
  } else if (file.exists(alignment_path) == FALSE){
    stop('Must supply valid file path to a sequence alignment.' )
  }
  
  ### Check whether the sequence format is valid
  if (sequence_format != "DNA" & sequence_format != "AA"){
    stop('Allowable sequence formats are "DNA" and "AA".' )
  }
  
  ### Convert the network method to all uppercase
  network.method <- toupper(network.method)
  
  ### Run the tree proportion metric
  if (network.method == "R" | network.method == "PHANGORN"){
    if (is.na(dist.ml.model) == TRUE){
      stop('Must supply model for calculating pairwise distance matrix. Look at help file for phangorn::dist.ml function (type "?dist.ml" into the command line) to see available substitution models.' )
    } else if (is.na(dist.ml.model) == FALSE){
      # Calculate NeighborNet network in R
      # Estimate a NeighborNet network from the distance matrix and order splits from strongest to weakest
      # Compute pairwise distances for the taxa using the specified model of sequence evolution
      mldist <- calculate.dna.pairwise.distance.matrix(alignment_path, sequence_format = sequence_format, substitution_model = dist.ml.model)
      # Create a NeighbourNet network from the alignment
      nnet <- neighborNet(mldist)
      # Extract the splits from the NeighborNet network
      unordered_nw_splits <- as.splits(nnet)
    }
  } else if (network.method == "SPLITSTREE" | network.method == "SPLITTREE" | network.method == "SPLITSTREE4" | network.method == "S"){
    if (is.na(splitstree_path) == TRUE){
      stop('Must supply file path to SplitsTree4 executable.' )
    } else if (file.exists(splitstree_path) == FALSE){
      stop('Must supply valid file path to SplitsTree4 executable.' )
    } else if (is.na(splitstree_path) == FALSE & file.exists(splitstree_path) == TRUE){
      unordered_nw_splits <- make.splitstree.neighbornet(alignment_path, splitstree_path, return.splits = TRUE)
    }
  } else {
    stop('Valid options for network.method are "SplitsTree4" (to calculate the Neighbor-Net in SplitsTree4) and "R" (to calculate the Neighbor-Net in R using phangorn).' )
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
  
  ## Return result
  # Return the tree proportion value
  return(tree_proportion)
}



#### Distance matrix and network estimation functions ####
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
