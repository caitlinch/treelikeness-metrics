# caitlinch/treelikeness_metrics/supp_code/suppA_func_NTLT.R
# Caitlin Cherryh 2023

# This file contains functions to apply the Network Treelikeness Test to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above) or SplitsTree (4.17.2 or above).

## Load required packages
library(ape) # for general tree/alignment wrangling, and the delta.plots function
library(phangorn) # for splits and networks, for midpoint rooting trees


## Network tree-likeness test (Huson and Bryant 2006)
network.treelikeness.test <- function(alignment_path, splitstree_path, sequence_format = "DNA", nexus.file.format = TRUE){
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
  if (nexus.file.format == TRUE){
    nexus_alignment_path = alignment_path
  } else {
    # Convert fasta to nexus
    nexus_alignment_path <- convert.to.nexus(alignment_path, sequence_format = "DNA", include_taxablock = TRUE)
  }
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
