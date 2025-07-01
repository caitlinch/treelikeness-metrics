# caitlinch/treelikeness_metrics/supp_code/suppA_func_NTLT.R
# Caitlin Cherryh 2023

# This file contains functions to apply the Network Treelikeness Test to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above) or SplitsTree (4.17.2 or above).

## Load required packages
library(ape) # for general tree/alignment wrangling, and the delta.plots function
library(phangorn) # for splits and networks, for midpoint rooting trees


## Wrapper function - calls NTLT, infers tree with IQ-Tree, returns results
network.treelikeness.test.wrapper <- function(
    row_number,
    alignment_dataframe,
    splitstree_path,
    iqtree2_path,
    iqtree2_num_threads = 1
){
  ## Wrapper to take dataframe and feed into NTLT function
  # Extract row
  row <- alignment_dataframe[row_number, ]
  # Create results file for the row
  row_results_file = gsub(
    "output_alignment.fa",
    "NTLT_results.csv",
    row$output_alignment_file
  )
  # Apply functions if the results file doesn't exist
  if (file.exists(row_results_file) == FALSE){
    # Apply NTLT function
    row_ntlt <- network.treelikeness.test(
      alignment_path = row$output_alignment_file,
      splitstree_path = splitstree_path,
      sequence_format = row$sequence_type,
      nexus.file.format = FALSE
    )
    # Infer tree with IQ-Tree
    row_iqtree_files <- infer.tree.IQTree2(
      alignment_path = row$output_alignment_file,
      substitution_model = row$alisim_gene_models,
      iqtree2_path = iqtree2_path,
      iqtree2_num_threads = iqtree2_num_threads
    )
    # Extract number of splits in the tree
    row_tree_splits <- extract.splits.IQTree2(row_iqtree_files[["iqtree_tree_file"]])
    # Check whether this alignment contains any identical sequences
    row_identical_check <- check.iqtree.log.for.identical.sequences(
      row$output_alignment_file,
      sequence_format = row$sequence_type
    )
    # Create output row
    output_row <- c(
      unlist(row),
      as.logical(row_identical_check[["identical_sequences_present"]]),
      row_identical_check[["number_unique_taxa"]],
      row_ntlt,
      row_tree_splits
    )
    names(output_row) <- c(
      names(row),
      "identical_sequences_present",
      "number_unique_taxa",
      names(row_ntlt),
      names(row_tree_splits)
    )
    # Remove row number ("X") column
    output_row <- output_row[2:length(output_row)]
    # Save results as csv (allows you to run in batches and extract files after)
    out_csv <- as.data.frame(rbind(output_row))
    rownames(out_csv) <- NULL
    write.csv(
      out_csv,
      file = row_results_file,
      row.names = FALSE
    )
  }
  # Create status row to output
  status_op <- paste0(
    "Row ",
    row_number,
    " complete: ",
    row_results_file
  )
  # Return results
  return(status_op)
}



## Infer tree with IQ-Tree2
infer.tree.IQTree2 <- function(alignment_path,
                               substitution_model,
                               iqtree2_path,
                               iqtree2_num_threads) {
  ## Infer tree for provided alignment using IQ-Tree
  # Check whether likelihood mapping or IQ-Tree have run before.
  # If one or both haven't run IQ-Tree to create the likelihood map
  iq_file <- paste0(alignment_path, ".iqtree")
  if (file.exists(iq_file) == FALSE) {
    call <- paste0(
      iqtree2_path,
      " -s ",
      alignment_path,
      " -m ",
      substitution_model,
      " -nt ",
      iqtree2_num_threads,
      " -redo -safe"
    )
    system(call)
  }
  # Return vector of IQ-Tree2 output paths
  iq_op <- c(
    "iqtree_file" = paste0(alignment_path, ".iqtree"),
    "iqtree_log_file" = paste0(alignment_path, ".log"),
    "iqtree_tree_file" = paste0(alignment_path, ".treefile"),
    "iqtree_mldist_file" = paste0(alignment_path, ".mldist")
  )
  return(iq_op)
}



## Extract splits in tree inferred by IQ-Tree2
extract.splits.IQTree2 <- function(iqtree_treefile){
  ## Check number of splits in tree
  # Open tree
  tree <- read.tree(iqtree_treefile)
  # Get number of tips in tree
  tree_num_tips <- Ntip(tree)
  # Extract all edges
  tree_branches <- as.data.frame(tree$edge)
  names(tree_branches) <- c("node1", "node2")
  tree_branches$edge_length <- tree$edge.length
  # Extract trivial branches i.e., branches where "node2" <= Ntip(tree)
  trivial_splits <- tree_branches[which(tree_branches$node2 <= tree_num_tips), ]
  # Extract nontrivial branches i.e., branches where "node2" > Ntip(tree)
  nontrivial_splits <- tree_branches[which(tree_branches$node2 > tree_num_tips), ]
  # Create output row
  tree_splits_op <- c(
    "tree_num_splits" = nrow(tree_branches),
    "tree_num_trivial_splits" = nrow(trivial_splits),
    "tree_num_nontrivial_splits" = nrow(nontrivial_splits)
  )
  return(tree_splits_op)
}



## Network tree-likeness test (Huson and Bryant 2006)
network.treelikeness.test <- function(
    alignment_path,
    splitstree_path,
    sequence_format = "DNA",
    nexus.file.format = TRUE
){
  ## Uses Splitstree4.17.1 software to implement the Network Treelikeness Test
  ##    described in Huson and Bryant (2006)
  ## Test steps:
  #   1. Infer a split network N from some data
  #         - Can be done in Splitstree using Networks -> NeighborNet
  #   2. Construct a confidence network for N with level 1 - alpha
  #         - Construct the set of splits with CIs excluding 0
  #         - Can be done in Splitstree4 using `Analysis -> Bootstrap...`
  #           then `Analysis -> Show Confidence Network`
  #   3. If the confidence network does not contain a tree, reject the null
  #      hypothesis that the data originated in a tree
  #         - Reject the null hypothesis if and only if this set is incompatible
  ## Construct and bootstrap a NeighborNet network
  if (nexus.file.format == TRUE){
    nexus_alignment_path = alignment_path
  } else {
    # Convert fasta to nexus
    nexus_alignment_path <- convert.to.nexus(
      alignment_path,
      sequence_format = "DNA",
      include_taxablock = TRUE)
  }
  # Name output path
  confidence_path <- paste0(alignment_path, "_confidence.nexus")
  output_path <- paste0(alignment_path, "_Splitstree_output.nex")
  # Run Splitstree4 if the confidence_path and output_path files do not exist
  if (
    (file.exists(confidence_path) == FALSE) |
    (file.exists(output_path) == FALSE)
  ){
    # Assemble the SplitsTree 4 command
    splitstree_command <- paste0(
      splitstree_path,
      " -g -x 'OPEN FILE=", nexus_alignment_path, ";",
      " ASSUME chartransform=Uncorrected_P",
      " HandleAmbiguousStates=Ignore",
      " Normalize=true;",
      " ASSUME disttransform=NeighborNet;",
      " bootstrap runs=100;",
      " confidence_splits level=95 file=", confidence_path, ";",
      " export file=", output_path,
      " REPLACE=yes;",
      " quit;'"
    )
    # Call SplitsTree 4
    system(splitstree_command)
  }
  ## Construct a confidence network using the bootstrap splits
  # Read in the nexus splits from the confidence network
  splits <- suppressWarnings(read.nexus.splits(confidence_path))
  ## Read in the splits from the confidence network
  # Find the starting line for the splits (line after "MATRIX"),
  #   and read in from that line down as tsv
  splits_start <- grep("MATRIX", readLines(confidence_path))
  splits_df <- read.delim(confidence_path, header = FALSE, skip = splits_start)
  # Remove rows that do not contain splits from data frame
  splits_df <- splits_df[1:length(splits), ]
  # Format data frame columns
  names(splits_df) <- c("split_and_size", "weight", "interval", "taxa")
  splits_df$taxa <- gsub(",", "", splits_df$taxa)
  ss <- gsub(" ", "", gsub("size=", "", gsub("\\]", "", gsub(
    "\\[", "", unlist(strsplit(splits_df$split_and_size, ","))
  ))))
  splits_df$split <- as.numeric(ss[c(TRUE, FALSE)])
  splits_df$size <- as.numeric(ss[c(FALSE, TRUE)])
  si <- gsub(" ", "", gsub("\\)", "", gsub("\\(", "", unlist(
    strsplit(splits_df$interval, ",")
  ))))
  splits_df$interval_start <- as.numeric(si[c(TRUE, FALSE)])
  splits_df$interval_end <- as.numeric(si[c(FALSE, TRUE)])
  # Reorder data frame columns
  splits_df <- splits_df[, c("split",
                             "size",
                             "interval",
                             "interval_start",
                             "interval_end",
                             "taxa")]
  # Add a column for splits with an interval that contain 0
  splits_df$interval_includes_0 <- unlist(
    lapply(
      1:nrow(splits_df),
      check.confidence.interval,
      splits_df = splits_df
    )
  )
  ## Determine Network Treelikeness Test result
  # First, determine which splits have confidence intervals excluding 0
  which_splits <- which(splits_df$interval_includes_0 == FALSE)
  # Second, use that set to determine the Network Treelikeness Test result
  if (length(which_splits) == 0){
    ntlt_result <- "Zero_splits_where_confidence_intervals_exclude_0"
  } else if (length(which_splits) == 1){
    # If there is only one split with confidence intervals excluding 0:
    # The set is compatible by default
    ntlt_result <- "Treelike"
  } else if (length(which_splits) > 1){
    # If there is one or more split with confidence intervals excluding 0,
    #   construct the set of splits with confidence intervals excluding 0
    test_df <- splits_df[which_splits, ]
    test_splits <- splits[which_splits]
    # A set of splits is compatible if all pairwise comparisons between
    #   splits are compatible
    # Check using the ape function is.compatible
    #   (requires splits to be in bitsplits format, which is also in ape)
    compatibility <- is.compatible.bitsplits(as.bitsplits(test_splits))
    # A set of splits is compatible if all pairwise comparisons between splits are compatible
    if (compatibility == FALSE){
      # Some pairwise comparisons between splits are incompatible:
      # the null hypothesis that data was originated in a tree is rejected
      ntlt_result <- "Non-treelike"
    } else if (compatibility == TRUE){
      # All pairwise comparisons between splits are compatible:
      # the null hypothesis that data was originated in a tree is accepted
      ntlt_result <- "Treelike"
    }
  }
  ## Create an output vector for the results
  ## Count number of splits
  output_vector <- c(
    "NetworkTreelikenessTest" = ntlt_result,
    "network_num_splits" = nrow(splits_df),
    "network_num_trivial_splits" = length(which(splits_df$size == 1)),
    "network_num_nontrivial_splits" = length(which(splits_df$size != 1)),
    "confidence_network_num_splits" = length(
      which(splits_df$interval_includes_0 == FALSE)
    ),
    "confidence_network_num_trivial_splits" = length(
      which(splits_df$size == 1 &
              splits_df$interval_includes_0 == FALSE)
    ),
    "confidence_network_num_nontrivial_splits" = length(
      which(splits_df$size != 1 &
              splits_df$interval_includes_0 == FALSE)
    )
  )
  ## Write the splits df to file
  write.csv(
    splits_df,
    file = gsub("output_alignment.fa", "splits_record.csv", alignment_path),
    row.names = FALSE
  )
  ## Return Network Treelikeness Test results
  return(output_vector)
}



check.confidence.interval <- function(split_row_number, splits_df){
  ## Function to take a single split and test if confidence interval includes 0
  splits_row <- splits_df[split_row_number, ]
  interval_start <- splits_row$interval_start
  interval_end <- splits_row$interval_end
  # Check if interval includes 0
  if (interval_start == 0 & interval_end == 0){
    # Both start and end == 0 - interval includes 0
    zero_check <- TRUE
  } else if (interval_start > 0 & interval_end > 0){
    # Both start and end > 0 - interval does not include 0
    zero_check <- FALSE
  } else if (interval_start < 0 & interval_end < 0){
    # Both start and end < 0 - interval does not include 0
    zero_check <- FALSE
  } else if (interval_start <=0 & interval_end >= 0){
    # Start < 0, end > 0 - interval includes 0
    zero_check <- TRUE
  }
  return(zero_check)
}



convert.to.nexus <- function(
    alignment_path,
    sequence_format = "DNA",
    include_taxablock = FALSE
){
  ### Convert fasta file to nexus file
  ###   (if there is no existing nexus file with the same name)
  ## Prepare parameters for file conversion
  # Name nexus file by simply appending "_converted.nex" to end of file name
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
  # Create a variable to specify whether to include a single DATA block
  # (datablock = TRUE) or separate TAXA and CHARACTER boxes (datablock = FALSE)
  if (include_taxablock == TRUE){
    datablock_bool = FALSE
  } else if (include_taxablock == FALSE){
    datablock_bool = TRUE
  }
  ## Convert to nexus using functions based on suffix
  if (suffix == "fasta" |
      suffix == "fa" |
      suffix == "fna" |
      suffix == "ffn" |
      suffix == "faa" |
      suffix == "frn" |
      suffix == "fas") {
    ## If the file is a fasta file, convert it to nexus file format
    ##    (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE) {
      # Read in the fasta data
      data <- read.FASTA(alignment_path, type = sequence_format)
      # Write out the nexus data
      if (include_taxablock == TRUE) {
        # Write the output as a nexus file with a taxa block
        #   (single data block = FALSE)
        write.nexus.data(
          data,
          file = nexus_alignment_path,
          format = nexus_format,
          datablock = FALSE,
          interleaved = FALSE
        )
      } else if (include_taxablock == FALSE) {
        # write the output as a nexus file without a taxa block -
        #   only a single datablock (single data block = TRUE)
        write.nexus.data(
          data,
          file = nexus_alignment_path,
          format = nexus_format,
          datablock = TRUE,
          interleaved = FALSE
        )
      }
    }
  } else if (suffix == "phy" | suffix == "phylip") {
    ## If the file is a phy file, convert it to nexus file format
    ##    (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE) {
      data <- read.phy(alignment_path)
      # Write out the nexus data
      if (include_taxablock == TRUE) {
        # write the output as a nexus file with a taxa block
        #   (single data block = FALSE)
        write.nexus.data(
          data,
          file = nexus_alignment_path,
          format = nexus_format,
          datablock = FALSE,
          interleaved = FALSE
        )
      } else if (include_taxablock == FALSE) {
        # write the output as a nexus file without a taxa block -
        #   only a single datablock (single data block = TRUE)
        write.nexus.data(
          data,
          file = nexus_alignment_path,
          format = nexus_format,
          datablock = TRUE,
          interleaved = FALSE
        )
      }
    }
  }
  ## Open the nexus file and delete the interleave = YES or INTERLEAVE = NO
  ##    part so IQ-TREE can read it
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



check.iqtree.log.for.identical.sequences <- function(alignment_path, sequence_format = "DNA") {
  ## Function to check whether every sequence in an alignment is unique, using the IQ-Tree log file
  # Open IQ-Tree log file
  iqtree_log_file <- paste0(alignment_path, ".log")
  log_lines <- readLines(iqtree_log_file)
  # Find number of sequences in original alignment
  num_taxa_line <- strsplit(log_lines[grepl("Alignment has", log_lines)], " ")[[1]]
  num_total_taxa <- num_taxa_line[3]
  # Check log file to see if identical sequences are present
  identical_check_1 <- grep(
    "identical sequences \\(see below\\) will be ignored for subsequent analysis",
    log_lines
  )
  identical_check_2.1 <- grep("is identical to", log_lines)
  identical_check_2.2 <- grep("but kept for subsequent analysis", log_lines)
  identical_check_3 <- grep(
    "WARNING: Your alignment contains too many identical sequences!",
    log_lines
  )
  identical_check_4.1 <- grep("\\(identical to", log_lines)
  identical_check_4.2 <- grep("\\) is ignored but added at the end", log_lines)
  # If one or more line indicating identical sequences is present, then
  #   identical sequences are present
  if ((length(identical_check_1) > 0) |
      (length(identical_check_2.1) > 0 &
       length(identical_check_2.2) > 0) |
      (length(identical_check_3) > 0) |
      (length(identical_check_4.1) > 0 &
       length(identical_check_4.2) > 0)) {
    identical_sequences_present = TRUE
  } else {
    identical_sequences_present = FALSE
  }
  # Determine the number of unique sequences
  # If no identical sequences are present, the number of unique sequences will
  #   be equal to the number of taxa
  if (identical_sequences_present == TRUE) {
    # Find the alignment with only unique sequences
    unique_seq_check_1 <- grep(
      "For your convenience alignment with unique sequences printed to",
      log_lines
    )[1]
    if (length(unique_seq_check_1) > 0) {
      unique_seq_path_1 <- gsub(
        " ",
        "",
        gsub(
          "For your convenience alignment with unique sequences printed to",
          "",
          log_lines[unique_seq_check_1]
        )
      )
    } else {
      unique_seq_path_1 <- NA
    }
    unique_seq_check_2 <- grep("Alignment was printed to ", log_lines)[1]
    if (length(unique_seq_check_2) > 0) {
      unique_seq_path_2 <- gsub(
        " ",
        "",
        gsub(
          "Alignment was printed to ",
          "",
          log_lines[unique_seq_check_2]
        )
      )
    } else {
      unique_seq_path_2 <- NA
    }
    # If one or more paths exist, open the phylip file containing the
    #     identical sequences
    if (length(unique_seq_check_1) > 0 &
        length(unique_seq_check_2) > 0 &
        is.na(unique_seq_check_1) == FALSE &
        is.na(unique_seq_check_2) == FALSE) {
      if (identical(unique_seq_path_1, unique_seq_path_2) == TRUE) {
        # If both paths are identical, doesn't matter which one is selected
        unique_seq_path <- unique_seq_path_1
      } else {
        # If paths are not identical, select the one that explicitly was stated
        #   to contain unique sequences
        unique_seq_path <- unique_seq_path_1
      }
      # Open the unique sequences file and determine how many unique taxa are present
      unique_seq_dna <- read.phyDat(
        unique_seq_path,
        format = "phylip",
        type = sequence_format
      )
      number_unique_taxa <- length(names(unique_seq_dna))
    } else if (is.na(unique_seq_check_1) == TRUE |
               is.na(unique_seq_check_2) == TRUE) {
      # If there is no unique sequences file, identify the number of unique
      #   sequences by opening and looking at the original alignment
      number_unique_taxa <- length(unique(as.list(read.dna(
        alignment_path, "fasta"
      ))))
    } else {
      # Backup condition for if something goes wrong
      number_unique_taxa = NA
    }
  } else if (identical_sequences_present == FALSE) {
    # If no paths to uniqueseq files exist, set num of identical sequences to NA
    number_unique_taxa <- num_total_taxa
  }
  # Prepare output
  unique_check_vector <- c(
    "identical_sequences_present" = identical_sequences_present,
    "number_total_taxa" = num_total_taxa,
    "number_unique_taxa" = number_unique_taxa
  )
  return(unique_check_vector)
}


