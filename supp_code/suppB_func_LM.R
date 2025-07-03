# caitlinch/treelikeness_metrics/supp_code/suppB_func_LM.R
# Caitlin Cherryh 2025

# This file contains functions to apply Likelihood Mapping to a single alignment
# Some functions require IQ-Tree2 (2.2-beta or above)

## Open packages
library(treestats) # for calculating treeness of quartets
library(Quartet) # for extracting quartets from phylogenetic tree
library(parallel) # for multithreading quartet treeness



## Wrapper function - calls LM in IQ-Tree, returns results
likelihood.mapping.wrapper <- function(
    row_number,
    alignment_dataframe,
    iqtree2_path,
    iqtree2_num_threads = 1
){
  ## Wrapper to take dataframe and feed into LM function
  # Extract row
  row <- alignment_dataframe[row_number, ]
  rownames(row) <- NULL
  # Create results file for the row
  row_results_file <- gsub(
    "output_alignment.fa",
    "LM_results.csv",
    row$output_alignment_file
  )
  quartet_treeness_file <- gsub(
    "output_alignment.fa",
    "quartet_treeness.csv",
    row$output_alignment_file
  )
  # Check whether LM results are present and complete
  run_LM <- check.LM.run.present(
    row_results_file = row_results_file,
    quartet_treeness_file = quartet_treeness_file
  )
  # Apply functions if the results file doesn't exist
  if (run_LM == TRUE){
    # Apply LM function
    row_LM <- likelihood.mapping(
      alignment_path = row$output_alignment_file,
      iqtree2_path = iqtree2_path,
      iqtree2_number_threads = iqtree2_num_threads,
      substitution_model = row$alisim_gene_models,
      number_of_taxa = row$num_taxa,
      sequence_format = row$sequence_type
    )
    # Calculate quartet treeness
    row_quartet_treeness <- quartet.treeness.wrapper(
      iqtree_tree_path = paste0(row$output_alignment_file, ".treefile"),
      num_sampled_quartets = row_LM[["total_num_quartets"]],
      num_threads = iqtree2_num_threads
    )
    # Save the dataframe of quartet treeness to output
    output_quartet_treeness_df <- cbind(row, row_quartet_treeness)
    write.csv(
      output_quartet_treeness_df,
      file = quartet_treeness_file,
      row.names = FALSE
    )
    # Create summary statistics from treeness df
    treeness_summary_stats <- summary(as.numeric(row_quartet_treeness$treeness))
    row_treeness_stats <- c(
      "quartet_treeness_min" = treeness_summary_stats[["Min."]],
      "quartet_treeness_1st_qu" = treeness_summary_stats[["1st Qu."]],
      "quartet_treeness_median" = treeness_summary_stats[["Median"]],
      "quartet_treeness_mean" = treeness_summary_stats[["Mean"]],
      "quartet_treeness_3rd_qu" = treeness_summary_stats[["3rd Qu."]],
      "quartet_treeness_max" = treeness_summary_stats[["Max."]],
      "quartet_treeness_num_sampled_quartets" = nrow(row_quartet_treeness),
      "quartet_treeness_total_num_quartets" = get.total.num.quartets(
        iqtree_tree_path = paste0(row$output_alignment_file, ".treefile")
        )
    )
    # Create output row
    output_row <- c(
      unlist(row),
      row_LM,
      row_treeness_stats
    )
    names(output_row) <- c(
      names(row),
      names(row_LM),
      names(row_treeness_stats)
    )
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


## Determine whether to run the LM (IQ-Tree)
check.LM.run.present <- function(row_results_file, quartet_treeness_file){
  if ((file.exists(row_results_file) == FALSE) |
      (file.exists(quartet_treeness_file) == FALSE)) {
    # One or both of the output files don't exist - run LM
    run_LM <- TRUE
  } else {
    # Check file size
    check_file_size <- file.info(row_results_file)[["size"]]
    if (check_file_size == 0) {
      # Empty file - run LM
      run_LM <- TRUE
    } else {
      # Open file and check all necessary rows are present
      check_row <- read.csv(row_results_file)
      # Check for results in column names
      col_names_check <- c(
        "uid",
        "num_resolved_quartets",
        "num_partly_resolved_quartets",
        "num_unresolved_quartets",
        "total_num_quartets",
        "proportion_resolved_quartets",
        "quartet_treeness_min",
        "quartet_treeness_1st_qu",
        "quartet_treeness_median",
        "quartet_treeness_mean",
        "quartet_treeness_3rd_qu",
        "quartet_treeness_max",
        "quartet_treeness_num_sampled_quartets",
        "quartet_treeness_total_num_quartets"
      )
      missing_cols <- length(which(!col_names_check %in% names(check_row)))
      # Check whether any of those columns are missing (i.e., missing_cols > 0)
      if (missing_cols == 0) {
        # All required rows are present - do not run LM
        run_LM <- FALSE
      } else {
        # One or more rows are missing - run LM
        run_LM <- TRUE
      }
    }
  }
  return(run_LM)
}



## Calculate the total number of quartets in a tree
get.total.num.quartets <- function(iqtree_tree_path){
  # Open tree
  tree <- read.tree(iqtree_tree_path)
  # Identify all quartets in tree (each row = 1 quartet)
  all_quartets <- t(Quartet::AllQuartets(Ntip(tree)))
  # Get total number of quartets in tree
  num_quartets <- nrow(all_quartets)
  return(num_quartets)
}



## Calculate treeness for quartets
quartet.treeness.wrapper <- function(
    iqtree_tree_path,
    num_sampled_quartets,
    num_threads = 1){
  # Open tree
  tree <- read.tree(iqtree_tree_path)
  # Identify all quartets in tree (each row = 1 quartet)
  all_quartets <- t(Quartet::AllQuartets(Ntip(tree)))
  # Reduce to the subset of quartets to sample
  if (num_sampled_quartets < nrow(all_quartets)){
    # Reduce number of quartets to sample
    quartets_to_sample <- sort(sample(
      1:nrow(all_quartets),
      size = num_sampled_quartets
    ))
    sample_quartets <- all_quartets[quartets_to_sample, ]
  } else {
    # Sample all quartets
    sample_quartets <- all_quartets
  }
  # Make sure num_threads is a number (may be "AUTO" for IQ-Tree2)
  mcl_num_threads = tryCatch(
    as.numeric(num_threads),
    warning = function(num_threads){num_threads = 1}
  )
  # Calculate treeness for each quartet
  all_quartet_df <- as.data.frame(do.call(
    rbind,
    mclapply(
      1:nrow(sample_quartets),
      quartet.treeness,
      quartets_df = sample_quartets,
      tree = tree,
      tip_names = tree$tip.label,
      mc.cores = mcl_num_threads
    )
  ))
  names(all_quartet_df) <- c(
    "quartet_num",
    "taxa1",
    "taxa2",
    "taxa3",
    "taxa4",
    "treeness"
  )
  return(all_quartet_df)
}


quartet.treeness <- function(quartet_number, quartets_df, tree, tip_names){
  # Calculate treeness for a single quartet
  quartet_tips <- tip_names[quartets_df[quartet_number, ]]
  quartet_tree <- ape::keep.tip(tree, tip = quartet_tips)
  quartet_treeness <- treestats::treeness(quartet_tree)
  quartet_results <- c(quartet_number, quartet_tips, quartet_treeness)
  return(quartet_results)
}


## Likelihood mapping (Strimmer and von Haeseler 1997)
likelihood.mapping <- function(alignment_path,
                               iqtree2_path,
                               iqtree2_number_threads = 1,
                               substitution_model = "MFP",
                               number_of_taxa = NA,
                               sequence_format = "DNA") {
  ### Function to call IQ-Tree and create a likelihood map for the alignment
  ## Create the likelihood map
  # Check whether likelihood mapping or IQ-Tree have run before.
  # If one or both haven't run IQ-Tree to create the likelihood map
  iq_file <- paste0(alignment_path, ".iqtree")
  map_file <- paste0(alignment_path, ".lmap.eps")
  if ((file.exists(iq_file) == FALSE) |
      (file.exists(map_file) == FALSE)) {
    number_of_quartets <- 25 * as.numeric(number_of_taxa)
    call <- paste0(
      iqtree2_path,
      " -s ",
      alignment_path,
      " -m ",
      substitution_model,
      " -nt ",
      iqtree2_number_threads,
      " -lmap ",
      number_of_quartets,
      " -redo -safe"
    )
    system(call)
  }
  # Check for identical taxa
  identical_check <- check.iqtree.log.for.identical.sequences(
    alignment_path,
    sequence_format = sequence_format
  )
  identical_sequences_present <- as.logical(identical_check[["identical_sequences_present"]])
  num_unique_taxa <- as.numeric(identical_check[["number_unique_taxa"]])
  ## Need four or more taxa to conduct likelihood mapping
  if (num_unique_taxa >= 4) {
    if (file.exists(iq_file) == TRUE) {
      # Extract results from likelihood map
      iq_log <- readLines(iq_file)
      ind <- grep("Number of fully resolved  quartets", iq_log)
      resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind], ":")[[1]][2], "\\(")[[1]][1])
      ind <- grep("Number of partly resolved quartets", iq_log)
      partly_resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind], ":")[[1]][2], "\\(")[[1]][1])
      ind <- grep("Number of unresolved", iq_log)
      unresolved_q <- as.numeric(strsplit(strsplit(iq_log[ind], ":")[[1]][2],"\\(")[[1]][1])
      total_q <- (resolved_q + partly_resolved_q + unresolved_q)
      prop_resolved <- resolved_q / total_q
      # Collate results into a vector
      lm_results <- c(resolved_q,
                      partly_resolved_q,
                      unresolved_q,
                      total_q,
                      prop_resolved)
    }
    else if (file.exists(iq_file) == FALSE) {
      # Create a vector noting that the .iqtree file does not exist
      lm_results <- rep("no_iqtree_file", 5)
    }
  } else {
    # Create a vector noting that there are insufficient unique taxa to create
    #   a likelihood map
    lm_results <- rep(
      paste0(num_unique_taxa, "_unique_taxa_no_likelihood_map"),
      5
    )
  }
  ## Rename vector of results
  names(lm_results) <- c(
    "num_resolved_quartets",
    "num_partly_resolved_quartets",
    "num_unresolved_quartets",
    "total_num_quartets",
    "proportion_resolved_quartets"
  )
  ## Return results
  return(lm_results)
}



check.iqtree.log.for.identical.sequences <- function(
    alignment_path,
    sequence_format = "DNA"
) {
  ## Function to check whether every sequence in an alignment is unique,
  ##    using the IQ-Tree log file
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
      log_lines)[1]
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
        gsub("Alignment was printed to ", "", log_lines[unique_seq_check_2])
      )
    } else {
      unique_seq_path_2 <- NA
    }
    # If one or more paths exist, open the phylip file containing the identical
    #   sequences
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
      unique_seq_dna <- read.phyDat(unique_seq_path, format = "phylip", type = sequence_format)
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
    # If no paths to uniqueseq files exist, number of identical sequences is NA
    number_unique_taxa <- num_total_taxa
  }
  # Prepare output
  output_vector <- c(identical_sequences_present,
                     num_total_taxa,
                     number_unique_taxa)
  names(output_vector) <- c("identical_sequences_present",
                            "number_total_taxa",
                            "number_unique_taxa")
  return(output_vector)
}
