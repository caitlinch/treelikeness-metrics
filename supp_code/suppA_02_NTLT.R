## Supplementary analysis A: Network Treelikeness Test

# caitlinch/treelikeness-metrics/supp_code/suppA_02_NTLT.R
# Caitlin Cherryh 2025

# This program will apply the Network Treelikeness Test to simulated alignments
# This program requires IQ-Tree2 (2.2-beta or above) and SplitsTree (4.17.2 or above).



#### 1. Set parameters ####
## Directories
# local_directory   <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory <- Directory where collated results from the treelikeness test statistics will be saved
# repo_directory    <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path    <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)
# splitstree_path <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores <- Number of parallel threads to use at once

## Control variables
# run_expA1 <- Whether to apply the treelikeness test statistics to the first set of alignments (logical)
# run_expA2 <- Whether to apply the treelikeness test statistics to the second set of alignments (logical)

run_location = "dayhoff"
if (run_location == "WSL"){
  # Directories
  local_directory <- ""
  results_directory <- paste0(local_directory, "suppA/")
  repo_directory <- ""

  # Executable paths
  iqtree2_path <- "iqtree2"
  splitstree_path <- ""

  # Run parameters
  num_cores <- 1
} else if (run_location == "dayhoff"){
  # Directories
  local_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/"
  results_directory <- local_directory
  repo_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/"

  # Executable paths
  iqtree2_path <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/software/iqtree-2.2.2.6-Linux/bin/iqtree2"
  splitstree_path <- "/mnt/data/dayhoff/home/u5348329/splitstree4/SplitsTree"

  # Run parameters
  num_cores <- 20
}

# Control variables
run_expA1 <- FALSE
run_expA2 <- TRUE

# Set number of threads for IQ-Tree2 and number of parallel processes
iqtree2_num_threads <- 4
num_parallel_processes <- ifelse(floor(num_cores / iqtree2_num_threads) == 0,
                                 1,
                                 floor(num_cores / iqtree2_num_threads))


#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "supp_code/suppA_func_NTLT.R"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder
# and apply the test statistics to the alignment within each directory
if (run_expA1 == TRUE){
  ## For experiment A1:
  # Extract all csv files from results dir
  csv_files = grep(".csv",
                   list.files(results_directory, include.dirs = FALSE),
                   value = TRUE)
  # Open output df and get names of alignments
  expA1_op_file <- paste0(results_directory,
                          grep(
                            "file_output_paths",
                            grep("expA1", csv_files, value = TRUE),
                            value = TRUE
                          ))
  expA1_op_df <- read.csv(expA1_op_file, stringsAsFactors = FALSE)
  # Call NTLT wrapper function on each row of the expA1_op_df
  mclapply(
    1:nrow(expA1_op_df),
    network.treelikeness.test.wrapper,
    alignment_dataframe = expA1_op_df,
    splitstree_path = splitstree_path,
    iqtree2_path = iqtree2_path,
    iqtree2_num_threads = iqtree2_num_threads,
    mc.cores = num_parallel_processes
  )
  # Collect all output files
  expA1_ntlt_csvs <- paste0(
    results_directory,
    "exp_A1/",
    grep(
      "NTLT_results.csv",
      grep("expA1_",
           list.files(
             paste0(results_directory, "exp_A1/"),
             recursive = TRUE),
           value = TRUE),
      value = TRUE)
  )
  # Read in all output files and collate
  expA1_ntlt_rows <- lapply(
    expA1_ntlt_csvs,
    read.csv
  )
  expA1_ntlt_df <- as.data.frame(do.call(rbind, expA1_ntlt_rows))
  write.csv(
    expA1_ntlt_df,
    file = paste0(results_directory, "expA1_NTLT_results_collated.csv"),
    row.names = FALSE
  )
}

if (run_expA2 == TRUE){
  ## For experiment A2:
  # Extract all csv files from results dir
  csv_files = grep(".csv",
                   list.files(results_directory, include.dirs = FALSE),
                   value = TRUE)
  # Open output df and get names of alignments
  expA2_op_file <- paste0(results_directory,
                          grep(
                            "file_output_paths",
                            grep("expA2", csv_files, value = TRUE),
                            value = TRUE
                          ))
  expA2_op_df <- read.csv(expA2_op_file, stringsAsFactors = FALSE)
  # Sort by number of trees and tree depth
  expA2_op_df <- expA2_op_df[order(
    expA2_op_df$total_alignment_length,
    expA2_op_df$num_trees
    ), ]
  # Call NTLT wrapper function on each row of the expA2_op_df
  mclapply(
    1:nrow(expA2_op_df),
    network.treelikeness.test.wrapper,
    alignment_dataframe = expA2_op_df,
    splitstree_path = splitstree_path,
    iqtree2_path = iqtree2_path,
    iqtree2_num_threads = iqtree2_num_threads,
    mc.cores = num_parallel_processes
  )
  # Collect all output files
  expA2_ntlt_csvs <- paste0(
    results_directory,
    "exp_A2/",
    grep(
      "NTLT_results.csv",
      grep("expA2_",
           list.files(
             paste0(results_directory, "exp_A2/"),
             recursive = TRUE),
           value = TRUE),
      value = TRUE)
  )
  # Read in all output files and collate
  expA2_ntlt_rows <- lapply(
    expA2_ntlt_csvs,
    read.csv
  )
  expA2_ntlt_df <- as.data.frame(do.call(rbind, expA2_ntlt_rows))
  write.csv(
    expA2_ntlt_df,
    file = paste0(results_directory, "expA2_NTLT_results_collated.csv"),
    row.names = FALSE
  )
}

