# caitlinch/treelikeness-metrics/code/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2023

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
  splitstree_path <-

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
run_expA1 <- TRUE
run_expA2 <- TRUE



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
  if (num_cores == 1){
    lapply(
      1:nrow(expA1_op_df),
      network.treelikeness.test.wrapper,
      alignment_dataframe = expA1_op_df,
      splitstree_path = splitstree_path,
      iqtree2_path = iqtree2_path,
      iqtree2_num_threads = 2
    )
  } else {
    mclapply(
      1:nrow(expA1_op_df),
      network.treelikeness.test.wrapper,
      alignment_dataframe = expA1_op_df,
      splitstree_path = splitstree_path,
      iqtree2_path = iqtree2_path,
      iqtree2_num_threads = 2,
      mc.cores = num_cores/2)
  }

  # Collect all output files
  expA1_ntlt_csvs <- paste0(
    results_directory,
    "exp_A1/",
    grep(
      "NTLT_results.csv",
      grep("expA1_",
           list.files(
             paste0(results_directory, "exp_A1/"),
             recursive = TRUE,
             value = TRUE),
           value = TRUE
      )
    )
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
  ## For experiment 1:
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
  # Call NTLT wrapper function on each row of the expA2_op_df
  if (num_cores == 1){
    lapply(
      1:nrow(expA2_op_df),
      network.treelikeness.test.wrapper,
      alignment_dataframe = expA2_op_df,
      splitstree_path = splitstree_path,
      iqtree2_path = iqtree2_path,
      iqtree2_num_threads = 2
    )
  } else {
    mclapply(
      1:nrow(expA2_op_df),
      network.treelikeness.test.wrapper,
      alignment_dataframe = expA2_op_df,
      splitstree_path = splitstree_path,
      iqtree2_path = iqtree2_path,
      iqtree2_num_threads = 2,
      mc.cores = num_cores/2)
  }

  # Collect all output files
  expA2_ntlt_csvs <- paste0(
    results_directory,
    "exp_A2/",
    grep(
      "NTLT_results.csv",
      grep("expA2_",
           list.files(
             paste0(results_directory, "exp_A2/"),
             recursive = TRUE,
             value = TRUE),
           value = TRUE
      )
    )
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

