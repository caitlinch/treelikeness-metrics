## Supplementary analysis B: Likelihood Mapping (proportion of resolved quartets)

# caitlinch/treelikeness-metrics/supp_code/suppB_02_LM.R
# Caitlin Cherryh 2025

# This program will apply likelihood mapping in IQ-Tree
# This program requires IQ-Tree2 (2.2-beta or above)



#### 1. Set parameters ####
## Directories
# local_directory   <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory <- Directory where collated results from the treelikeness test statistics will be saved
# repo_directory    <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path    <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)

## Run parameters
# num_cores <- Number of parallel threads to use at once

## Control variables
# run_expB1 <- Whether to apply the treelikeness test statistics to the first set of alignments (logical)
# run_expB2 <- Whether to apply the treelikeness test statistics to the second set of alignments (logical)

run_location = "dayhoff"
if (run_location == "WSL"){
  # Directories
  local_directory <- ""
  results_directory <- paste0(local_directory, "suppA/")
  repo_directory <- ""

  # Executable paths
  iqtree2_path <- "iqtree2"

  # Run parameters
  num_cores <- 1
} else if (run_location == "dayhoff"){
  # Directories
  local_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppB/"
  results_directory <- local_directory
  repo_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/"

  # Executable paths
  iqtree2_path <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/software/iqtree-2.2.2.6-Linux/bin/iqtree2"

  # Run parameters
  num_cores <- 20
}

# Control variables
run_expB1 <- TRUE
run_expB2 <- TRUE

# Set number of parallel processes
# 2 threads for IQ-Tree, so number of processes = num_cores/2 (rounded down, must be >=1)
mcl_cores <- ifelse(floor(num_cores/2) == 0, 1, floor(num_cores/2))



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "supp_code/suppB_func_LM.R"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder
# and apply the test statistics to the alignment within each directory
if (run_expB1 == TRUE){
  ## For experiment B1:
  # Extract all csv files from results dir
  csv_files = grep(".csv",
                   list.files(results_directory, include.dirs = FALSE),
                   value = TRUE)
  # Open output df and get names of alignments
  expB1_op_file <- paste0(results_directory,
                          grep(
                            "file_output_paths",
                            grep("expB1", csv_files, value = TRUE),
                            value = TRUE
                          ))
  expB1_op_df <- read.csv(expB1_op_file, stringsAsFactors = FALSE)
  # Call LM wrapper function on each row of the expB1_op_df
  mclapply(
    1:nrow(expB1_op_df),
    likelihood.mapping.wrapper,
    alignment_dataframe = expB1_op_df,
    iqtree2_path = iqtree2_path,
    iqtree2_num_threads = 2,
    mc.cores = mcl_cores
  )
  # Collect all output files
  expB1_lm_csvs <- paste0(
    results_directory,
    "exp_B1/",
    grep(
      "LM_results.csv",
      grep("expB1_",
           list.files(
             paste0(results_directory, "exp_B1/"),
             recursive = TRUE),
           value = TRUE),
      value = TRUE)
  )
  # Read in all output files and collate
  expB1_lm_rows <- lapply(
    expB1_lm_csvs,
    read.csv
  )
  expB1_lm_df <- as.data.frame(do.call(rbind, expB1_lm_rows))
  write.csv(
    expB1_lm_df,
    file = paste0(results_directory, "expB1_LM_results_collated.csv"),
    row.names = FALSE
  )
}

if (run_expB2 == TRUE){
  ## For experiment B2:
  # Extract all csv files from results dir
  csv_files = grep(".csv",
                   list.files(results_directory, include.dirs = FALSE),
                   value = TRUE)
  # Open output df and get names of alignments
  expB2_op_file <- paste0(results_directory,
                          grep(
                            "file_output_paths",
                            grep("expB2", csv_files, value = TRUE),
                            value = TRUE
                          ))
  expB2_op_df <- read.csv(expB2_op_file, stringsAsFactors = FALSE)
  # Call LM wrapper function on each row of the expB2_op_df
  mclapply(
    1:nrow(expB2_op_df),
    likelihood.mapping.wrapper,
    alignment_dataframe = expB2_op_df,
    iqtree2_path = iqtree2_path,
    iqtree2_num_threads = 2,
    mc.cores = mcl_cores
  )
  # Collect all output files
  expB2_lm_csvs <- paste0(
    results_directory,
    "exp_B2/",
    grep(
      "LM_results.csv",
      grep("expB2_",
           list.files(
             paste0(results_directory, "exp_B2/"),
             recursive = TRUE),
           value = TRUE),
      value = TRUE)
  )
  # Read in all output files and collate
  expB2_lm_rows <- lapply(
    expB2_lm_csvs,
    read.csv
  )
  expB2_lm_df <- as.data.frame(do.call(rbind, expB2_lm_rows))
  write.csv(
    expB2_lm_df,
    file = paste0(results_directory, "expB2_LM_results_collated.csv"),
    row.names = FALSE
  )
}


