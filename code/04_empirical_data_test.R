# caitlinch/treelikeness-metrics/code/04_empirical_data_test.R
# Caitlin Cherryh 2023

## Script summary:
# This program will apply various tests for treelikeness to empirical alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

#### 1. Set parameters ####
## Directories
# output_directory                <- Directory where alignments will be saved/treelikeness metrics will be run
# replicate_alignment_directory   <- Location of bootstrap replicate alignments, generated using Alisim in IQ-Tree
# repo_directory                  <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree (requires version 2.2-beta or above)
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "local"
if (run_location == "local"){
  # Directories
  output_directory                <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/05_empirical_treelikeness_results/"
  replicate_alignment_directory   <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/04_bs_replicate_alignments/"
  repo_directory                  <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2_path      <- "iqtree2"
  splitstree_path   <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  output_directory                <- "/data/caitlin/treelikeness_metrics/"
  replicate_alignment_directory   <- ""
  repo_directory                  <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path      <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path   <- "/home/caitlin/splitstree4/SplitsTree"
  
  # Run parameters
  num_cores <- 30
}

alignment_ids <- c("WEA17", "WEA17F")
control_parameters <- list(edit.replicate.alignments = FALSE,
                           apply.treelikeness.tests = FALSE)


#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_empirical.R"))
source(paste0(repo_directory, "code/func_parametric_bootstraps.R"))
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Calculate the number of parallel processes to run using mc.apply
mclapply_num_cores <- num_cores/10
if (mclapply_num_cores < 1){
  # Minimum of one process running at once
  mclapply_num_cores <- 1
} else {
  # Take floor of the decimal to run a conservative number of processes 
  #     (e.g. if 31 cores and each IQ-Tree run uses 10, you can run 3 alignments at once)
  mclapply_num_cores <- floor(mclapply_num_cores)
}


#### 3. Add gaps to the alignments ####




#### 4. Apply tests for treelikeness to each alignment ####
if (control_parameters$apply.treelikeness.tests == TRUE){
  ## Extract list of alignments
  all_files <- list.files(replicate_alignment_directory, recursive = T)
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17_dir <- paste0(output_directory, "WEA17/")
  if (dir.exists(wea17_dir) == FALSE){ dir.create(wea17_dir) }
  wea17_copy <- paste0(wea17_dir, basename(wea17_al))
  file.copy(from = wea17_al, to = wea17_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  # Apply the wrapper function to calculate treelikeness of bootstrap replicates
  wea17_df <- bootstrap.wrapper(wea17_replicates, 
                                output_directory = output_directory, 
                                splitstree_path = splitstree_path, 
                                iqtree2_path = iqtree2_path, 
                                iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                num_iqtree2_threads = "10", 
                                sequence_format = "AA", 
                                redo = FALSE, 
                                number_parallel_cores = mclapply_num_cores)
  
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17f_dir <- paste0(output_directory, "WEA17F/")
  if (dir.exists(wea17f_dir) == FALSE){ dir.create(wea17f_dir) }
  wea17f_copy <- paste0(wea17f_dir, basename(wea17f_al))
  file.copy(from = wea17f_al, to = wea17f_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17f_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  
  # Apply the wrapper function
  wea17f_df <- bootstrap.wrapper(wea17f_replicates, 
                                 output_directory = output_directory, 
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = "10", 
                                 sequence_format = "AA", 
                                 redo = FALSE, 
                                 number_parallel_cores = mclapply_num_cores)
}



#### 5. Calculate p-values for each test statistic ####





