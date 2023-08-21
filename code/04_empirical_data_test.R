# caitlinch/treelikeness-metrics/code/04_empirical_data_test.R
# Caitlin Cherryh 2023

## Script summary:
# This program will apply various tests for treelikeness to empirical alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

#### 1. Set parameters ####
## Directories
# run_directory           <- Directory where alignments will be saved/treelikeness metrics will be run
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "soma"
if (run_location == "local"){
  # Directories
  output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/04_empirical_treelikeness_results/"
  replicate_alignment_folder <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical_generate_replicate_alignments/"
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2_path <- "iqtree2"
  splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  run_directory <- "/data/caitlin/treelikeness_metrics/"
  replicate_alignment_folder <- ""
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  
  # Run parameters
  num_cores <- 30
}

alignment_ids <- c("WEA17", "WEA17_F")
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



#### 5. Apply tests for treelikeness to each alignment ####
## For WEA17
# Extract the list of alignments using the id (WEA17)

# Apply the wrapper function
wea17_df <- bootstrap.wrapper(bs_rep_al_paths, output_directory = output_directory, 
                                                splitstree_path = splitstree_path, iqtree2_path = iqtree2_path, 
                                                num_iqtree2_threads = "AUTO", sequence_format = "AA", 
                                                redo = FALSE, number_parallel_cores = num_cores)

## For WEA17_filtered
# Extract the list of alignments using the id (WEA17F)

# Apply the wrapper function
wea17_filtered_df <- bootstrap.wrapper(bs_rep_al_paths, output_directory = output_directory, 
                                                         splitstree_path = splitstree_path, iqtree2_path = iqtree2_path, 
                                                         num_iqtree2_threads = "AUTO", sequence_format = "AA", 
                                                         redo = FALSE, number_parallel_cores = num_cores)



#### 4. Calculate p-values for each test statistic ####





