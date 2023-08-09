# caitlinch/treelikeness-metrics/code/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2023

# This program will apply various tests for treelikeness to simulated alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).



#### 1. Set parameters ####
## Directories
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory       <- Directory where collated results from the treelikeness test statistics will be saved
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)
# fast_TIGER_path         <- Path to fast TIGER executable
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

## Control variables
# run_exp1             <- Whether to apply the treelikeness test statistics to the first set of alignments (logical)
# run_exp3             <- Whether to apply the treelikeness test statistics to the second set of alignments (logical)

run_location = "soma"
if (run_location == "local"){
  # Directories
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  results_directory <- paste0(local_directory, "01_results/")
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2_path <- "iqtree2"
  splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
  fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  results_directory <- local_directory
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  
  # Run parameters
  num_cores <- 30
}

# Control variables
run_exp1 <- FALSE
run_exp3 <- TRUE



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_3/"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the test statistics to the alignment within each directory
if (run_exp1 == FALSE){
  ## For experiment 1:
  # Extract all file names from results folder
  results_files <- list.files(results_directory)
  # Open output df and get names of alignments
  exp1_op_file <- paste0(results_directory, grep("rerun", grep("exp1", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  exp1_op_df <- read.csv(exp1_op_file, stringsAsFactors = FALSE)
  # Exp1 encountering errors in all cores. Not running properly. Remove all alignments with substitution rate 1e-04 and 0.001 (too many identical sequences)
  exp1_op_df <- exp1_op_df[(exp1_op_df$tree_depth != 1e-04 & exp1_op_df$tree_depth != 1e-03),]
  # Get list of alignments
  exp1_als <- exp1_op_df$output_alignment_file
  # Apply treelikeness metrics to all alignments 
  mclapply(exp1_als, treelikeness.metrics.simulations,
           iqtree2_path, splitstree_path, 
           phylogemetric_path, fast_TIGER_path, 
           supply_number_of_taxa = FALSE, number_of_taxa = NA, 
           num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
           iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
           num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE, 
           run_splitstree_for_tree_proportion = TRUE, sequence_format = "DNA", 
           apply.TIGER = FALSE, redo = FALSE, 
           save_timers = TRUE,
           mc.cores = num_cores)
  
  # Collect and collate results
  exp1_list <- mclapply(exp1_als, collate.treelikeness.results, experiment_number = 1, mc.cores = num_cores)
  # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
  keep_indexes <- which(!sapply(exp1_list, is.null))
  exp1_list_filtered <- exp1_list[keep_indexes]
  # Save output dataframe
  exp1_df <- as.data.frame(do.call("rbind", exp1_list_filtered))
  exp1_df_name <- paste0(results_directory, "exp1_treelikeness_metrics_collated_results.csv")
  write.csv(exp1_df, exp1_df_name, row.names = FALSE)
}

if (run_exp3 == TRUE){
  ## For experiment 2:
  # Extract all file names from results folder
  results_files <- list.files(results_directory)
  # Open output df and get names of alignments
  exp3_op_file <- paste0(results_directory, grep("rerun", grep("exp3", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  exp3_op_df <- read.csv(exp3_op_file, stringsAsFactors = FALSE)
  # Get list of alignments
  exp3_als <- exp3_op_df$output_alignment_file
  # Apply treelikeness metrics to all alignments 
  mclapply(exp3_als, treelikeness.metrics.simulations,
           iqtree2_path, splitstree_path, 
           phylogemetric_path, fast_TIGER_path, 
           supply_number_of_taxa = FALSE, number_of_taxa = NA, 
           num_iqtree2_threads = 1, num_iqtree2_scf_quartets = 100, 
           iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", 
           num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE, 
           run_splitstree_for_tree_proportion = TRUE, sequence_format = "DNA", 
           apply.TIGER = TRUE, redo = FALSE,
           mc.cores = num_cores)
  
  # Collect and collate results
  exp3_list <- mclapply(exp3_als, collate.treelikeness.results, experiment_number = 3, mc.cores = num_cores)
  # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
  keep_indexes <- which(!sapply(exp3_list, is.null))
  exp3_list_filtered <- exp3_list[keep_indexes]
  # Save output dataframe
  exp3_df <- as.data.frame(do.call("rbind", exp3_list_filtered))
  exp3_df_name <- paste0(results_directory, "exp3_treelikeness_metrics_collated_results.csv")
  write.csv(exp3_df, exp3_df_name, row.names = FALSE)
}


