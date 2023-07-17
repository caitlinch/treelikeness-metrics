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
# run_exp2             <- Whether to apply the treelikeness test statistics to the second set of alignments (logical)
# rerun_missing_runs   <- Whether or not to rerun alignments that do not have treelikeness test metrics complete (logical)
# rerun_experiment_ids <- Which experiments to check for incomplete runs and rerun missing alignments. 
#                           To check both exp1 and exp2, set <- c("exp1", "exp2")

run_location = "soma"
if (run_location == "local"){
  # Directories
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  results_directory <- paste0(local_directory, "01_results/")
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2.2.2_path <- "/Users/caitlincherryh/Documents/Executables/iqtree-2.2.2-MacOSX/bin/iqtree2"
  
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  results_directory <- local_directory
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2.2.2_path <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  
  # Run parameters
  num_cores <- 30
}



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the updated scfs (scfl) to the alignment within each directory

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
# Apply new scfl to all alignments 
exp1_scfl_list <- mclapply(exp1_als, recalculate.scf,
                           iqtree2_path = iqtree2.2.2_path,  num_iqtree2_threads = "AUTO", 
                           num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "MFP",
                           mc.cores = num_cores)
# Save output dataframe
exp1_df <- as.data.frame(do.call("rbind", exp1_scfl_list))
exp1_df_name <- paste0(results_directory, "exp1_iqtree2.2.2_scfl_results.csv")
write.csv(exp1_df, exp1_df_name, row.names = FALSE)

## For experiment 2:
# Extract all file names from results folder
results_files <- list.files(results_directory)
# Open output df and get names of alignments
exp2_op_file <- paste0(results_directory, grep("rerun", grep("exp2", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
exp2_op_df <- read.csv(exp2_op_file, stringsAsFactors = FALSE)
# Get list of alignments
exp2_als <- exp2_op_df$output_alignment_file
# Apply new scfl to all alignments 
exp2_scfl_list <- mclapply(exp2_als, recalculate.scf,
                           iqtree2_path = iqtree2.2.2_path,  num_iqtree2_threads = "AUTO", 
                           num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "MFP",
                           mc.cores = num_cores)
# Save output dataframe
exp2_df <- as.data.frame(do.call("rbind", exp2_scfl_list))
exp2_df_name <- paste0(results_directory, "exp2_iqtree2.2.2_scfl_results.csv")
write.csv(exp2_df, exp2_df_name, row.names = FALSE)


