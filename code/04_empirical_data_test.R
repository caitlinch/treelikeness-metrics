# caitlinch/treelikeness-metrics/code/04_empirical_data_test.R
# Caitlin Cherryh 2023

## Script summary:
# This program will apply various tests for treelikeness to empirical alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

## Empirical dataset:
# Paper: 
# Data repository: 


#### 1. Set parameters ####
## Directories
# run_directory           <- Directory where alignments will be saved/treelikeness metrics will be run
# data_directory          <- Location of empirical dataset
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)
# fast_TIGER_path         <- Path to fast TIGER executable
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "soma"
if (run_location == "local"){
  # Directories
  run_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical/"
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
  run_directory <- "/data/caitlin/treelikeness_metrics/"
  data_directory <- ""
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  
  # Run parameters
  num_cores <- 30
}

# For testing treelikeness metrics
alignment_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical/ITS_GB/ITS_GB.nxs"
supply_number_of_taxa = TRUE
number_of_taxa = 37
num_iqtree2_threads = "AUTO"
num_iqtree2_scf_quartets = 100
iqtree_substitution_model = "JC"
distance_matrix_substitution_method = "JC69"
num_phylogemetric_threads = NA
tree_proportion_remove_trivial_splits = TRUE
run_splitstree_for_tree_proportion = FALSE
sequence_format = "DNA"
apply.TIGER = TRUE
redo = FALSE




#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Find the folders of simulated alignments
exp_folder <- paste0(run_directory, "emp/")



#### 3. Apply tests for treelikeness to each gene within the alignment ####
# Open output df and get names of alignments
emp_op_file <- paste0(run_directory, grep("rerun", grep("emp", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
emp_op_df <- read.csv(emp_op_file, stringsAsFactors = FALSE)
# Get list of alignments
emp_als <- emp_op_df$output_alignment_file
# Apply treelikeness metrics to all alignments 
mclapply(emp_als, treelikeness.metrics.empirical,
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
emp_list <- mclapply(emp_als, collate.treelikeness.results, experiment_number = 1, mc.cores = num_cores)
# Save output dataframe
emp_df <- as.data.frame(do.call("rbind", emp_list))
emp_df_name <- paste0(run_directory, "emp_treelikeness_metrics_collated_results.csv")
write.csv(emp_df, emp_df_name, row.names = FALSE)





