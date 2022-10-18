# /caitlinch/treelikeness-metrics/code/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will the tree proportion test



#### 1. Set parameters ####
## Directories
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory       <- Directory where collated results from the treelikeness test statistics will be saved
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Run parameters
# num_cores               <- Number of parallel threads to use at once


run_location = "soma"
if (run_location == "local"){
  # Directories
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  results_directory <- paste0(local_directory, "01_results/")
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  
  # Run parameters
  num_cores <- 1
  
} else if (run_location == "soma"){
  # Directories
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  results_directory <- local_directory
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  
  # Run parameters
  num_cores <- 20
}


#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_tree_proportion.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))


#### 3. Apply tree proportion
# Set experiment id 
exp_ids <- c("exp1", "exp2")

# Extract all filenames from results folder
results_files <- list.files(results_directory)

for (e in exp_ids){
  # Get csv files
  e_params_file <- paste0(results_directory, grep("rerun", grep(e, grep("parameters", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  e_results_file <- paste0(results_directory, grep("rerun", grep(e, grep("treelikeness_metrics_collated_results", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  e_op_file <- paste0(results_directory, grep("rerun", grep(e, grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  # Open output paths csv
  e_op_df <- read.csv(e_op_file, stringsAsFactors = FALSE)
  
  # Extract the list of all alignments
  all_alignments <- e_op_df$output_alignment_file
  test_als <- all_alignments[1:10]
  
  # Apply the tree proportion function
  # tree.proportion(alignment_path, sequence_format = "DNA", remove_trivial_splits = TRUE, network.method = "SplitsTree4" , splitstree_path = NA, dist.ml.model = NA)
  tp_list <- mclapply(all_alignments, tree.proportion.output.csv, sequence_format = "DNA", remove_trivial_splits = TRUE, 
                      network.method = "SplitsTree4", splitstree_path = splitstree_path, dist.ml.model = NA, 
                      mc.cores = num_cores)
  
  # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
  keep_indexes <- which(!sapply(tp_list, is.null))
  tp_list_filtered <- tp_list[keep_indexes]
  # Save output dataframe
  tp_df <- as.data.frame(do.call("rbind", tp_list_filtered))
  tp_df_name <- paste0(results_directory, e, "_TreeProportion_collated_results.csv")
  write.csv(tp_df, tp_df_name, row.names = FALSE)
}




#### 4. Collate tree proportion results
for (e in exp_ids){
  # Get list of all files
  all_uids <- e_op_df$uid
  all_dirs <- dirname(e_op_df$output_alignment_file)
  tp_output_files <- paste0(all_dirs, "/", all_uids, "_TreeProportion_results.csv")
  complete_tp_files <- tp_output_files[which(file.exists(tp_output_files)==TRUE)]
  tp_op_list <- lapply(complete_tp_files, read.csv)
  tp_op_df <- do.call(rbind, tp_op_list)
  tp_op_df_name <- paste0(results_directory, e, "_treeProportion_completeRuns_collated_results.csv")
  write.csv(tp_op_df, tp_op_df_name, row.names = FALSE)
}
