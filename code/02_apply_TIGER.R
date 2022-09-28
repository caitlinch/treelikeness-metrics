# /caitlinch/treelikeness-metrics/code/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will apply various tests for treelikeness to simulated alignments
# This program requires fast TIGER



#### 1. Set parameters ####
## Directories
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory       <- Directory where collated results from the treelikeness test statistics will be saved
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# fast_TIGER_path         <- Path to fast TIGER executable


## Run parameters
# num_cores               <- Number of parallel threads to use at once


run_location = "soma"
if (run_location == "local"){
  # Directories
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  results_directory <- paste0(local_directory, "01_results/")
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  results_directory <- local_directory
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  
  # Run parameters
  num_cores <- 30
}


#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))


#### 3. Apply TIGER
# Set experiment id 
e <- "exp2"

# Extract all filenames from results folder
results_files <- list.files(results_directory)

# Get csv files
e_params_file <- paste0(results_directory, grep("rerun", grep(e, grep("parameters", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
e_results_file <- paste0(results_directory, grep("rerun", grep(e, grep("treelikeness_metrics_collated_results", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
e_op_file <- paste0(results_directory, grep("rerun", grep(e, grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
# Open output paths csv
e_op_df <- read.csv(e_op_file, stringsAsFactors = FALSE)

# Extract the list of all alignments
all_alignments <- e_op_df$output_alignment_file

# Apply the TIGER function
# tiger.empirical(alignment_path, fast_TIGER_path, sequence_format = "DNA")
tiger_list <- lapply(all_alignments, tiger.empirical, fast_TIGER_path, sequence_format = "DNA")

# Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
keep_indexes <- which(!sapply(tiger_list, is.null))
tiger_list_filtered <- tiger_list[keep_indexes]
# Save output dataframe
tiger_df <- as.data.frame(do.call("rbind", tiger_list_filtered))
tiger_df_name <- paste0(results_directory, e, "_treelikeness_metrics_collated_results.csv")
write.csv(tiger_df, tiger_df_name, row.names = FALSE)

