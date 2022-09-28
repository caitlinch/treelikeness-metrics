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
  iqtree2_path <- "iqtree2.2-beta"
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
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
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

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_2/"))



#### 3. Apply TIGER
# Set experiment id 
e <- "exp2"

# Extract all filenames from results folder
results_files <- list.files(results_directory)

# Open parameters csv
e_params_file <- paste0(results_directory, grep("rerun", grep(e, grep("parameters", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
e_params_df <- read.csv(e_params_file, stringsAsFactors = FALSE)
# Open results csv
e_results_file <- paste0(results_directory, grep("rerun", grep(e, grep("treelikeness_metrics_collated_results", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
e_results_df <- read.csv(e_results_file, stringsAsFactors = FALSE)
# Open output paths csv
e_op_file <- paste0(results_directory, grep("rerun", grep(e, grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
e_op_df <- read.csv(e_op_file, stringsAsFactors = FALSE)

# Extract the list of all alignments
all_alignments <- e_op_df$output_alignment_file

# Check that all unique ids have a match
all_uids_present = setequal(e_results_df$uid, e_params_df$uid)
# If some uids are missing, report and print which
if (all_uids_present == FALSE){
  # Get both lists of unique identifiers
  params_ids <- e_params_df$uid
  results_ids <- e_results_df$uid
  # Check which is longer and save the list of uids missing from one dataframe
  if (length(params_ids) > length(results_ids)){
    # Set output id
    output_id <- "missing.from.results"
    
    # Get ids of alignments missing from results df (unrun - need to run) and save
    missing_ids <- params_ids[!(params_ids %in% results_ids)]
    missing_ids_file <- paste0(results_directory, e, "_uids_", output_id, ".txt")
    write(missing_ids, file = missing_ids_file)
    
    # Make dataframe consisting of only missing alignments that need running and save
    missing_als_df <- e_op_df[e_op_df$uid %in% missing_ids, ]
    missing_als_file <- paste0(results_directory, e, "_parameters_rerun_", output_id, ".csv")
    write.csv(missing_als_df, file = missing_als_file, row.names = TRUE)
    
    if (rerun_missing_runs == TRUE){
      # Set which alignments to rerun (remove alignments with certain parameters that are too slow or have insufficient information to run)
      if (e == "exp1"){
        # For experiment 1, remove rows with substitution rates that are too low
        missing_als_df  <- missing_als_df[(missing_als_df$tree_depth != 1e-04 & missing_als_df$tree_depth != 1e-03), ]
      }
      
      # Extract vector of alignments to rerun
      rerun_al_paths <- missing_als_df$output_alignment_file
      
      # Run missing alignments
      mclapply(rerun_al_paths, treelikeness.metrics.simulations,
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
      
      # Collect and collate ALL results
      if (e == "exp1"){
        e_list <- mclapply(all_alignments, collate.treelikeness.results, experiment_number = 1, mc.cores = num_cores)  
      } else if (e == "exp2"){
        e_list <- mclapply(all_alignments, collate.treelikeness.results, experiment_number = 2, mc.cores = num_cores)
      }
      
      # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
      keep_indexes <- which(!sapply(e_list, is.null))
      e_list_filtered <- e_list[keep_indexes]
      # Save output dataframe
      e_df <- as.data.frame(do.call("rbind", e_list_filtered))
      e_df_name <- paste0(results_directory, e, "_treelikeness_metrics_collated_results.csv")
      write.csv(e_df, e_df_name, row.names = FALSE)
      
    } else if (length(results_ids) > length(params_ids)){
      # Set output id
      output_id <- "missing.from.params"
      
      # Get ids of alignments missing from params df (safety check - shouldn't be possible) and save
      missing_ids <- results_ids[!(results_ids %in% params_ids)]
      missing_ids_file <- paste0(results_directory, e, "_uids_", output_id, ".txt") 
      write(missing_its, file = missing_ids_file)
      
    } # end (rerun_missing_runs == TRUE)
    
  } # end (length(results_ids) > length(params_ids))
  
} # end (all_uids_present == FALSE)
