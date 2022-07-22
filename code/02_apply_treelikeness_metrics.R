# /caitlinch/treelikeness-metrics/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
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
#                      <- To check both exp1 and exp2, set <- c("exp1", "exp2")

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
  num_cores <- 50
}

# Control variables
run_exp1 <- FALSE
run_exp2 <- FALSE
rerun_missing_runs <- TRUE
rerun_experiment_ids <- c("exp2")


#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the test statistics to the alignment within each directory

if (run_exp1 == FALSE){
  ## For experiment 1:
  # Open output df and get names of alignments
  exp1_op_file <- paste0(results_directory, grep("rerun", grep("exp1", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  exp1_op_df <- read.csv(exp1_op_file, stringsAsFactors = FALSE)
  # Exp1 encountering errors in all cores. Not running properly. Remove all alignments with substitution rate 1e-04 and 0.001 (too many identical sequences)
  exp1_op_df <- exp1_op_df[(exp1_op_df$tree_depth != 1e-04 & exp1_op_df$tree_depth != 1e-03),]
  # Get list of alignments
  exp1_als <- exp1_op_df$output_alignment_file
  # Apply treelikeness metrics to all alignments 
  mclapply(exp1_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
           supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
           num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
           distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
           tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
           sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = TRUE,
           redo = FALSE,
           mc.cores = num_cores)
  
  # Collect and collate results
  exp1_list <- mclapply(exp1_als, collate.treelikeness.results, mc.cores = num_cores)
  # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
  keep_indexes <- which(!sapply(exp1_list, is.null))
  exp1_list_filtered <- exp1_list[keep_indexes]
  # Save output dataframe
  exp1_df <- as.data.frame(do.call("rbind", exp1_list_filtered))
  exp1_df_name <- paste0(results_directory, "exp1_treelikeness_metrics_collated_results.csv")
  write.csv(exp1_df, exp1_df_name, row.names = FALSE)
}

if (run_exp2 == TRUE){
  ## For experiment 2:
  # Open output df and get names of alignments
  exp2_op_file <- paste0(results_directory, grep("rerun", grep("exp2", grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  exp2_op_df <- read.csv(exp2_op_file, stringsAsFactors = FALSE)
  # Get list of alignments
  exp2_als <- exp2_op_df$output_alignment_file
  # Apply treelikeness metrics to all alignments 
  mclapply(exp2_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
           supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
           num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
           distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
           tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
           sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = FALSE,
           redo = TRUE,
           mc.cores = num_cores)
  
  # Collect and collate results
  exp2_list <- mclapply(exp2_als, collate.treelikeness.results, mc.cores = num_cores)
  # Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
  keep_indexes <- which(!sapply(exp2_list, is.null))
  exp2_list_filtered <- exp2_list[keep_indexes]
  # Save output dataframe
  exp2_df <- as.data.frame(do.call("rbind", exp2_list_filtered))
  exp2_df_name <- paste0(results_directory, "exp2_treelikeness_metrics_collated_results.csv")
  write.csv(exp2_df, exp2_df_name, row.names = FALSE)
}



#### 4. Print list of alignments that do not have treelikeness results files ####
# Compare ids with those in parameters csv to determine if there are any incomplete/missing alignments

# Extract all filenames from results folder
results_files <- list.files(results_directory)

## For experiments 1 and 2 (simulations)
exp_ids <- rerun_experiment_ids
for (e in exp_ids){
  # Open parameters csv
  e_params_file <- paste0(results_directory, grep("rerun", grep(e, grep("parameters", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  e_params_df <- read.csv(e_params_file, stringsAsFactors = FALSE)
  # Open results csv
  e_results_file <- paste0(results_directory, grep("rerun", grep(e, grep("treelikeness_metrics_collated_results", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  e_results_df <- read.csv(e_results_file, stringsAsFactors = FALSE)
  # Open output paths csv
  e_op_file <- paste0(results_directory, grep("rerun", grep(e, grep("file_output_paths", results_files, value = TRUE), value = TRUE), value = TRUE, invert = TRUE))
  e_op_df <- read.csv(e_op_file, stringsAsFactors = FALSE)
  
  # For experiment 1, remove rows with substitution rates that are too low
  if (e == "exp1"){
    e_params_df  <- e_params_df[(e_params_df$tree_depth != 1e-04 & e_params_df$tree_depth != 1e-03), ]
    e_results_df <- e_results_df[(e_results_df$tree_depth != 1e-04 & e_results_df$tree_depth != 1e-03), ]
    e_op_df      <- e_op_df[(e_op_df$tree_depth != 1e-04 & e_op_df$tree_depth != 1e-03), ]
  }
  
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
        # Set which alignments to rerun
        rerun_al_paths <- missing_als_df$output_alignment_file
        # Run missing alignments
        mclapply(rerun_al_paths, treelikeness.metrics.simulations, 
                 iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                 supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
                 num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                 distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                 tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
                 sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = FALSE,
                 redo = FALSE,
                 mc.cores = num_cores)
 
        # Collect and collate results
        e_list <- mclapply(rerun_al_paths, collate.treelikeness.results, mc.cores = num_cores)
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
  
} # end for (e in exp_ids)





