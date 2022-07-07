# /caitlinch/treelikeness-metrics/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).



#### 1. Set parameters ####
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# fast_TIGER_path         <- Path to fast TIGER executable.
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above
# num_cores               <- Number of parallel threads to use at once

# run_exp1 <- Whether to apply the treelikeness test statistics to the first set of alignments
# run_exp2 <- Whether to apply the treelikeness test statistics to the second set of alignments

run_location = "local"
if (run_location == "local"){
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  iqtree2_path <- "iqtree2.2-beta"
  splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
  fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
  num_cores <- 1
} else if (run_location == "soma"){
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  repo_directory <- "/data/caitlin/treelikeness_metrics/code/"
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  num_cores <- 50
}

run_exp1 <- FALSE
run_exp2 <- FALSE
collect_missing_runs <- FALSE
rerun_missing_runs <- TRUE



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/"))



#### 3. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the test statistics to the alignment within each directory

if (run_exp1 == FALSE){
  # For experiment 1:
  exp1_all_files <- paste0(exp_folders[1], list.files(exp_folders[1], recursive = TRUE))
  exp1_aln_files <- grep("_output_alignment", exp1_all_files, value = TRUE)
  exp1_als <- grep(".fa.", exp1_aln_files, value = TRUE, invert = TRUE)
  # Exp1 encountering errors in all cores. Not running properly. Remove all alignments with substitution rate 1e-04 and 0.001 (too many identical sequences)
  exp1_als <- grep("1e-04", exp1_als, value = TRUE, invert = TRUE)
  exp1_als <- grep("0.001", exp1_als, value = TRUE, invert = TRUE)
  # Apply treelikeness metrics to all alignments with substitution rate of 1e-02 (0.01) or higher
  exp1_list <- mclapply(exp1_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                        supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
                        num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                        distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                        tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
                        sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = TRUE,
                        redo = FALSE,
                        mc.cores = num_cores)
  exp1_df <- as.data.frame(do.call("rbind", exp1_list))
  exp1_df_name <- paste0(local_directory, "exp1_treelikeness_metrics_results.csv")
  write.csv(exp1_df, exp1_df_name, row.names = FALSE)
}

if (run_exp2 == TRUE){
  # For experiment 2:
  exp2_all_files <- paste0(exp_folders[2], list.files(exp_folders[2], recursive = TRUE))
  exp2_aln_files <- grep("_output_alignment", exp2_all_files, value = TRUE)
  exp2_als <- grep(".fa.", exp2_aln_files, value = TRUE, invert = TRUE)
  exp2_list <- mclapply(exp2_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                        supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
                        num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                        distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                        tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
                        sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = FALSE,
                        redo = TRUE,
                        mc.cores = num_cores)
  exp2_df <- as.data.frame(do.call("rbind", exp2_list))
  exp2_df_name <- paste0(local_directory, "exp2_treelikeness_metrics_results.csv")
  write.csv(exp2_df, exp2_df_name, row.names = FALSE)
}



#### 4. Print list of alignments that do not have treelikeness results files ####
if (collect_missing_runs == TRUE){
  # Collect experiment 1 missing files
  exp1_all_files <- paste0(exp_folders[1], list.files(exp_folders[1], recursive = TRUE))
  exp1_aln_files <- grep("_output_alignment", exp1_all_files, value = TRUE)
  exp1_als <- grep(".fa.", exp1_aln_files, value = TRUE, invert = TRUE)
  exp1_als <- grep("1e-04", exp1_als, value = TRUE, invert = TRUE) # Remove all alignments with substitution rate 1e-04 (too many identical sequences)
  exp1_als <- grep("0.001", exp1_als, value = TRUE, invert = TRUE) # Remove all alignments with substitution rate 0.001 (too many identical sequences)
  exp1_apply_metrics <- exp1_als[which(file.exists(gsub("output_alignment.fa", "treelikeness_results.csv", exp1_als)) == FALSE)]
  # Output list of experiment 1 missing files
  write(exp1_apply_metrics, paste0(local_directory, "exp1_missing_results.csv"))
  
  # Collect experiment 2 missing files
  exp2_all_files <- paste0(exp_folders[2], list.files(exp_folders[2], recursive = TRUE))
  exp2_aln_files <- grep("_output_alignment", exp2_all_files, value = TRUE)
  exp2_als <- grep(".fa.", exp2_aln_files, value = TRUE, invert = TRUE)
  exp2_apply_metrics <- exp2_als[which(file.exists(gsub("output_alignment.fa", "treelikeness_results.csv", exp2_als)) == FALSE)]
  # Output list of experiment 2 missing files
  write(exp2_apply_metrics, paste0(local_directory, "exp2_missing_results.csv"))
}



#### 5. Cross check parameters csv with treelikeness results csv ####
# Compare ids with those in parameters csv to determine if there are any incomplete/missing alignments

# Check for results df
results_dir <- paste0(local_directory, "01_results/")
if (dir.exists(results_dir) == FALSE){dir.create(results_dir)}

# Extract all filenames from results folder
results_files <- list.files(results_dir)

## For experiments 1 and 2 (simulations)
exp_ids <- c("exp1", "exp2")
for (e in exp_ids){
  # Open parameters csv
  e_params_file <- paste0(results_dir, grep(e, grep("parameters", results_files, value = TRUE), value = TRUE))
  e_params_df <- read.csv(e_params_file, stringsAsFactors = FALSE)
  # Open results csv
  e_results_file <- paste0(results_dir, grep(e, grep("treelikeness_metrics_results", results_files, value = TRUE), value = TRUE))
  e_results_df <- read.csv(e_results_file, stringsAsFactors = FALSE)
  # Open output paths csv
  e_op_file <- paste0(results_dir, grep(e, grep("file_output_paths", results_files, value = TRUE), value = TRUE))
  e_op_df <- read.csv(e_op_file, stringsAsFactors = FALSE)
  
  # For experiment 1, remove rows with substitution rates that are too low
  if (e == "exp1"){
    e_params_df <- e_params_df[(e_params_df$tree_depth != 1e-04 & e_params_df$tree_depth != 1e-03), ]
    e_results_df <- e_results_df[(e_results_df$tree_depth != 1e-04 & e_results_df$tree_depth != 1e-03), ]
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
      # Get ids of alignments missing from params df (safety check - shouldn't be possible) and save
      missing_ids <- params_ids[!(params_ids %in% results_ids)]
      missing_ids_file <- paste0(results_dir, e, "_uids_", output_id, ".txt") 
      write(missing_its, file = missing_ids_file)
    } else if (length(results_ids) > length(params_ids)){
      # Set output id
      output_id <- "missing.from.params"
      # Get ids of alignments missing from results df (unrun - need to run) and save
      missing_ids <- results_ids[!(results_ids %in% params_ids)]
      missing_ids_file <- paste0(results_dir, e, "_uids_", output_id, ".txt") 
      write(missing_its, file = missing_ids_file)
      # Make dataframe consisting of only missing alignments that need running and save
      missing_als_df <- e_op_df[e_op_df$uid %in% missing_ids, ]
      missing_als_file <- paste0(results_dir, e, "_parameters_rerun_", output_id, ".csv")
      write.csv(missing_als_df, file = missing_als_file, row.names = TRUE)
      
      if (rerun_missing_runs == TRUE){
        # Run missing alignments
        mclapply(missing_als_df$output_alignment_file, treelikeness.metrics.simulations, 
                 iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                 supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
                 num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                 distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                 tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
                 sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = TRUE,
                 redo = FALSE,
                 mc.cores = num_cores)
        # Collate all alignments 
        e_list <- mclapply(missing_als_df$output_alignment_file, treelikeness.metrics.simulations, 
                           iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                           supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
                           num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                           distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                           tree_proportion_remove_trivial_splits = TRUE, run_splitstree_for_tree_proportion = TRUE,
                           sequence_format = "DNA", return_collated_data = TRUE, apply.TIGER = TRUE,
                           redo = FALSE,
                           mc.cores = num_cores)
        e_rerun_df <- as.data.frame(do.call("rbind", e_list))
        e_rerun_df_name <- paste0(local_directory, e, "_treelikeness_metrics_results_collated.csv")
        write.csv(e_rerun_df, e_rerun_df_name, row.names = FALSE)
      } # end (rerun_missing_runs == TRUE)
      
    } # end (length(results_ids) > length(params_ids))
    
  } # end (all_uids_present == FALSE)
  
} # end for (e in exp_ids)





