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

run_location = "dayhoff"
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
} else if (run_location == "dayhoff"){
  # Directories
  local_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_metrics/"
  results_directory <- local_directory
  repo_directory <- "/mnt/data/dayhoff/home/u5348329/treelikeness_metrics/"

  # Executable paths
  iqtree2_path <- "/mnt/data/dayhoff/home/u5348329/treelikeness_metrics/software/iqtree-2.4.0-Linux-intel/bin/iqtree2"
  splitstree_path <- "/mnt/data/dayhoff/home/u5348329/splitstree4/SplitsTree"
  phylogemetric_path <- "/mnt/data/dayhoff/home/u5348329/.local/bin/phylogemetric"
  fast_TIGER_path <- "/mnt/data/dayhoff/home/u5348329/treelikeness_metrics/software/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"

  # Run parameters
  num_cores <- 30
}

# Control variables
run_exp1 <- TRUE
run_exp3 <- FALSE



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
if (run_exp1 == TRUE){
  ## For experiment 1:
  # Extract all file names from results folder
  results_files <- list.files(results_directory)
  # Open output df and get names of alignments
  exp1_op_file <- paste0(results_directory,
                         grep(
                           "rerun",
                           grep(
                             "exp1",
                             grep("file_output_paths", results_files, value = TRUE),
                             value = TRUE
                           ),
                           value = TRUE,
                           invert = TRUE
                         ))
  exp1_op_df <- read.csv(exp1_op_file, stringsAsFactors = FALSE)
  # Reduce to only first 20 reps for each set of simulation parameters
  exp1_op_df <- exp1_op_df[which(exp1_op_df$num_reps <= 20), ]
  # Extract alignments
  exp1_als <- exp1_op_df$output_alignment_file
  # # Apply treelikeness metrics to all alignments
  # mclapply(
  #   exp1_als,
  #   treelikeness.metrics.simulations,
  #   iqtree2_path,
  #   splitstree_path,
  #   phylogemetric_path,
  #   fast_TIGER_path,
  #   supply_number_of_taxa = FALSE,
  #   number_of_taxa = NA,
  #   num_iqtree2_threads = 1,
  #   num_iqtree2_scf_quartets = 100,
  #   iqtree_substitution_model = "JC",
  #   distance_matrix_substitution_method = "JC69",
  #   num_phylogemetric_threads = NA,
  #   tree_proportion_remove_trivial_splits = TRUE,
  #   run_splitstree_for_tree_proportion = TRUE,
  #   sequence_format = "DNA",
  #   apply.TIGER = TRUE,
  #   redo = TRUE,
  #   mc.cores = num_cores
  # )
  ## FORCE RERUN OF MISSING REPS
  missing_reps <- c(1712, 1742, 1772, 1802, 1832, 1862, 1892, 1922, 1952, 1982,
                    2012, 2042, 2072, 2102, 2132, 2162, 2192, 2222, 2252, 2282,
                    2312, 2342, 2372, 2402, 2432, 2462, 2492, 2522, 2552, 2582,
                    2612, 2642, 2672, 2702, 2732, 2762, 2792, 2822, 2852, 2882,
                    2912, 2942, 2972, 3002, 3032, 3062, 3092, 3122, 3152, 3182,
                    3212, 3242, 3272, 3302, 3332, 3362, 3392, 3422, 3452, 3482,
                    3512, 3542, 3572, 3602, 3632, 3662, 3692, 3722, 3752, 3782,
                    3812, 3842, 3872, 3902, 3932, 3962, 3992, 4022, 4052, 4082,
                    4112, 4142, 4172, 4202, 4232, 4262, 4292, 4322, 4352, 4382,
                    4412, 4442, 4472, 4502, 4532, 4562, 4592, 4622, 4652, 4682,
                    4712, 4742, 4772, 4802, 4832, 4862, 4892, 4922, 4952, 4982,
                    5012, 5042, 5072, 5102, 5132, 5162, 5192, 5222, 5252, 5282,
                    5312, 5342, 5372, 5402, 5432, 5462, 5492, 5522, 5552, 5582,
                    5612, 5642, 5672, 5702, 5732, 5762, 5792, 5822, 5852, 5882,
                    5912, 5942, 5972, 6002, 6032, 6062, 6092, 6122, 6152, 6182,
                    6212, 6242, 6272, 6302, 6332, 6362, 6392, 6422, 6452, 6482,
                    6512, 6542, 6572, 6602, 6632, 6662, 6692, 6722, 6752, 6782,
                    6812, 6842, 6872, 6902, 6932, 6962, 6992, 7022, 7052, 7082,
                    7112, 7142, 7172, 7202, 7232, 7262, 7292, 7322, 7352, 7382,
                    7412, 7442, 7472)
  mclapply(
    exp1_als[missing_reps],
    treelikeness.metrics.simulations,
    iqtree2_path,
    splitstree_path,
    phylogemetric_path,
    fast_TIGER_path,
    supply_number_of_taxa = FALSE,
    number_of_taxa = NA,
    num_iqtree2_threads = 1,
    num_iqtree2_scf_quartets = 100,
    iqtree_substitution_model = "JC",
    distance_matrix_substitution_method = "JC69",
    num_phylogemetric_threads = NA,
    tree_proportion_remove_trivial_splits = TRUE,
    run_splitstree_for_tree_proportion = TRUE,
    sequence_format = "DNA",
    apply.TIGER = TRUE,
    redo = TRUE,
    mc.cores = num_cores
  )
  # Identify unrun files
  exp1_results_files <- gsub("_output_alignment.fa", "_treelikeness_results.csv", exp1_als)
  exp1_results_file_sizes <- unlist(lapply(exp1_results_files, function(x){file.size(x)}))
  missing_result_ids <- which(exp1_results_file_sizes == 0)
  if (length(missing_result_ids)){
    print(paste0("Alignment indexes needing reruns:", paste(missing_result_ids, collapse = ", ")))
  }
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


