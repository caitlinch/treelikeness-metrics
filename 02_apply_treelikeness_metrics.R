# /caitlinch/treelikeness_metrics/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).



#### 1. Open packages ####
library(parallel)



#### 2. Set parameters ####
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# fast_TIGER_path         <- Path to fast TIGER executable.
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above
# num_cores               <- Number of parallel threads to use at once

run_location = "soma"
if (run_location == "local"){
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
  iqtree2_path <- "iqtree2.2-beta"
  splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
  fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
  num_cores <- 30
} else if (run_location == "soma"){
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  repo_directory <- "/data/caitlin/treelikeness_metrics/code/"
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  num_cores <- 1
}


#### 3. Source functions from caitlinch/treelikeness_metrics ####
source(paste0(repo_directory, "func_metrics.R"))



#### 4. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the test statistics to the alignment within each directory
# Sample function call:
# treelikeness.metrics.simulations(alignment_path, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
#                                  supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", 
#                                  num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
#                                  distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
#                                  tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA",
#                                  return_collated_data = FALSE)

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/"))

# For experiment 1:
exp1_all_files <- paste0(exp_folders[1], list.files(exp_folders[1], recursive = TRUE))
exp1_aln_files <- grep("_output_alignment", exp1_all_files, value = TRUE)
exp1_als <- grep(".fa.", exp1_aln_files, value = TRUE, invert = TRUE)
exp1_list <- mclapply(exp1_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                      supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                      iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA, 
                      tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA", return_collated_data = TRUE,
                      mc.cores = num_cores)
exp1_df <- as.data.frame(do.call("rbind", exp1_list))

# For experiment 2:
exp2_all_files <- paste0(exp_folders[2], list.files(exp_folders[2], recursive = TRUE))
exp2_aln_files <- grep("_output_alignment", exp2_all_files, value = TRUE)
exp2_als <- grep(".fa.", exp2_aln_files, value = TRUE, invert = TRUE)
exp2_list <- mclapply(exp2_als, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
                      supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100,
                      iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
                      tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA", return_collated_data = TRUE,
                      mc.cores = num_cores)
exp2_df <- as.data.frame(do.call("rbind", exp2_list))


# Print list of alignments that do not have treelikeness results files
exp1_apply_metrics <- exp1_als[which(file.exists(gsub("output_alignment.fa", "treelikeness_results.csv", exp1_als)) == FALSE)]
write(exp1_apply_metrics, paste0(local_directory, "exp1_missing_results.csv"))
exp2_apply_metrics <- exp2_als[which(file.exists(gsub("output_alignment.fa", "treelikeness_results.csv", exp2_als)) == FALSE)]
write(exp2_apply_metrics, paste0(local_directory, "exp2_missing_results.csv"))



# #### Tests ####
# # test params
# alignment_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00004_0050_001_1/exp1_00004_0050_001_1_output_alignment.fa"
# iqtree2_path <- "iqtree2.2-beta"
# splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
# phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
# fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
# supply_number_of_taxa = FALSE
# number_of_taxa = NA
# num_iqtree2_threads = "AUTO"
# num_iqtree2_scf_quartets = 100
# iqtree_substitution_model = "JC"
# distance_matrix_substitution_method = "JC69"
# num_phylogemetric_threads = NA
# tree_proportion_remove_trivial_splits = TRUE
# sequence_format = "DNA"
# return_collated_data = TRUE
# 
# 
# treelikeness.metrics.simulations("/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00005_1000_009_1/exp1_00005_1000_009_1_output_alignment.fa",
#                                  iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
#                                  supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
#                                  iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA, 
#                                  tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA", return_collated_data = TRUE)
