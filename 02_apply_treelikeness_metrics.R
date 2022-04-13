# /caitlinch/treelikeness_metrics/02_apply_treelikeness_metrics.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).



#### 1. Open packages ####



#### 2. Set parameters ####
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# fast_TIGER_path         <- Path to fast TIGER executable.
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
iqtree2_path <- "iqtree2.2-beta"
splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"



#### 3. Source functions from caitlinch/treelikeness_metrics ####
source(paste0(repo_directory, "func_metrics.R"))



#### 4. Apply tests for treelikeness to each simulated alignment ####
# For each experiment, get the list of directories within that experiment folder and apply the test statistics to the alignment within each directory
# Sample function call:
# treelikeness.metrics.simulations(alignment_path, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
#                                  num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
#                                  distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
#                                  tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA")

# Find the folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/"))

# For experiment 1:
exp1_runs <- paste0(exp_folders[1], list.files(exp_folders[1]), "/") # should be list of alignments
exp1_list <- lapply(exp1_runs, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                    supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                    iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA, 
                    tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA")
exp1_df <- as.data.frame(do.call("rbind", exp1_list))

# For experiment 2:
exp2_runs <- paste0(exp_folders[2], list.files(exp_folders[2]), "/")  # should be list of alignments
exp2_list <- lapply(exp2_runs, treelikeness.metrics.simulations, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path, 
                    supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, 
                    iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA, 
                    tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA")
exp2_df <- as.data.frame(do.call("rbind", exp2_list))




#### Tests ####
# test params
iqtree2_path <- "iqtree2.2-beta"
splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
supply_number_of_taxa = FALSE
number_of_taxa = NA
num_iqtree2_threads = "AUTO"
num_iqtree2_scf_quartets = 100
iqtree_substitution_model = "JC"
distance_matrix_substitution_method = "JC69"
num_phylogemetric_threads = NA
tree_proportion_remove_trivial_splits = TRUE
sequence_format = "DNA"


