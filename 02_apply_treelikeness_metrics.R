## Open packages


## Parameters for applying metrics
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# fast_TIGER_path         <- Path to fast TIGER executable.
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above
# netmake_path            <- Path to netmake executable within SPECTRE
# netme_path              <- Path to netme executable within SPECTRE

local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
iqtree2_path <- "iqtree2.2-beta"
fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
splitstree_path <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
netmake_path <- "/Applications/Spectre.app/Contents/MacOS/netmake"
netme_path <- "/Applications/Spectre.app/Contents/MacOS/netme"


## Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "func_metrics.R"))


## Find the three folders of simulated alignments
exp_folders <- paste0(local_directory, c("exp_1/", "exp_2/", "exp_3/"))

# For each experiment, get the list of folders within that experiment folder
exp1_runs <- paste0(exp_folders[1], list.files(exp_folders[1]), "/")
replicate_folder <- exp1_runs[1]
replicate_folder <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00010_0010_001_1/"

# test params
num_iqtree2_threads = "AUTO"
num_iqtree2_scf_quartets = 100
iqtree_substitution_model = "JC"
delta_plot_substitution_method = "JC69"
num_phylogemetric_threads = NA
sequence_format = "DNA"

treelikeness.metrics.simulations(replicate_folder, iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
                                 num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100, iqtree_substitution_model = "JC", 
                                 delta_plot_substitution_method = "JC69", num_phylogemetric_threads = NA, sequence_format = "DNA")

