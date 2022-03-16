## Open packages


## Parameters for applying metrics
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# total_alignment_length  <- Total length of concatenated alignments (we chose 10000).
# sequence_type           <- Sequence type for simulation (we chose "DNA").
# taxa_vec                <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# tree_depth              <- Single value or vector of values for total tree length
# r_vec                   <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models      <- model of sequence evolution for Alisim 
# alisim_gene_tree_length <- gene-specific tree length for Alisim

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
