#### Tests ####
# test params
alignment_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00004_0050_001_1/exp1_00004_0050_001_1_output_alignment.fa"
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
return_collated_data = TRUE


# treelikeness.metrics.simulations(alignment_path,
#                                  iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
#                                  supply_number_of_taxa = FALSE, number_of_taxa = NA, num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100,
#                                  iqtree_substitution_model = "JC", distance_matrix_substitution_method = "JC69", num_phylogemetric_threads = NA,
#                                  tree_proportion_remove_trivial_splits = TRUE, sequence_format = "DNA", return_collated_data = TRUE)

# # To check number of completed runs:
d <- "/data/caitlin/treelikeness_metrics/exp_1/"
all_files <- list.files(d, recursive = TRUE)
all_files <- grep("1e-04", all_files, value = TRUE, invert = TRUE)
all_files <- grep("0.001", all_files, value = TRUE, invert = TRUE)
aln_files <- grep("_output_alignment", all_files, value = TRUE)
als <- grep(".fa.",aln_files, value = TRUE, invert = TRUE)
col_files <- grep("collated", all_files, value = TRUE)
tl_files <- grep("treelikeness_results", all_files, value = TRUE)
print(length(als))
print(length(tl_files))
print(paste0(length(tl_files)/length(als)*100, " % Experiment 1 alignments completed"))
# To get unrun alignments:
al_dirs <- list.dirs(d, full.names = FALSE)[2:length(list.dirs(d, full.names = FALSE))]
al_dirs <- grep("1e-04", al_dirs, value = TRUE, invert = TRUE)
al_dirs <- grep("0.001", al_dirs, value = TRUE, invert = TRUE)
all_tl_files <- paste0(d, al_dirs, "/", al_dirs, "_treelikeness_results.csv")
missing_tl_files <- all_tl_files[!file.exists(all_tl_files)]
missing_als <- gsub("_treelikeness_results.csv", "_output_alignment.fa", missing_tl_files)
write(missing_als, file = "/data/caitlin/treelikeness_metrics/exp1_alignments_to_run.txt")
# To check number of completed runs:
d <- "/data/caitlin/treelikeness_metrics/exp_2/"
all_files <- list.files(d, recursive = TRUE)
aln_files <- grep("_output_alignment", all_files, value = TRUE)
als <- grep(".fa.",aln_files, value = TRUE, invert = TRUE)
col_files <- grep("collated", all_files, value = TRUE)
tl_files <- grep("treelikeness_results", all_files, value = TRUE)
print(length(als))
print(length(tl_files))
print(paste0(length(tl_files)/length(als)*100, " % Experiment 2 alignments completed"))
# To get unrun alignments:
al_dirs <- list.dirs(d, full.names = FALSE)[2:length(list.dirs(d, full.names = FALSE))]
all_tl_files <- paste0(d, al_dirs, "/", al_dirs, "_treelikeness_results.csv")
missing_tl_files <- all_tl_files[!file.exists(all_tl_files)]
missing_als <- gsub("_treelikeness_results.csv", "_output_alignment.fa", missing_tl_files)
write(missing_als, file = "/data/caitlin/treelikeness_metrics/exp2_alignments_to_run.txt")

