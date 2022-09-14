# /caitlinch/treelikeness-metrics/code/03_empirical_data_treelikeness.R
# Caitlin Cherryh 2022

# This program will prepare a mitochondrial DNA data set for analysis and apply various tests for treelikeness
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

# Paper:
#   Emilie J Richards, Jeremy M Brown, Anthony J Barley, Rebecca A Chong, Robert C Thomson 2018,
#       "Variation Across Mitochondrial Gene Trees Provides Evidence for Systematic Error: 
#       How Much Gene Tree Variation Is Biological?", Systematic Biology, 67:5, Pages 847–860, 
#       https://doi.org/10.1093/sysbio/syy013

# Data available from: 
#    Emilie J Richards, Jeremy M Brown, Anthony J Barley, Rebecca A Chong, Robert C Thomson 2018,
#       "Data from: Variation across mitochondrial gene trees provides evidence for systematic error: 
#       how much gene tree variation is biological?", Dryad, Dataset, https://doi.org/10.5061/dryad.hj07m



#### 1. Set parameters ####
## Directories
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run
# results_directory       <- Directory where collated results from the treelikeness test statistics will be saved
# richards2011_directory  <- Directory containing mtDNA alignment files downloaded from the DataDryad (nexus_alignment_files.tar.gz)
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim)
# fast_TIGER_path         <- Path to fast TIGER executable
# phylogemetric_path      <- Path to phylogemetric executable
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "local"
if (run_location == "local"){
  # Directories
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  results_directory <- paste0(local_directory, "01_results/")
  richards2011_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_Richards2018/nexus_files/"
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
  richards2011_directory <- "/data/caitlin/treelikeness_metrics/Richards2018/"
  repo_directory <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  splitstree_path <- "/home/caitlin/splitstree4/SplitsTree"
  phylogemetric_path <- "/home/caitlin/.local/bin/phylogemetric"
  fast_TIGER_path <- "/data/caitlin/linux_executables/fast_TIGER/fast_TIGER"
  
  # Run parameters
  num_cores <- 50
}



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_empirical_datasets.R"))



#### 3. Prepare Richards2018 alignments for analysis ####
## Extract basic details about the alignments
# Specify the csv file containing basic information about the richards2018 nexus alignments
alignment_dim_file <- paste0(results_directory, "exp3_richards2018_alignment_dimensions.csv")
# If the file containing basic information about the alignments exists, open it. If not, create it.
if (file.exists(alignment_dim_file) == TRUE){
  alignment_dims <- read.csv(alignment_dim_file)
} else {
  all_alignments <- list.files(richards2011_directory, recursive = T)
  # Extend to full file path
  if (length(all_alignments) > 0){
    all_alignments <- paste0(richards2011_directory, all_alignments)
  }
  # Remove any misaligned files
  all_alignments <- grep("misalignment", all_alignments, value = T, invert = T)
  # Extract the number of characters and taxa in each alignment
  alignment_dim_list <- lapply(all_alignments, alignment.dimensions.nex)
  alignment_dims <- as.data.frame(do.call(rbind, alignment_dim_list))
  names(alignment_dims) <- c("alignment_path", "num_taxa", "num_sites")
  # Create new column with a unique identifier for each alignment
  alignment_dims$uid <- paste0(gsub("\\.nex", "", basename(alignment_dims$alignment_path)), "_copy")
  # Add dataset id as column
  alignment_dims$dataset <- "Richards2018"
  # Add column for the clade and gene
  alignment_dims$clade <- unlist(lapply(strsplit(alignment_dims$uid, "_"), `[[`, 1))
  alignment_dims$gene <- unlist(lapply(strsplit(alignment_dims$uid, "_"), `[[`, 2))
  # Add a column for whether the species have been shuffled
  alignment_dims$shuffled_taxa <- "FALSE"
  # Add columns noting the number and percent of changed taxa labels
  alignment_dims$number_unchanged_taxa <- alignment_dims$num_taxa
  alignment_dims$percent_unchanged_taxa <- 100
  # Reorder columns
  alignment_dims <- alignment_dims[, c("uid", "dataset", "clade", "gene", "num_taxa", "num_sites", "shuffled_taxa", "number_unchanged_taxa", "percent_unchanged_taxa",  "alignment_path")]
  # Sort so that the genes are in order of increasing number of taxa
  alignment_dims <- alignment_dims[order(as.numeric(alignment_dims$num_taxa), alignment_dims$gene), ]
  # Save the nice data frame of information
  write.csv(alignment_dims, alignment_dim_file, row.names = F)
}

## Shuffle the names in the alignments
# For each alignment, take a row from the alignment_dims csv file. Take that alignment, copy it, and shuffle the sites.
# Return a summary row about the new alignment



#### 4. Construct parameters dataframe for empirical alignments ####
# Create filename for Oaks 2011 parameters csv
oaks_csv_file <- paste0(results_directory, "Oaks2011_parameters.csv")

# Find or create the Oaks 2011 parameters csv
if (file.exists(oaks_csv_file)){
  oaks_df <- read.csv(oaks_csv_file)
} else {
  ## Find the genes for analysis
  # Collect all alignments from Oaks 2011 dataset
  all_files <- list.files(oaks_directory, recursive = TRUE)
  # Get all genes (will have the phrase "gene_alignments" in file name)
  all_genes <- grep("gene_alignments", all_files, value = TRUE)
  # Make list of the names of the genes you want to extract
  gene_names <- c("cmos", "CYTB", "ND2", "ND3")
  # Identify the files for those genes
  run_genes <- grep(paste(gene_names, collapse = "|"), all_genes, value = TRUE)
  # Extract the partitions of interest for those genes
  partition_genes <- grep("_12", run_genes, value = TRUE, invert = TRUE)
  # Extract the subsets of interest for those genes
  subset_genes <- c(grep("8taxa", partition_genes, value = TRUE), grep("23taxa", partition_genes, value = TRUE), grep("79taxa", partition_genes, value = TRUE))
  # Add location to get full file path to each gene 
  subset_gene_paths <- paste0(oaks_directory, subset_genes)
  
  ## Construct a dataframe containing information for analysis
  # Identify whether each gene is mtDNA or nDNA
  gene_type <- gsub("_gene_alignments", "", unlist(strsplit(subset_genes, "/"))[c(T,F)])
  # Get just the file name of each alignment file (include no part of the file path)
  gene_file_name <- basename(subset_genes)
  # Get the name of each gene
  gene_split <- strsplit(gene_file_name, "_")
  gene_name <- unlist(lapply(gene_split, `[[`, 1))
  # Identify the number of taxa in each gene
  gene_num_taxa <- as.numeric(gsub("taxa.fa", "", grep("taxa", unlist(strsplit(gene_file_name, "_")), value = TRUE)))
  # Identify the partitioning scheme for each gene
  codon_partitions <- gene_file_name
  codon_partitions[grep("_1_", gene_file_name)] <- "1"
  codon_partitions[grep("_2_", gene_file_name)] <- "2"
  codon_partitions[grep("_3_", gene_file_name)] <- "3"
  codon_partitions[grep("_1_|_2_|_3_", gene_file_name, invert = TRUE)] <- "All"
  # Assemble into a dataframe
  oaks_df <- data.frame(row_id = 1:length(subset_gene_paths),
                        uid = paste0(gene_name, "_", gene_type, "_", gene_num_taxa, "_", codon_partitions), 
                        gene_name = gene_name, 
                        num_taxa = gene_num_taxa, 
                        codon_position = codon_partitions, 
                        DNA_type = gene_type, 
                        old_alignment_file = subset_gene_paths)
  # Save dataframe
  write.csv(oaks_df, file = oaks_csv_file, row.names = FALSE)
}



#### 5. Copy alignments and save the parameters for each alignment ####
# Make a new folder for the Oaks 2011 analyses
copy_directory <-  paste0(results_directory, "Oaks2011/")
if (dir.exists(copy_directory) == FALSE){dir.create(copy_directory)}
# If the columns for output alignments and output csvs do not exist, add them
if (length(grep("output_alignment_file", names(oaks_df))) == 0){
  # Make a new column in the dataframe. This column is the location where each alignment will be copied and stored
  oaks_df$output_alignment_file <- paste0(copy_directory, oaks_df$uid, "/", oaks_df$uid, "_output_alignment.fa")
  oaks_df$parameters_path <- paste0(copy_directory, oaks_df$uid, "/", oaks_df$uid, "_parameters.csv")
  # Save Oaks 2011 dataframe with the two new columns
  write.csv(oaks_df, file = oaks_csv_file, row.names = FALSE)
}
# If the output alignment paths do not exist, create them by copying the alignments to their new home
if ((FALSE %in% file.exists(oaks_df$output_alignment_file)) == TRUE){
  # Feed each row into the function to copy the alignment and save the corresponding alignment parameters 
  lapply(1:nrow(oaks_df), copy.empirical.alignment, data_df = oaks_df)
}



#### 6. Apply tests for treelikeness to each empirical alignment ####
# Get list of all the Oaks 2011 alignments to run
all_oaks_alignments <- oaks_df$output_alignment_file
# Remove the alignments with either a .log or .iqtree file
oaks_alignments_to_run <- unique(c(all_oaks_alignments[!file.exists(paste0(all_oaks_alignments, ".log"))], 
                                   all_oaks_alignments[!file.exists(paste0(all_oaks_alignments, ".iqtree"))] ))
# Remove the alignments that are not partitioned by codon position
oaks_alignments_to_run <- grep("All", oaks_alignments_to_run, value = TRUE, invert = TRUE)
# Apply treelikeness metrics to all alignments 
if (length(oaks_alignments_to_run) > 0){
  mclapply(oaks_alignments_to_run, treelikeness.metrics.empirical,
           iqtree2_path, splitstree_path, phylogemetric_path, fast_TIGER_path,
           supply_number_of_taxa = FALSE, number_of_taxa = NA,
           num_iqtree2_threads = "AUTO", num_iqtree2_scf_quartets = 100,
           iqtree_substitution_model = "MFP", distance_matrix_substitution_method = "F81",
           num_phylogemetric_threads = NA, tree_proportion_remove_trivial_splits = TRUE,
           run_splitstree_for_tree_proportion = TRUE, sequence_format = "DNA",
           apply.TIGER = TRUE, redo = FALSE,
           mc.cores = num_cores)
}

# Collect and collate results
oaks_list <- mclapply(all_oaks_alignments, collate.empirical.treelikeness.results, mc.cores = num_cores)
# Remove NULL objects in list (indicates treelikeness metrics csv does not exist for this alignment)
keep_indexes <- which(!sapply(oaks_list, is.null))
oaks_list_filtered <- oaks_list[keep_indexes]
# Save output dataframe
oaks_tl_df <- as.data.frame(do.call("rbind", oaks_list_filtered))
oaks_df_name <- paste0(results_directory, "Oaks2011_treelikeness_metrics_collated_results.csv")
write.csv(oaks_tl_df, oaks_df_name, row.names = FALSE)
