# caitlinch/treelikeness-metrics/code/04_empirical_data_test.R
# Caitlin Cherryh 2023

## Script summary:
# This program will apply various tests for treelikeness to genes from empirical alignments
# This program requires IQ-Tree2 (2.2-beta or above), fast TIGER, phylogemetric, and SplitsTree (4.17.2 or above).

#### 1. Set parameters ####
## Directories
# output_directory                <- Directory where alignments will be saved/treelikeness metrics will be run
# replicate_alignment_directory   <- Location of bootstrap replicate alignments, generated using Alisim in IQ-Tree
# repo_directory                  <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

## Executable paths
# iqtree2_path            <- Path to IQ-Tree (requires version 2.2-beta or above)
# splitstree_path         <- Path to SplitsTree 4 version 4.17.2 or above

## Run parameters
# num_cores               <- Number of parallel threads to use at once

run_location = "local"
if (run_location == "local"){
  # Directories
  empirical_data_directory        <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/03_empirical_tree_estimation/"
  output_directory                <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/05_empirical_treelikeness_results/"
  replicate_alignment_directory   <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/04_bs_replicate_alignments/"
  repo_directory                  <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  
  # Executable paths
  iqtree2_path      <- "iqtree2"
  splitstree_path   <- "/Applications/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
  
  # Run parameters
  num_cores <- 1
} else if (run_location == "soma"){
  # Directories
  output_directory                <- "/data/caitlin/treelikeness_metrics/empirical_treelikeness_output/"
  replicate_alignment_directory   <- "/data/caitlin/treelikeness_metrics/empirical_alignments/"
  repo_directory                  <- "/data/caitlin/treelikeness_metrics/"
  
  # Executable paths
  iqtree2_path      <- "/data/caitlin/executables/iqtree-2.2.2-Linux/bin/iqtree2"
  splitstree_path   <- "/home/caitlin/splitstree4/SplitsTree"
  
  # Run parameters
  num_cores <- 30
}

control_parameters <- list(edit.replicate.alignments = FALSE,
                           apply.treelikeness.tests = FALSE,
                           calculate.empirical.p_values = FALSE,
                           per.gene.analysis = TRUE)



#### 2. Prepare analyses ####
# Open packages
library(parallel)

# Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "code/func_empirical.R"))
source(paste0(repo_directory, "code/func_parametric_bootstrap.R"))
source(paste0(repo_directory, "code/func_metrics.R"))
source(paste0(repo_directory, "code/func_data_analysis.R"))

if (control_parameters$apply.treelikeness.tests == TRUE){
  # Calculate the number of parallel processes to run using mc.apply
  mclapply_num_cores <- num_cores/10
  if (mclapply_num_cores < 1){
    # Minimum of one process running at once
    mclapply_num_cores <- 1
  } else {
    # Take floor of the decimal to run a conservative number of processes 
    #     (e.g. if 31 cores and each IQ-Tree run uses 10, you can run 3 alignments at once)
    mclapply_num_cores <- floor(mclapply_num_cores)
  }
}



#### 3. Apply tests for treelikeness to each gene in the 2 empirical datasets ####
if (control_parameters$per.gene.analysis == TRUE){
  ## Split each alignment into individual genes
  # List all files in the empirical data directory
  data_files <- paste0(empirical_data_directory, list.files(empirical_data_directory, recursive = TRUE))
  data_files <- grep("genes", data_files, value = T)
  # Split Whelan 2017 alignment into separate genes
  wea17_alignment_file <- grep("filtered", grep("WEA17.fa", data_files, value = T), value = T, invert = TRUE)
  wea17_partition_file <- grep("filtered", grep("partitions", grep("Whelan2017", data_files, value = T), value = T), value = T, invert = T)
  wea17_gene_directory <- paste0(dirname(wea17_alignment_file), "/")
  wea17_gene_paths <- split.partitions(alignment_file = wea17_alignment_file, partition_file = wea17_partition_file, gene_output_directory = wea17_gene_directory)
  # Split Whelan 2017 filtered alignment into separate genes
  wea17f_alignment_file <- grep("fa", grep("WEA17F", data_files, value = T), value = T)
  wea17f_partition_file <- grep("partitions", grep("Whelan2017_filtered", data_files, value = T), value = T)
  wea17f_gene_directory <- paste0(dirname(wea17f_alignment_file), "/")
  wea17f_gene_paths <- split.partitions(alignment_file = wea17f_alignment_file, partition_file = wea17f_partition_file, gene_output_directory = wea17f_gene_directory)
  
  ## Assemble dataframe of genes
  gene_df <- data.frame(dataset = c(rep("WEA17", length(wea17_gene_paths)), rep("WEA17F", length(wea17f_gene_paths))), 
                        gene = c(gsub(".fa", "", basename(wea17_gene_paths)), gsub(".fa", "", basename(wea17f_gene_paths))),
                        alignment_path = c(wea17_gene_paths, wea17f_gene_paths))
  gene_df$ID <- paste0(gene_df$dataset, ".", gene_df$gene)
  
  ## For each gene: estimate tree, then apply metrics: tree proportion, sCF, and delta plots
  # Prepare output directory
  gene_output_directory <- paste0(output_directory, "genes/")
  if (dir.exists(gene_output_directory) == FALSE){dir.create(gene_output_directory)}
  # Apply metrics
  gene_metrics <- lapply(2:nrow(gene_df), treelikeness.metrics.with.parametric.bootstrap, 
                         df = gene_df, tl_output_directory = gene_output_directory, 
                         splitstree_path = splitstree_path, iqtree2_path = iqtree2_path, 
                         num_iqtree2_threads = "3", sequence_format = "AA", 
                         redo = FALSE, number_parallel_cores = 1)
  
  ## Collate and output the p-values
  pvalue_df <- as.data.frame(do.call(rbind, gene_metrics))
  write.csv(pvalue_df, file = paste0(repo_directory, "gene_treelikeness_pvalues.csv"), row.names = F)
}


