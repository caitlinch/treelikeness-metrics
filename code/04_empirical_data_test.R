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



#### 3. Add gaps to the replicate alignments ####
if (control_parameters$edit.replicate.alignments == TRUE){
  ## Extract list of alignments
  all_files <- list.files(replicate_alignment_directory, recursive = T)
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Apply the gaps from the empirical alignment to the replicate alignments
  lapply(wea17_replicates, copy.alignment.gaps, template_alignment_path = wea17_al)
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Apply the gaps from the empirical alignment to the replicate alignments
  lapply(wea17f_replicates, copy.alignment.gaps, template_alignment_path = wea17f_al)
}



#### 4. Apply tests for treelikeness to each alignment ####
if (control_parameters$apply.treelikeness.tests == TRUE){
  ## Extract list of alignments
  all_files <- list.files(replicate_alignment_directory, recursive = T)
  
  ## For WEA17
  # Extract the list of alignments using the id (WEA17F)
  wea17_files <- paste0(replicate_alignment_directory, grep("WEA17_", all_files, value = T))
  wea17_al <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T, invert = T)]
  wea17_replicates <- wea17_files[grep("bs_rep", basename(wea17_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17_dir <- paste0(output_directory, "WEA17/")
  if (dir.exists(wea17_dir) == FALSE){ dir.create(wea17_dir) }
  wea17_copy <- paste0(wea17_dir, basename(wea17_al))
  file.copy(from = wea17_al, to = wea17_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  # Apply the wrapper function to calculate treelikeness of bootstrap replicates
  wea17_df <- bootstrap.wrapper(wea17_replicates, 
                                output_directory = output_directory, 
                                splitstree_path = splitstree_path, 
                                iqtree2_path = iqtree2_path, 
                                iqtree_model = "'Q.insect+R8{0.1639,0.0315,0.1812,0.1846,0.1454,0.4618,0.1171,0.7116,0.1742,1.1628,0.1403,2.0884,0.0676,3.6840,0.0103,6.4538}'", 
                                num_iqtree2_threads = "10", 
                                sequence_format = "AA", 
                                redo = FALSE, 
                                number_parallel_cores = mclapply_num_cores)
  
  ## For WEA17_filtered
  # Extract the list of alignments using the id (WEA17F)
  wea17f_files <- paste0(replicate_alignment_directory, grep("WEA17F_", all_files, value = T))
  wea17f_al <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T, invert = T)]
  wea17f_replicates <- wea17f_files[grep("bs_rep", basename(wea17f_files), ignore.case = T)]
  # Copy original alignment and prepare for treelikeness metric runs
  wea17f_dir <- paste0(output_directory, "WEA17F/")
  if (dir.exists(wea17f_dir) == FALSE){ dir.create(wea17f_dir) }
  wea17f_copy <- paste0(wea17f_dir, basename(wea17f_al))
  file.copy(from = wea17f_al, to = wea17f_copy)
  # Calculate the tree proportion for the original alignments
  treelikeness.metrics.empirical(wea17f_copy,
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = num_cores, 
                                 sequence_format = "AA", 
                                 redo = FALSE)
  
  # Apply the wrapper function
  wea17f_df <- bootstrap.wrapper(wea17f_replicates, 
                                 output_directory = output_directory, 
                                 splitstree_path = splitstree_path, 
                                 iqtree2_path = iqtree2_path, 
                                 iqtree_model = "'Q.insect+I{0.0659}+R6{0.2004,0.1133,0.1734,0.3827,0.1439,0.6654,0.2108,1.1677,0.1549,2.2361,0.0508,4.3870}'", 
                                 num_iqtree2_threads = "10", 
                                 sequence_format = "AA", 
                                 redo = FALSE, 
                                 number_parallel_cores = mclapply_num_cores)
}



#### 5. Calculate p-values for each test statistic ####
if (control_parameters$calculate.empirical.p_values == TRUE){
  ## Find file with parametric bootstrap results
  all_output <- list.files(output_directory, recursive = TRUE)
  all_csv <- grep("csv", all_output, value = T)
  bs_csv_file <- paste0(output_directory, grep("empirical", grep("collated", grep("treelikeness_metrics", all_csv, value = T), value = T), value = T))
  # Open csv file
  bs_df <- read.csv(bs_csv_file, stringsAsFactors = F)
  # Add identification columns for faceting and filtering
  bs_df$rep_type <- factor(bs_df$unique_id,
                           levels = bs_df$unique_id,
                           labels = rep(c(rep("Bootstrap replicate", 100), "Alignment"), 2),
                           ordered = TRUE)
  # Add a column specifying which alignment the row is
  bs_df$alignment_id <- factor(bs_df$unique_id,
                               levels = bs_df$unique_id,
                               labels = rep(c("WEA17", "WEA17F"), each = 101),
                               ordered = TRUE)
  
  ## Calculate p-values for "WEA17" alignment
  wea17_df <- bs_df[grep("WEA17F", bs_df$unique_id, invert = T), ]
  wea17_df$unique_id[grep("rep", wea17_df$unique_id, invert = T)] <- "alignment"
  wea17_p_value_df <- calculate.all.p_values(output_df = wea17_df, 
                                             test_statistic_names = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  wea17_p_value_df$dataset <- "WEA17"
  wea17_p_value_df$statistic_type  <- "p-value"
  wea17_p_value_df <- wea17_p_value_df[c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "p_value")]
  names(wea17_p_value_df) <- c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "statistic_value")
  # Extract Network Treelikeness Test values
  wea17_ntlt <- bs_df[which(bs_df$rep_type == "Alignment" & bs_df$alignment_id == "WEA17"),]$NetworkTreelikenessTest
  wea17_ntlt_stat_value <- length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17" & bs_df$NetworkTreelikenessTest == "Treelike"))/
    length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17"))
  wea17_ntlt_vector <- c("dataset" = "WEA17", "test_statistic" = "NetworkTreelikenessTest", 
                         "test_statistic_value" = wea17_ntlt, "statistic_type" = "proportion_treelike_alignments", 
                         "statistic_value" = wea17_ntlt_stat_value) 
  # Bind the tree proportion results
  wea17_p_value_df <- rbind(wea17_p_value_df, wea17_ntlt_vector)
  
  ## Calculate p-values for "WEA17F" alignment
  wea17f_df <- bs_df[grep("WEA17F", bs_df$unique_id), ]
  wea17f_df$unique_id[grep("rep", wea17f_df$unique_id, invert = T)] <- "alignment"
  wea17f_p_value_df <- calculate.all.p_values(output_df = wea17f_df, 
                                              test_statistic_names = c("LM_proportion_resolved_quartets", "sCF_mean", "mean_delta_plot_value", "Cunningham_test", "tree_proportion"))
  wea17f_p_value_df$dataset <- "WEA17F"
  wea17f_p_value_df$statistic_type  <- "p-value"
  wea17f_p_value_df <- wea17f_p_value_df[c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "p_value")]
  names(wea17f_p_value_df) <- c("dataset", "test_statistic", "test_statistic_value", "statistic_type", "statistic_value")
  # Extract Network Treelikeness Test values
  wea17f_ntlt <- bs_df[which(bs_df$rep_type == "Alignment" & bs_df$alignment_id == "WEA17F"),]$NetworkTreelikenessTest
  wea17f_ntlt_stat_value <- length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F" & bs_df$NetworkTreelikenessTest == "Treelike"))/
    length(which(bs_df$rep_type == "Bootstrap replicate" & bs_df$alignment_id == "WEA17F"))
  wea17f_ntlt_vector <- c("dataset" = "WEA17F", "test_statistic" = "NetworkTreelikenessTest", 
                          "test_statistic_value" = wea17f_ntlt, "statistic_type" = "proportion_treelike_alignments", 
                          "statistic_value" = wea17f_ntlt_stat_value) 
  # Bind the tree proportion results
  wea17f_p_value_df <- rbind(wea17f_p_value_df, wea17f_ntlt_vector)
  
  ## Combine and save dataframes
  collated_p_value_df <- rbind(wea17_p_value_df, wea17f_p_value_df)
  collated_p_value_file <- paste0(output_directory, "empirical_collated_p_values.csv")
  write.csv(collated_p_value_df, file = collated_p_value_file, row.names = F)
}



#### 6. Apply tests for treelikeness to each gene ####
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
  
  ## Apply metrics
  gene_metrics <- lapply(1:nrow(gene_df), treelikeness.metrics.without.bootstrap, 
                         df = gene_df, tl_output_directory = "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/05_empirical_treelikeness_results/genes_treelikeness_metrics/", 
                         splitstree_path = splitstree_path, iqtree2_path = iqtree2_path, 
                         num_iqtree2_threads = "3", sequence_format = "AA", 
                         redo = FALSE)
}


