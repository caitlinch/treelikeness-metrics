## Supplementary analysis A: Network Treelikeness Test

# caitlinch/treelikeness-metrics/supp_code/suppA_01_simulations
# Caitlin Cherryh 2025

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above) and ms.


#### 01. Set parameters ####
## CONTROL PARAMETERS
# parameter.values  <- Control flag for creating the simulation parameter values (TRUE to create objects containing simulation parameter values)
# run.experiment.A1 <- Control flag for experiment A1 (TRUE to run code to generate simulations for Experiment A1)
# run.experiment.A2 <- Control flag for experiment A2 (TRUE to run code to generate simulations for Experiment A2)

## DIRECTORY PATHS
# simulation_directory    <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions).
# ms_path                 <- Path to ms executable
# iqtree2_path            <- Path to IQ-Tree2 executable (version 2.2-beta or later to ensure Alisim is included).
# number_parallel_threads <- Number of threads to run simultaneously in mclapply when generating alignments

## SIMULATION PARAMETERS
# total_alignment_length      <- Total length of concatenated alignments in base pairs (we chose 10000) for the random tree analyses (experiment 1).
# gene_length                 <- Length of each gene generated in ms. Total alignment length will be length of each gene multiplied by number of gene trees.
#                                 For Total Alignment Length = 10000, use 100 gene trees of 100 bp each.
# sequence_type               <- Sequence type for simulation (we chose "DNA").
# taxa_vec                    <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                    <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# tree_depth_random_sims      <- One or more values for the tree length of randomly generated trees (in substitutions per site).
# tree_depth_coalescent_sims  <- One or more values for the tree length of coalescent trees (in coalescent units).
# speciation_rates            <- One of more values for the speciation rate for the introgression simulations (experiment 3)
# number_gene_trees           <- Number of gene trees to generate for coalescent simulations
# r_vec                       <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models          <- Model of sequence evolution for Alisim
# alisim_gene_tree_length     <- Gene-specific tree length for Alisim


## CONTROL PARAMETERS
parameter.values    <- FALSE
run.experiment.A1   <- FALSE
run.experiment.A2   <- TRUE

## DIRECTORY PATHS
run_location = "dayhoff"
if (run_location == "dayhoff"){
  repo_directory          <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/"
  simulation_directory    <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/suppA/"
  iqtree2_path            <- "/mnt/data/dayhoff/home/u5348329/treelikeness_supp/software/iqtree-2.2.2.6-Linux/bin/iqtree2"
  number_parallel_threads <- 10
} else if (run_location == "WSL"){
  repo_directory          <- ""
  simulation_directory    <- "suppA/"
  iqtree2_path             <- "iqtree2"
  number_parallel_threads <- 1
  if (dir.exists(simulation_directory) == FALSE){dir.create(simulation_directory)}
}

## SIMULATION PARAMETERS
if (parameter.values == TRUE){
  # Experiment SA1
  SA1_alignment_length            <- 10000
  SA1_num_taxa                    <- seq(10, 100, 10)
  SA1_num_trees                   <- 1
  SA1_sequence_type               <- "DNA"
  SA1_tree_depth                  <- c(0.01, 0.1, 1)
  SA1_reps                        <- 1:100
  SA1_alisim_gene_models          <- "JC"
  SA1_alisim_gene_tree_length     <- NA
  # Experiment SA2
  SA2_num_taxa                    <- 100
  SA2_num_trees                   <- 1:10
  SA2_sequence_type               <- "DNA"
  SA2_tree_depth                  <- 1.0
  SA2_reps                        <- 1:100
  SA2_alisim_gene_models          <- "JC"
  SA2_alisim_gene_tree_length     <- NA
  SA2_alignment_length            <- 10000 * 1:10
}

#### 02. Prepare packages and functions ####
# Load packages
library(ape)
library(phytools)
library(parallel)

# Load functions to simulate alignments
source(paste0(repo_directory, "supp_code/suppA_func_random_tree_simulations.R"))


#### 4. Generate simulations ####
## Experiment A1: Random trees ##
# Generate x random trees with y taxa
#     Total alignment length = 10,000 bp
#     Number of trees ranges from 0 to 10,000 (whole number divisors of 10000)
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 100 (intervals of 10)
# Simulate DNA along each tree with Alisim, using the topology-unlinked partition model
if (run.experiment.A1 == TRUE){
  # Create folder to store results of this experiment, if it doesn't already exist
  expA1_dir <- paste0(simulation_directory, "exp_A1/")
  if(!file.exists(expA1_dir)){dir.create(expA1_dir)}
  # Create file path for parameters csv
  expA1_df_path <- paste0(simulation_directory, "expA1_parameters.csv")
  # Open or create parameters dataframe
  if (file.exists(expA1_df_path) == TRUE){
    expA1_params <- read.csv(expA1_df_path)
  } else {
    # Create matrix with parameters for generating each simulated alignment
    expA1_params <- expand.grid(
      "num_reps" = SA1_reps,
      "num_taxa" = SA1_num_taxa,
      "num_trees" = SA1_num_trees,
      "tree_depth" = SA1_tree_depth,
      "total_alignment_length"= SA1_alignment_length
    )
    # Add a unique identifier (uid):
    # expA1_{num. trees}_{num. taxa}_{rep. num.}_{tree_depth}_{alnmt. length}
    expA1_params$uid <- paste0(
      "expA1_",
      sprintf("%05d", expA1_params$num_trees),
      "_",
      sprintf("%04d", expA1_params$num_taxa),
      "_",
      sprintf("%03d", expA1_params$num_reps),
      "_",
      expA1_params$tree_depth,
      "_",
      paste0(expA1_params$total_alignment_length/1000, "kbp")
    )
    # Add parameters for Alisim
    expA1_params$alisim_gene_models <- SA1_alisim_gene_models
    expA1_params$alisim_gene_tree_length <- SA1_alisim_gene_tree_length
    # Add other parameters
    expA1_params$sequence_type <- SA1_sequence_type
    # Add names for the tree file, partition file and output alignment file for each simulated alignment
    expA1_params$tree_file <- paste0(expA1_params$uid, "_random_trees.phy")
    expA1_params$partition_file <- paste0(expA1_params$uid, "_partitions.nex")
    expA1_params$output_alignment_file <- paste0(expA1_params$uid, "_output_alignment")
    # Write expA1_params dataframe to file as a csv
    write.csv(expA1_params, file = expA1_df_path, row.names = TRUE)
  }
  # Iterate through each row in the parameters dataframe
  if (number_parallel_threads == 1) {
    expA1_op_list <- lapply(
      1:nrow(expA1_params),
      random.trees.generate.alignment,
      output_directory = expA1_dir,
      iqtree2_path = iqtree2_path,
      experiment_params = expA1_params
    )
  } else {
    expA1_op_list <- mclapply(
      1:nrow(expA1_params),
      random.trees.generate.alignment,
      output_directory = expA1_dir,
      iqtree2_path = iqtree2_path,
      experiment_params = expA1_params,
      mc.cores = number_parallel_threads
    )
  }
  # Change output file names from list to dataframe
  expA1_op_df <- as.data.frame(do.call(rbind, expA1_op_list))
  expA1_op_df_path <- paste0(simulation_directory, "expA1_file_output_paths.csv")
  write.csv(expA1_op_df, file = expA1_op_df_path, row.names = TRUE)
}


## Experiment A2: Random trees ##
# Generate x random trees with y taxa
#     Total alignment length = N*10000bp, where N is the number of taxa
#     Number of trees ranges from 0 to 10,000 (whole number divisors of 10000)
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 100 (intervals of 10)
# Simulate DNA along each tree with Alisim, using the topology-unlinked partition model
if (run.experiment.A2 == TRUE){
  # Create folder to store results of this experiment, if it doesn't already exist
  expA2_dir <- paste0(simulation_directory, "exp_A2/")
  if(!file.exists(expA2_dir)){dir.create(expA2_dir)}
  # Create file path for parameters csv
  expA2_df_path <- paste0(simulation_directory, "expA2_parameters.csv")
  # Open or create parameters dataframe
  if (file.exists(expA2_df_path) == TRUE){
    expA2_params <- read.csv(expA2_df_path)
  } else {
    # Create matrix with parameters for generating each simulated alignment
    expA2_params <- expand.grid(
      "num_reps" = SA2_reps,
      "num_taxa" = SA2_num_taxa,
      "num_trees" = SA2_num_trees,
      "tree_depth" = SA2_tree_depth,
      "total_alignment_length" = SA2_alignment_length
    )
    # Add a unique identifier (uid):
    # expA1_{num. trees}_{num. taxa}_{rep. num.}_{tree_depth}_{alnmt. length}
    expA2_params$uid <- paste0(
      "expA2_",
      sprintf("%05d", expA2_params$num_trees),
      "_",
      sprintf("%04d", expA2_params$num_taxa),
      "_",
      sprintf("%03d", expA2_params$num_reps),
      "_",
      expA2_params$tree_depth,
      "_",
      paste0(expA2_params$total_alignment_length/1000, "kbp")
    )
    # Add parameters for Alisim
    expA2_params$alisim_gene_models <- SA2_alisim_gene_models
    expA2_params$alisim_gene_tree_length <- SA2_alisim_gene_tree_length
    # Add other parameters
    expA2_params$sequence_type <- SA2_sequence_type
    # Add names for the tree file, partition file and output alignment file for each simulated alignment
    expA2_params$tree_file <- paste0(expA2_params$uid, "_random_trees.phy")
    expA2_params$partition_file <- paste0(expA2_params$uid, "_partitions.nex")
    expA2_params$output_alignment_file <- paste0(expA2_params$uid, "_output_alignment")
    # Write expA2_params dataframe to file as a csv
    write.csv(expA2_params, file = expA2_df_path, row.names = TRUE)
  }
  # Iterate through each row in the parameters dataframe
  if (number_parallel_threads == 1) {
    expA2_op_list <- lapply(
      1:nrow(expA2_params),
      random.trees.generate.alignment,
      output_directory = expA2_dir,
      iqtree2_path = iqtree2_path,
      experiment_params = expA2_params
    )
  } else {
    expA2_op_list <- mclapply(
      1:nrow(expA2_params),
      random.trees.generate.alignment,
      output_directory = expA2_dir,
      iqtree2_path = iqtree2_path,
      experiment_params = expA2_params,
      mc.cores = number_parallel_threads
    )
  }
  # Change output file names from list to dataframe
  expA2_op_df <- as.data.frame(do.call(rbind, expA2_op_list))
  expA2_op_df_path <- paste0(simulation_directory, "expA2_file_output_paths.csv")
  write.csv(expA2_op_df, file = expA2_op_df_path, row.names = TRUE)
}


