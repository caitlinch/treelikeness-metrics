# caitlinch/treelikeness-metrics/code/01_simulations.R
# Caitlin Cherryh 2023

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above) and ms.



#### 1. Set parameters ####
## CONTROL PARAMETERS
# parameter.values                <- Control flag for creating the simulation parameter values (TRUE to create objects containing simulation parameter values)
# run.experiment.1                <- Control flag for experiment 1 (TRUE to run code to generate simulations for Experiment 1)
# run.experiment.3                <- Control flag for experiment 3 (TRUE to run code to generate simulations for Experiment 3)

## DIRECTORY PATHS
# simulation_directory            <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory                  <- Location of caitlinch/treelikeness-metrics github repository (for access to functions).
# ms_path                         <- Path to ms executable 
# iqtree2_path                    <- Path to IQ-Tree2 executable (version 2.2-beta or later to ensure Alisim is included). 
# number_parallel_threads         <- Number of threads to run simultaneously in mclapply when generating alignments

## SIMULATION PARAMETERS
# total_alignment_length          <- Total length of concatenated alignments in base pairs (we chose 10000) for the random tree analyses (experiment 1).
# gene_length                     <- Length of each gene generated in ms. Total alignment length will be length of each gene multiplied by number of gene trees. 
#                                     For Total Alignment Length = 10000, use 100 gene trees of 100 bp each.
# sequence_type                   <- Sequence type for simulation (we chose "DNA").
# taxa_vec                        <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                        <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# tree_depth_random_sims          <- One or more values for the tree length of randomly generated trees (in substitutions per site).
# tree_depth_coalescent_sims      <- One or more values for the tree length of coalescent trees (in coalescent units).
# speciation_rates                <- One of more values for the speciation rate for the introgression simulations (experiment 3)
# number_gene_trees               <- Number of gene trees to generate for coalescent simulations
# r_vec                           <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models              <- Model of sequence evolution for Alisim 
# alisim_gene_tree_length         <- Gene-specific tree length for Alisim
# conversion_depth_subs_per_site     <- Tree depth for ML tree - used to convert depth of ms gene trees from coalescent units to substitutions per site


## CONTROL PARAMETERS
parameter.values            <- FALSE
run.experiment.1            <- FALSE
run.experiment.3            <- TRUE

## DIRECTORY PATHS
run_location = "soma"
if (run_location == "local"){
  simulation_directory    <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  repo_directory          <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"
  ms_path                 <- "ms"
  iqtree2_path            <- "iqtree2"
  number_parallel_threads <- 1
} else if (run_location == "soma"){
  simulation_directory    <- "/data/caitlin/treelikeness_metrics/"
  repo_directory          <- "/data/caitlin/treelikeness_metrics/"
  ms_path                 <- "/data/caitlin/executables/msdir/ms"
  iqtree2_path            <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  number_parallel_threads <- 20
}

## SIMULATION PARAMETERS
if (parameter.values == TRUE){
  total_alignment_length          <- 10000
  gene_length                     <- 500
  sequence_type                   <- "DNA"
  taxa_vec                        <- c(5,10,20,50,100)
  num_reps                        <- 10
  tree_depth_random_sims          <- c(0.01, 0.1, 1)
  tree_depth_coalescent_sims      <- c(5, 50, 500) # where bounds for coalescent tree depth are in millions of years (see TreeSim doco)
  speciation_rates                <- c(0.1, 1) # for generating Yule tree for the introgression experiments (experiment 3)
  number_gene_trees               <- 200
  r_vec                           <- seq(0, 0.5, 0.05)
  alisim_gene_models              <- "JC"
  alisim_gene_tree_length         <- NA
  conversion_depth_subs_per_site  <- 0.1
  
  # Determine number of trees - all whole numbers that are a divisor of the total alignment length
  number_of_trees <- divisors(total_alignment_length)
  # Set number of taxa equal to taxa_vec
  number_of_taxa <- taxa_vec
  # Create a list of all replicate numbers using the num_reps value
  number_of_replicates <- 1:num_reps
}



#### 2. Open packages ####
library(ape)
library(phytools)
library(parallel)



#### 3. Prepare analyses ####
source(paste0(repo_directory, "code/func_simulating_alignments.R"))



#### 4. Generate simulations ####
## Experiment 1: Random trees ##
# Generate x random trees with y taxa 
#     Total alignment length = 10,000 bp
#     Number of trees ranges from 0 to 10,000 (whole number divisors of 10000)
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 1000
# Simulate DNA along each tree with Alisim, using the topology-unlinked partition model
if (run.experiment.1 == TRUE){
  # Create folder to store results of this experiment, if it doesn't already exist
  exp1_dir <- paste0(simulation_directory, "exp_1/")
  if(!file.exists(exp1_dir)){dir.create(exp1_dir)}
  # Create file path for parameters csv
  exp1_df_path <- paste0(simulation_directory, "exp1_parameters.csv")
  
  if (file.exists(exp1_df_path) == TRUE){
    exp1_params <- read.csv(exp1_df_path)
  } else {
    # Create matrix with parameters for generating each simulated alignment
    exp1_params <- expand.grid("num_reps" = number_of_replicates, "num_taxa" = number_of_taxa, "num_trees" = number_of_trees, "tree_depth" = tree_depth_random_sims)
    # Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`_`tree_depth`
    exp1_params$uid <- paste0("exp1_",sprintf("%05d", exp1_params$num_trees), "_", sprintf("%04d", exp1_params$num_taxa), "_",
                              sprintf("%03d", exp1_params$num_reps), "_", exp1_params$tree_depth)
    # Add parameters for Alisim
    exp1_params$alisim_gene_models <- alisim_gene_models
    exp1_params$alisim_gene_tree_length <- alisim_gene_tree_length
    # Add other parameters
    exp1_params$total_alignment_length <- total_alignment_length
    exp1_params$sequence_type <- sequence_type
    # Add names for the tree file, partition file and output alignment file for each simulated alignment
    exp1_params$tree_file <- paste0(exp1_params$uid, "_random_trees.phy")
    exp1_params$partition_file <- paste0(exp1_params$uid, "_partitions.nex")
    exp1_params$output_alignment_file <- paste0(exp1_params$uid, "_output_alignment")
    
    # Write exp1_params dataframe to file as a csv
    write.csv(exp1_params, file = exp1_df_path, row.names = TRUE)
  }
  
  # Iterate through each row in the parameters dataframe
  # Run all reps:
  #   lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
  # Run single rep:
  #   lapply(1, random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
  
  if (number_parallel_threads == 1){
    exp1_op_list <- lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
  } else {
    exp1_op_list <- mclapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params,
                             mc.cores = number_parallel_threads)
  }
  
  # Change output file names from list to dataframe
  exp1_op_df <- as.data.frame(do.call(rbind, exp1_op_list))
  exp1_op_df_path <- paste0(simulation_directory, "exp1_file_output_paths.csv")
  write.csv(exp1_op_df, file = exp1_op_df_path, row.names = TRUE)
}



## Experiment 3: Mimicking ILS and introgression with a Yule tree##
# Generate gene trees using ms to simulate ILS with a single introgression event
#   - Fix alignment at total_alignment_length
#   - Simulate gene trees for tree containing a single introgression event in ms. Vary proportion of introgressed DNA from 0 to 1
#   - Simulate DNA along each tree with Alisim
#   - Concatenate alignments
if (run.experiment.3 == TRUE){
  # Create folder to store results of this experiment, if it doesn't already exist
  exp3_dir <- paste0(simulation_directory, "exp_3/")
  if (dir.exists(exp3_dir) == FALSE){dir.create(exp3_dir)}
  
  # Create file path for parameters csv
  exp3_df_path <- paste0(simulation_directory, "exp3_parameters.csv")
  
  if (file.exists(exp3_df_path) == TRUE){
    exp3_params <- read.csv(exp3_df_path)
  } else {
    exp3_params <- expand.grid("num_reps" = number_of_replicates, "num_taxa" = number_of_taxa, "num_trees" = number_gene_trees, 
                               "tree_depth_coalescent" = tree_depth_coalescent_sims, "recombination_value" = r_vec, "recombination_type" = c("Ancient","Recent"),
                               "speciation_rate" = speciation_rates)
    # Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`_`tree_depth`_`recombination proportion`_`introgression event type`
    exp3_params$uid <- paste0("exp3_",sprintf("%05d", exp3_params$num_trees), "_", sprintf("%04d", exp3_params$num_taxa), "_",
                              sprintf("%03d", exp3_params$num_reps), "_", exp3_params$tree_depth, "_", exp3_params$recombination_value,
                              "_", exp3_params$recombination_type, "_", exp3_params$speciation_rate)
    # Add parameters for Alisim
    exp3_params$alisim_gene_models <- alisim_gene_models
    exp3_params$alisim_gene_tree_length <- alisim_gene_tree_length
    # Add other parameters
    exp3_params$tree_depth_subs_per_sites <- conversion_depth_subs_per_site
    exp3_params$total_alignment_length <- number_gene_trees * gene_length
    exp3_params$sequence_type <- sequence_type
    # Add name for the partition file and output alignment file for each simulated alignment
    exp3_params$partition_file <- paste0(exp3_params$uid, "_partitions.nex")
    exp3_params$output_alignment_file <- paste0(exp3_params$uid, "_output_alignment")
    
    # Remove any rows that have 5 taxa and an ancient introgression event
    # Due to the way ancient introgression events are structured (take place at time where four taxa exist, between two non-sister taxa), a recent and an ancient
    #     introgression event will be identical for a tree with 5 taxa
    remove_rows <- which(exp3_params$num_taxa == 5 & exp3_params$recombination_type == "Ancient")
    keep_rows <- setdiff(1:nrow(exp3_params), remove_rows)
    exp3_params <- exp3_params[keep_rows, ]
    row.names(exp3_params) <- 1:nrow(exp3_params)
    
    # Write exp3_params dataframe to file as a csv
    write.csv(exp3_params, file = exp3_df_path, row.names = TRUE)
  }
  
  # Errors testing: 900 (Ancient) and 6666 (Recent)
  
  # Iterate through each row in the parameters dataframe and generate an alignment for each set of parameters
  # Run all reps: 
  #   lapply(1:nrow(exp3_params), ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp3_params)
  # Run single rep:
  #   lapply(1, ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp3_params, select.sister = FALSE)
  if (number_parallel_threads == 1){
    exp3_op_list<- lapply(1:nrow(exp3_params), ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path,
                          experiment_params_df = exp3_params, select.sister = FALSE, scale.gene.trees = TRUE)
  } else {
    exp3_op_list <- mclapply(1:nrow(exp3_params), ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path,
                             experiment_params_df = exp3_params, select.sister = FALSE, scale.gene.trees = TRUE, 
                             mc.cores = number_parallel_threads)
  }
  
  # Change output file names from list to dataframe
  exp3_op_df <- as.data.frame(do.call(rbind, exp3_op_list))
  exp3_op_df_path <- paste0(simulation_directory, "exp3_file_output_paths.csv")
  write.csv(exp3_op_df, file = exp3_op_df_path, row.names = TRUE)
}


