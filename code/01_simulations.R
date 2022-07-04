# /caitlinch/treelikeness_metrics/01_simulations.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above) and ms.

# If running on soma, remember to activate the empirical_treelikeness conda environment first:
#     source /home/rob/anaconda3/etc/profile.d/conda.sh
#     conda activate empirical_treelikeness



#### 1. Open packages ####
library(ape)
library(phytools)
library(parallel)



#### 2. Set parameters ####
# local_directory                 <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory                  <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# ms_path                         <- Path to ms executable 
# iqtree2_path                    <- Path to IQ-Tree2 executable (version 2.2-beta or later to ensure Alisim is included). 
# number_parallel_threads         <- Number of threads to run simultaneously in mclapply when generating alignments
# total_alignment_length          <- Total length of concatenated alignments in base pairs (we chose 10000) for the random tree analyses.
# gene_length                     <- Length of each gene generated in ms. Total alignment length will be length of each gene multiplied by number of gene trees. 
#                                     For Total Alignment Length = 10000, use 100 gene trees of 100 bp each.
# sequence_type                   <- Sequence type for simulation (we chose "DNA").
# taxa_vec                        <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                        <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# tree_depth_random               <- One or more values for the tree length of randomly generated trees (in substitutions per site).
# tree_depth_coalescent           <- One or more values for the tree length of coalescent trees (in coalescent units).
# number_gene_trees               <- Number of gene trees to generate for coalescent simulations
# r_vec                           <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models              <- Model of sequence evolution for Alisim 
# alisim_gene_tree_length         <- Gene-specific tree length for Alisim

run_location = "soma"
if (run_location == "local"){
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
  ms_path <- "ms"
  iqtree2_path <- "iqtree2.2-beta"
  number_parallel_threads <- 1
} else if (run_location == "soma"){
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  repo_directory <- "/data/caitlin/treelikeness_metrics/code/"
  ms_path <- "/data/caitlin/executables/msdir/ms"
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.2.0-Linux/bin/iqtree2"
  number_parallel_threads <- 20
}

total_alignment_length <- 10000
gene_length <- 1000
sequence_type <- "DNA"
taxa_vec <- c(5,10,20,50,100)
num_reps <- 10
tree_depth_random <- c(0.01, 0.1, 1)
tree_depth_coalescent <- c(0.1, 1, 10, 100) # where bounds for coalescent tree depth are 0.1 (minimum) and 100 (maximum) in coalescent units
number_gene_trees <- 1000
r_vec <- seq(0, 0.5, 0.05)
alisim_gene_models <- "JC"
alisim_gene_tree_length <- NA



#### 3. Source functions from caitlinch/treelikeness_metrics and prepare variables ####
source(paste0(repo_directory, "func_simulating_alignments.R"))

## Prepare variables that remain stable for all experiments using the input parameters
# Determine number of trees - all whole numbers that are a divisor of the total alignment length
number_of_trees <- divisors(total_alignment_length)
# Set number of taxa equal to taxa_vec
number_of_taxa <- taxa_vec
# Create a list of all replicate numbers using the num_reps value
number_of_replicates <- 1:num_reps



#### 4. Generate simulations ####
## Experiment 1: Random trees ##
# Generate x random trees with y taxa 
#     Total alignment length = 10,000 bp
#     Number of trees ranges from 0 to 10,000 (whole number divisors of 10000)
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 1000
# Simulate DNA along each tree with Alisim, using the topology-unlinked partition model

# Create folder to store results of this experiment, if it doesn't already exist
exp1_dir <- paste0(local_directory, "exp_1/")
if(!file.exists(exp1_dir)){dir.create(exp1_dir)}

# Create matrix with parameters for generating each simulated alignment
exp1_params <- expand.grid("num_reps" = number_of_replicates, "num_taxa" = number_of_taxa, "num_trees" = number_of_trees, "tree_depth" = tree_depth_random)
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
exp1_df_path <- paste0(local_directory, "exp1_parameters.csv")
write.csv(exp1_params, file = exp1_df_path, row.names = TRUE)

# Iterate through each row in the parameters dataframe
# Run all reps:
#   lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
# Run single rep:
# lapply(1, random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)

if (number_parallel_threads == 1){
  exp1_op_list <- lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
} else {
  exp1_op_list <- mclapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params,
           mc.cores = number_parallel_threads)
}

# Change output file names from list to dataframe
exp1_op_df <- as.data.frame(do.call(rbind, exp1_op_list))
exp1_op_df_path <- paste0(local_directory, "exp1_file_output_paths.csv")
write.csv(exp1_op_df, file = exp1_op_df_path, row.names = TRUE)



## Experiment 2: Mimicking ILS and introgression ##
# Generate gene trees using ms to simulate ILS with a single introgression event
#   - Fix alignment at total_alignment_length
#   - Simulate gene trees for tree containing a single introgression event in ms. Vary proportion of introgressed DNA from 0 to 1
#   - Simulate DNA along each tree with Alisim
#   - Concatenate alignments

# Create folder to store results of this experiment, if it doesn't already exist
exp2_dir <- paste0(local_directory, "exp_2/")
if (dir.exists(exp2_dir) == FALSE){dir.create(exp2_dir)}

exp2_params <- expand.grid("num_reps" = number_of_replicates, "num_taxa" = number_of_taxa, "num_trees" = number_gene_trees, 
                           "tree_depth" = tree_depth_coalescent, "recombination_value" = r_vec, "recombination_type" = c("Ancient","Recent"))
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`_`tree_depth`_`recombination proportion`_`introgression event type`
exp2_params$uid <- paste0("exp2_",sprintf("%05d", exp2_params$num_trees), "_", sprintf("%04d", exp2_params$num_taxa), "_",
                          sprintf("%03d", exp2_params$num_reps), "_", exp2_params$tree_depth, "_", exp2_params$recombination_value,
                          "_", exp2_params$recombination_type)
# Add parameters for Alisim
exp2_params$alisim_gene_models <- alisim_gene_models
exp2_params$alisim_gene_tree_length <- alisim_gene_tree_length
# Add other parameters
exp2_params$total_alignment_length <- number_gene_trees * gene_length
exp2_params$sequence_type <- sequence_type
# Add name for the partition file and output alignment file for each simulated alignment
exp2_params$partition_file <- paste0(exp2_params$uid, "_partitions.nex")
exp2_params$output_alignment_file <- paste0(exp2_params$uid, "_output_alignment")

# Remove any rows that have 5 taxa and an ancient introgression event
# Due to the way ancient introgression events are structured (take place at time where four taxa exist, between two non-sister taxa), a recent and an ancient
#     introgression event will be identical for a tree with 5 taxa
remove_rows <- which(exp2_params$num_taxa == 5 & exp2_params$recombination_type == "Ancient")
keep_rows <- setdiff(1:nrow(exp2_params), remove_rows)
exp2_params <- exp2_params[keep_rows, ]


# Write exp2_params dataframe to file as a csv
exp2_df_path <- paste0(local_directory, "exp2_parameters.csv")
write.csv(exp2_params, file = exp2_df_path, row.names = TRUE)

# Iterate through each row in the parameters dataframe and generate an alignment for each set of parameters
# Run all reps: 
#   lapply(1:nrow(exp2_params), ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params)
# Run single rep:
#   lapply(1, ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params, select.sister = FALSE)
if (number_parallel_threads == 1){
  exp2_op_list<- lapply(1:nrow(exp2_params), ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params)
} else {
  exp2_op_list <- mclapply(1:nrow(exp2_params), ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params,
           mc.cores = number_parallel_threads)
}

# Change output file names from list to dataframe
exp2_op_df <- as.data.frame(do.call(rbind, exp2_op_list))
exp2_op_df_path <- paste0(local_directory, "exp2_file_output_paths.csv")
write.csv(exp2_op_df, file = exp2_op_df_path, row.names = TRUE)


## Ideas for extension:
#   -  Repeat above experiments but add random noise
#   -  Repeat above experiments but add alignment error

