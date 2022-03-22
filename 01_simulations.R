# /caitlinch/treelikeness_metrics/01_simulations.R
# Caitlin Cherryh 2022

# This program will simulate alignments with varying levels of treelikeness
# This program requires IQ-Tree2 (2.2-beta or above) and ms.



#### 1. Open packages ####
library(ape)
library(phytools)



#### 2. Set parameters ####
# local_directory                 <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory                  <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# ms_path                         <- Path to ms executable 
# iqtree2_path                    <- Path to IQ-Tree2 executable (version 2.2-beta or later to ensure Alisim is included). 
# total_alignment_length          <- Total length of concatenated alignments (we chose 10000).
# sequence_type                   <- Sequence type for simulation (we chose "DNA").
# taxa_vec                        <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                        <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# tree_depth_random               <- One or more values for the tree length of randomly generated trees.
# tree_depth_coalescent_bounds    <- Two values for minimum and maximum tree length of coalescent trees.
# r_vec                           <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models              <- model of sequence evolution for Alisim 
# alisim_gene_tree_length         <- gene-specific tree length for Alisim

run_location = "soma"
if (run_location == "local"){
  local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
  repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
  ms_path <- "ms"
  iqtree2_path <- "iqtree2.2-beta"
} else if (run_location == "soma"){
  local_directory <- "/data/caitlin/treelikeness_metrics/"
  repo_directory <- "/data/caitlin/treelikeness_metrics/code/"
  ms_path <- "/data/caitlin/executables/msdir/ms"
  iqtree2_path <- "/data/caitlin/linux_executables/iqtree-2.1.2-Linux/bin/iqtree2"
}

total_alignment_length <- 10000
sequence_type <- "DNA"
taxa_vec <- c(10,20,50,100,200,500,1000)
num_reps <- 10
tree_depth_random <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5)
tree_depth_coalescent_bounds <- c(0.25, 70) # where bounds for coalescent tree depth are 0.25 (minimum) and 70 (maximum) in coalescent units
r_vec <- seq(0, 0.5, 0.05)
alisim_gene_models <- "JC"
alisim_gene_tree_length <- NA

## Prepare variables that remain stable for all experiments using the input parameters
# Determine number of trees - all whole numbers that are a divisor of the total alignment length
number_of_trees <- divisors(total_alignment_length)
# Set number of taxa equal to taxa_vec
number_of_taxa <- taxa_vec
# Create a list of all replicate numbers using the num_reps value
number_of_replicates <- 1:num_reps
# Prepare vector for depth of coalescent trees
tree_depth_coalescent <- seq(floor(tree_depth_coalescent_bounds[1]),ceiling(tree_depth_coalescent_bounds[2]), (ceiling(tree_depth_coalescent_bounds[2])/10))
tree_depth_coalescent[1] <- tree_depth_coalescent_bounds[1] # Replace first value with actual first value


#### 3. Source functions from caitlinch/treelikeness_metrics ####
source(paste0(repo_directory, "func_simulating_alignments.R"))



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
exp1_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees, tree_depth_random)
names(exp1_params) <- c("num_reps", "num_taxa", "num_trees", "tree_depth")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
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

# Iterate through each row in the parameters dataframe
# lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)
lapply(1, random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path, experiment_params = exp1_params)


## Experiment 2: ILS ##
# Generate gene trees using ms to simulate ILS
#     Total alignment length = 10,000 bp
#     Number of trees ranges from 0 to 10,000  (whole number divisors of 10000)
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 1000
# Simulate DNA along each tree with Alisim, using the topology-unlinked partition model

# Create folder to store results of this experiment, if it doesn't already exist
exp2_dir <- paste0(local_directory, "exp_2/")
if(!file.exists(exp2_dir)){dir.create(exp2_dir)}

# Prepare parameters for experiments
# This experiment differs because the number of NNI moves is a variable included in the expand.grid function
# Create matrix with parameters for generating each simulated alignment
exp2_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees, tree_depth_coalescent)
names(exp2_params) <- c("num_reps", "num_taxa", "num_trees", "tree_depth")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
exp2_params$uid <- paste0("exp2_",sprintf("%05d", exp2_params$num_trees), "_", sprintf("%04d", exp2_params$num_taxa), "_",
                          sprintf("%03d", exp2_params$num_reps), "_", exp2_params$tree_depth)
# Add parameters for Alisim
exp2_params$alisim_gene_models <- alisim_gene_models
exp2_params$alisim_gene_tree_length <- alisim_gene_tree_length
exp2_params$recombination_value <- 0
exp2_params$recombination_type <- NA
# Add other parameters
exp2_params$total_alignment_length <- total_alignment_length
exp2_params$sequence_type <- sequence_type
# Add name for the partition file and output alignment file for each simulated alignment
exp2_params$partition_file <- paste0(exp2_params$uid, "_partitions.nex")
exp2_params$output_alignment_file <- paste0(exp2_params$uid, "_output_alignment")

# Iterate through each row in the parameters dataframe
# lapply(1:nrow(exp2_params), ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params)
lapply(1, ms.generate.alignment, output_directory = exp2_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp2_params)


## Experiment 3: Mimicking introgression ##
# fix alignment at total_alignment_length
# simulate tree containing a single introgression event in ms and vary proportion of introgressed DNA from 0 to 1
# simulate DNA along each tree with Alisim
# Concatenate alignments

# Create folder to store results of this experiment, if it doesn't already exist
exp3_dir <- paste0(local_directory, "exp_3/")
if (dir.exists(exp3_dir) == FALSE){dir.create(exp3_dir)}

exp3_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees, tree_depth_coalescent, r_vec, c("Ancient","Recent"))
names(exp3_params) <- c("num_reps", "num_taxa", "num_trees", "tree_depth", "recombination_value", "recombination_type")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
exp3_params$uid <- paste0("exp3_",sprintf("%05d", exp3_params$num_trees), "_", sprintf("%04d", exp3_params$num_taxa), "_",
                          sprintf("%03d", exp3_params$num_reps), "_", exp3_params$tree_depth, "_", exp3_params$recombination_value,
                          "_", exp3_params$recombination_type)
# Add parameters for Alisim
exp3_params$alisim_gene_models <- alisim_gene_models
exp3_params$alisim_gene_tree_length <- alisim_gene_tree_length
# Add other parameters
exp3_params$total_alignment_length <- total_alignment_length
exp3_params$sequence_type <- sequence_type
# Add name for the partition file and output alignment file for each simulated alignment
exp3_params$partition_file <- paste0(exp3_params$uid, "_partitions.nex")
exp3_params$output_alignment_file <- paste0(exp3_params$uid, "_output_alignment")

# Iterate through each row in the parameters dataframe and generate an alignment for each set of parameters
# lapply(1:nrow(exp3_params), ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp3_params)
lapply(1, ms.generate.alignment, output_directory = exp3_dir, ms_path = ms_path, iqtree2_path = iqtree2_path, experiment_params_df = exp3_params)


## Experiment 4: Repeat above experiments but adding random noise ##



## Experiment 5: Repeat above experiments adding alignment error ##



