# Open packages
library(ape)
library(phytools)

# Parameters for all simulations
# local_directory         <- Directory where alignments will be saved/treelikeness metrics will be run.
# repo_directory          <- Location of caitlinch/treelikeness_metrics github repository (for access to functions).
# iqtree2_path            <- Path to IQ-Tree2.2-beta executable (this is the IQ-Tree2 release containing Alisim). 
# total_alignment_length  <- Total length of concatenated alignments (we chose 10000).
# sequence_type           <- Sequence type for simulation (we chose "DNA").
# taxa_vec                <- Number of taxa to simulate (we chose 10,20,50,100,200,500, and 1000).
# num_reps                <- Number of replicates to run for each set of simulation conditions (we chose 10). Must be >= 1.
# r_vec                   <- Values of introgression (we chose from 0 to 1 in intervals of 0.05)
# alisim_gene_models      <- model of sequence evolution for Alisim 
# alisim_gene_tree_length <- gene-specific tree length for Alisim

local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
ms_path <- "ms"
iqtree2_path <- "iqtree2.2-beta"
total_alignment_length <- 10000
sequence_type <- "DNA"
taxa_vec <- c(10,20,50,100,200,500,1000)
num_reps <- 10
r_vec <- seq(0, 1, 0.05)
alisim_gene_models <- NA
alisim_gene_tree_length <- NA

# test params
output_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/testing_metrics/testing_ms/"
ntaxa = 10
ntrees = 20
replicate_number = 1

## Source functions from caitlinch/treelikeness_metrics
source(paste0(repo_directory, "func_simulating_alignments.R"))



## Prepare parameters that remain stable for all experiments
# Prepare parameters for experiments
# The number of trees should be the divisors for total_alignment_length (as only whole numbers of trees are possible)
number_of_trees <- divisors(total_alignment_length)
number_of_taxa <- taxa_vec
number_of_replicates <- 1:num_reps



## Experiment 1: Random trees
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
exp1_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees)
names(exp1_params) <- c("num_reps", "num_taxa", "num_trees")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
exp1_params$uid <- paste0("exp1_",sprintf("%05d", exp1_params$num_trees), "_", sprintf("%04d", exp1_params$num_taxa), "_",
                          sprintf("%03d", exp1_params$num_reps))
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
lapply(1:nrow(exp1_params), random.trees.generate.alignment, output_directory = exp1_dir, iqtree2_path = iqtree2_path,
       experiment_params = exp1_params)



## Experiment 2: ILS
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
exp2_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees, r_vec)
names(exp2_params) <- c("num_reps", "num_taxa", "num_trees", "recombination_value")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
exp2_params$uid <- paste0("exp2_",sprintf("%05d", exp2_params$num_trees), "_", sprintf("%04d", exp2_params$num_taxa), "_",
                          sprintf("%03d", exp2_params$num_reps))
# Add parameters for Alisim
exp2_params$alisim_gene_models <- alisim_gene_models
exp2_params$alisim_gene_tree_length <- alisim_gene_tree_length
# Add other parameters
exp2_params$total_alignment_length <- total_alignment_length
exp2_params$sequence_type <- sequence_type
# Add name for the partition file and output alignment file for each simulated alignment
exp2_params$partition_file <- paste0(exp2_params$uid, "_partitions.nex")
exp2_params$output_alignment_file <- paste0(exp2_params$uid, "_output_alignment")

# Iterate through each row in the parameters dataframe
lapply(1:nrow(exp2_params), ILS.generate.alignment, output_directory = exp2_dir, iqtree2_path = iqtree2_path,
       experiment_params = exp2_params)


#### Run one replicate of experiment 2: ILS ####
# Generate gene trees in ms from random starting tree
ms_output_files <- ms.generate.trees(ntaxa, ntrees, output_directory, ms_path, replicate_number)
gene_trees_file <- ms_output_files[3]
# Generate the partition file
gene_partition_file <- gsub("starting_tree", "partitions", ms_output_files[[1]])
partition.random.trees(ntrees, total_alignment_length, sequence_type, models = NA, rescaled_tree_lengths = NA, output_filepath = gene_partition_file)
# Generate DNA along gene trees
output_alignment_file <- gsub("starting_tree\\.txt", "output_alignment.fasta", ms_output_files[[1]])
alisim.topology.unlinked.partition.model(iqtree_path = iqtree2_path, output_alignment_path = output_alignment_file, partition_file_path = gene_partition_file, 
                                         trees_path = gene_trees_file, output_format = "fasta", sequence_type = "DNA")

## Experiment 3: Mimicking introgression

# fix alignment at total_alignment_length

# simulate tree containing a single introgression event in ms and vary proportion of introgressed DNA from 0 to 1

# simulate DNA along each tree with Alisim

# Concatenate alignments


## Experiment 4: Repeat above experiments but adding random noise

## Experiment 5: Repeat above experiments adding alignment error