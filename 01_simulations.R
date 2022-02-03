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
# alisim_gene_models      <- model of sequence evolution for Alisim 
# alisim_gene_tree_length <- gene-specific tree length for Alisim

local_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness_metrics/"
iqtree2_path <- "iqtree2.2-beta"
total_alignment_length <- 10000
sequence_type <- "DNA"
taxa_vec <- c(10,20,50,100,200,500,1000)
num_reps <- 10
alisim_gene_models <- NA
alisim_gene_tree_length <- NA



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
#     Number of trees ranges from 0 to 10,000 in intervals of 100
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



## Experiment 2: Related trees
# Generate x related trees with y taxa 
#     Start by generating 1 random tree with y taxa, then perform a single NNI move on the starting tree until x related trees have been generated
#     Total alignment length = 10,000 bp
#     Number of trees ranges from 0 to 10,000 in intervals of 100
#     Length of alignment for each tree is total alignment length divided by the number of trees
#     Number of taxa varies from 10 to 1000

# Create folder to store results of this experiment, if it doesn't already exist
exp2_dir <- paste0(local_directory, "exp_2/")
if(!file.exists(exp2_dir)){dir.create(exp2_dir)}

# Prepare parameters for experiments
# This experiment differs because the number of NNI moves is a variable included in the expand.grid function
number_of_NNI_moves = 1
# Create matrix with parameters for generating each simulated alignment
exp2_params <- expand.grid(number_of_replicates, number_of_taxa, number_of_trees, number_of_NNI_moves)
names(exp2_params) <- c("num_reps", "num_taxa", "num_trees", "num_NNI_moves")
# Add a unique identifier (uid) of the form: experiment_`number of trees`_`number of taxa`_`replicate number`
exp2_params$uid <- paste0("exp2_",sprintf("%05d", exp2_params$num_trees), "_", sprintf("%04d", exp2_params$num_taxa), "_",
                          sprintf("%03d", exp2_params$num_reps))
# Add parameters for Alisim
exp2_params$alisim_gene_models <- alisim_gene_models
exp2_params$alisim_gene_tree_length <- alisim_gene_tree_length
# Add other parameters
exp2_params$total_alignment_length <- total_alignment_length
exp2_params$sequence_type <- sequence_type
# Add names for the tree file, partition file and output alignment file for each simulated alignment
exp2_params$tree_file <- paste0(exp2_params$uid, "_random_trees.phy")
exp2_params$partition_file <- paste0(exp2_params$uid, "_partitions.nex")
exp2_params$output_alignment_file <- paste0(exp2_params$uid, "_output_alignment")

# Iterate through each row in the parameters dataframe
lapply(1:nrow(exp2_params), ILS.generate.alignment, output_directory = exp2_dir, iqtree2_path = iqtree2_path,
       experiment_params = exp2_params)


#### Generate one ms command ####
# generate random tree
ntaxa = 5
ntrees = 5
# Generate random coalescent tree 
t <- rcoal(ntaxa)
ms_coal_ints <- calculate.ms.coalescent.times(t$Nnode, coalescent.intervals(t))
# Determine order of coalescence by getting which taxa are in a clade derived from each node
nodes <- (ntaxa+1):(ntaxa+t$Nnode)
# Extract information about all clades from tree
node_df <- do.call(rbind.data.frame, lapply(nodes, extract.clade.from.node, tree = t, coalescent_times = ms_coal_ints))
names(node_df) <- c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "coalescence_time", "removed_taxa", "ms_input")

# Determine ms command line using coalescent intervals (earlier coalescent events to the left, tips towards the right)
coalescence_command <- paste0("./ms ", ntaxa, " ", ntrees, " -T -I ", ntaxa," ", paste(rep(1, ntaxa), collapse = " "), " ",
                              "-ej 2.326 3 1 -ej 0.294 2 1 -ej 0.166 5 1 -ej 0.004 4 2")

coalescence_command <- paste0("./ms ", ntaxa, " ", ntrees, " -T -I ", ntaxa," ", paste(rep(1, ntaxa), collapse = " "), " ",
                              "-ej ", ms_coal_ints[4], " 2 1 ", "-ej ", ms_coal_ints[3], " 3 2 ", "-ej ", ms_coal_ints[2], " 1 4 ","-ej ", ms_coal_ints[1], " 5 2 ")




## Experiment 3: Mimicking introgression

# fix alignment at total_alignment_length

# simulate tree containing a single introgression event in ms and vary proportion of introgressed DNA from 0 to 1

# simulate DNA along each tree with Alisim

# Concatenate alignments


## Experiment 4: Repeat above experiments but adding random noise

## Experiment 5: Repeat above experiments adding alignment error