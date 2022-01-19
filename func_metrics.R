#### caitlinch/treelikeness_metrics/func_metrics.R
## This file contains functions to apply treelikeness metrics to a single alignment

## Load required packages
library(ape)

# here's paths for different programs needed for test statistics:
iqtree2_path <- "iqtree2.2-beta"

# here's a file path to a test alignment (one tree, 10000bp, 20 taxa - should be treelike):
al_tl_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00001_0020_001_output_alignment.fa"
alignment_path <- al_tl_path
sequence_format = "DNA"
substitution_model = "raw"

# here's a few test alignments with 20 taxa each, with either 1, 10, 100, 1000, or 10000 trees:
test_paths <- paste0("/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/", 
                     c("exp1_00001_0020_001_output_alignment.fa", "exp1_00010_0020_001_output_alignment.fa",
                       "exp1_00100_0020_001_output_alignment.fa", "exp1_01000_0020_001_output_alignment.fa"))



## Site concordance factors (Minh et. al. 2020)
scf <- function(alignment_path, iqtree2_path, iqtree2_number_threads, number_scf_quartets = 100){
  # Function to calculate the site concordance factors for an alignment, given a maximum likelihood tree
  
  ## Check that the treefile already exists: if it doesn't, run IQ-Tree and create it
  if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
    # Given an alignment, estimate the maximum likelihood tree with IQ-Tree2
    call <- paste0(iqtree2_path," -s ",alignment_path," -nt ", iqtree2_number_threads, " -redo -safe")
    system(call)
  }
  ## Check if the site concordance factors have already been calculated: if they have not, calculate them
  if (file.exists(paste0(alignment_path,".treefile.cf.stat")) == FALSE){
    # Create the command and call it in the system
    # for sCF: iqtree -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -nt 10
    treefile <- paste0(alignment_path,".treefile")
    call <- paste0(iqtree_path," -t ",treefile," -s ",alignment_path," --scf ",num_scf_quartets," -nt ","1"," -redo -safe")
    system(call) # call IQ-tree!
  }
  ## Retrieve the site concordance factors from the output table
  scf_table <- read.table(paste0(alignment_path,".treefile.cf.stat"), header = TRUE, sep = "\t")
  scf_extracts <- list(mean_scf = round(mean(scf_table$sCF), digits = 2), 
                       median_scf = round(median(scf_table$sCF), digits = 2), 
                       all_scfs = scf_table$sCF, 
                       branch_ids = scf_table$ID)
  ## Return the site concordance factor results
  return(scf_results)
}



## Delta plots (Holland et. al. 2002)
mean.delta.plot.value <- function(alignment_path, sequence_format = "DNA", substitution_model = "raw"){
  # This function takes an alignment, calculates a distance matrix for the alignment, and the applies the
  # `ape` function `delta.plot`. We take the mean delta plot value as the test statistic. 
  
  ## Open the alignment as a DNAbin object
  alignment <- read.FASTA(alignment_path, type = sequence_format)
  ## Calculate a distance matrix of pairwise distances from DNA sequences using a model of DNA substitution
  # Default model of DNA substitution is the default model for the `ape` function `dist.dna` ("K80")
  pdm <- dist.dna(alignment, model = substitution_model)
  ## Call ape::delta.plot function
  # Set the number of intervals for the delta plot
  dp_intervals = 100
  # Make a delta.plot based on the pairwise distance matrix
  dp <- delta.plot(pdm, k = dp_intervals, plot = FALSE)
  ## To calculate the mean delta q from ALL quartets:
  # Create two vectors, one containing the counts and one containing the midpoint of each interval
  # To determine the midpoint of each interval, first find the intervals (e.g. for k = 2, there will be 2 intervals: 0-0.5 and 0.5-1),
  #     and the midpoint of each interval will be 0.25 and 0.75 (the mean of the start and endpoint of each interval)
  interval_midpoint = (seq(0,0.999,1/(dp_intervals)) + (0.5 * (0 + seq(0,1,1/(dp_intervals))[2])))
  interval_count = dp$counts
  # To calculate the mean delta_q, calculate a weighted mean from the 
  mean_dq <- weighted.mean(interval_midpoint, interval_count)
  ## To calculate the mean delta bar:
  # Calculate the mean delta bar (delta bar = the mean delta value for each observation/taxa)
  mean_db <- mean(dp$delta.bar)
  ## Return values to outside function
  # Return the mean delta bar (the mean delta q value across all taxa)
  return(mean_db)
}

