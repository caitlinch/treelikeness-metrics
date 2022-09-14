# /caitlinch/treelikeness_metrics/code/func_empirical_datasets.R
# Caitlin Cherryh 2022

# This file contains functions for manipulating empirical alignments or working with empirical datasets

library(ape)
library(phylotools)



remove.unnecessary.taxa <- function(gene_path, taxa_to_remove, id){
  # Small function to take a single gene, remove only the relevant taxa, and save the reduced matrix
  
  # Make new file name for gene
  gene_new_id <- gsub(".fa", "_", gene_path)
  gene_new_path <- paste0(gene_new_id, id, ".fa")
  
  # Make a new alignment with only the desired taxa included
  rm.sequence.fasta(infile = gene_path, outfile = gene_new_path, to.rm = taxa_to_remove)
}


alignment.dimensions.nex <- function(alignment_path){
  # Small function to take a nexus alignment path and return the dimensions of that alignment
  
  # Open nexus file
  n <- read.nexus.data(alignment_path)
  # Get the dimensions of the matrix
  n_taxa <- length(n)
  n_chars <- unique(lapply(n, length))[[1]]
  # Create output vector
  op_vector <- c(basename(alignment_path), n_taxa, n_chars)
  # Return output
  return(op_vector)
}