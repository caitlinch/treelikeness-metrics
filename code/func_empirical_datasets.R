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
  op_vector <- c(alignment_path, n_taxa, n_chars)
  # Return output
  return(op_vector)
}



process.alignment.dimensions.row <- function(row){
  # Small function to take a row from the alignment_dimensions csv, shuffle the taxa labels, update the row and return the updated row
  
  # Copy the information into a new row
  new_row <- row
  # Change the uid from "`clade`_`gene`_copy" to "`clade`_`gene`_shuffled" (where clade and gene are properties of the alignment)
  new_row$uid <- gsub("copy", "shuffled_all", new_row$uid)
  # Change the shuffled_taxa column in the new row to TRUE
  new_row$shuffled_taxa <- "TRUE"
  
  # Copy the alignment and shuffle the taxa of the new alignment randomly
  new_alignment_path <- shuffle.taxa(new_row$alignment_path)
  
  # Return output
  return(new_row)
}



shuffle.taxa <- function(alignment_path){
  # Small function to take a file path to a nexus alignment, copy the file, and shuffle the taxa randomly
  
  # Split the alignment path at the file extension
  alignment_path_chunks <- strsplit(basename(alignment_path), "\\.")[[1]]
  # Create a file path for the new shuffled alignment
  shuffled_alignment <- paste0(dirname(alignment_path), "/", alignment_path_chunks[1], "_shuffled_all", ".", alignment_path_chunks[2])
  # Open the copied alignment
  n <- read.nexus.data(alignment_path)
  # Extract the taxa names
  taxa_names <- names(n)
  # Shuffle the taxa names 
  shuffled_taxa <- sample(taxa_names)
  # Replace the taxa names with the shuffled taxa
  shuffled_n <- n
  names(shuffled_n) <- shuffled_taxa
  # Save the nexus file to the shuffled_alignment path
  write.nexus.data(shuffled_n, file = shuffled_alignment, format = "dna", interleaved = FALSE)
  
  # Identify the number of taxa that remain in the same position
  kept_taxa <- length(which(names(n) == names(shuffled_n)))
  percent_kept_taxa <- kept_taxa/length(n)*100
  
  # Assemble the output vector
  op_vector <- c(shuffled_alignment, kept_taxa, percent_kept_taxa)
  names(op_vector) = c("shuffled_alignment_path", "number_unchanged_taxa", "percent_unchanged_taxa")
  
  # Return the file path of the shuffled alignment 
  return(op_vector)
}




