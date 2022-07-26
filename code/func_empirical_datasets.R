# /caitlinch/treelikeness_metrics/code/func_empirical_datasets.R
# Caitlin Cherryh 2022

# This file contains functions for manipulating empirical alignments or working with empirical datasets



library(phylotools)



remove.unnecessary.taxa <- function(gene_path, taxa_to_remove, id){
  # Small function to take a single gene, remove only the relevant taxa, and save the reduced matrix
  
  # Make new file name for gene
  gene_new_id <- gsub(".fa", "_", gene_path)
  gene_new_path <- paste0(gene_new_id, id, ".fa")
  
  # Make a new alignment with only the desired taxa included
  rm.sequence.fasta(infile = gene_path, outfile = gene_new_path, to.rm = taxa_to_remove)
}