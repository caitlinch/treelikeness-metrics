# caitlinch/treelikeness_metrics/code/func_parametric_bootstraps.R
# Caitlin Cherryh 2023

# This file contains functions to generate parametric bootstraps from an empirical alignment
# Some functions require IQ-Tree2 (2.2-beta or above)

#### Parametric bootstrap wrapper ####
complete.one.bootstrap.replicate <- function(row_id, bootstrap_df, iqtree2_path){
  # Perform one full bootstrap replicate 
  
  ## Prepare for analysis
  # Extract row
  bs_row <- bootstrap_df[row_id, ]
  # Make replicate folder
  dir.create(bs_row$replicate_folder)
  
  ## Create alignment
  
  
  ## Run treelikeness tests
  
  ## Assemble output
  
  
}



#### Generate alignment ####
generate.mimic.alignment <- function(iqtree_path, output_alignment_path, trees_path, 
                                     output_format = "fasta", sequence_type = "AA"){
  ## This function uses the topology-unlinked partition model in Alisim to generate a sequence alignment
  #     containing multiple concatenated genes, each with its own tree topology and branch lengths
  
  # Assemble function call 
  function_call <- paste0(iqtree_path, " --alisim ", output_alignment_path, " -t ", trees_path, 
                          " --seqtype ", sequence_type, " -af ", output_format)
  # Invoke the OS command and call IQ-Tree
  system(function_call)
  
  # Print completion statement
  print("Alisim (IQ-Tree2) run complete")
}


