#### caitlinch/treelikeness_metrics/func_metrics.R
## This file contains functions to apply treelikeness metrics to a single alignment

## Load required packages
library(ape) # for general tree/alignment wrangling, and the delta.plots function
library(ips) # to determine the indices of the parsimony informative sites
library(phangorn) # for splits and networks

# here's paths for different programs needed for test statistics:
iqtree2_path <- "iqtree2"
fast_TIGER_path <- "/Users/caitlincherryh/Documents/Executables/fast_TIGER-0.0.2/DAAD_project/fast_TIGER"
phylogemetric_path <- "/Users/caitlincherryh/Documents/Executables/phylogemetric/phylogemetric_executable"
splitstree_path <- "/Users/caitlincherryh/Documents/Executables/SplitsTree/SplitsTree.app/Contents/MacOS/JavaApplicationStub"
reticulation_index_path <- "/Users/caitlincherryh/Documents/Executables/Coalescent_simulation_and_gene_flow_detection-master/"

# here's a file path to a test alignment (one tree, 10000bp, 20 taxa - should be treelike):
al_tl_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/exp1_00001_0020_001_output_alignment.fa"

# here's a few test alignments with 20 taxa each, with either 1, 10, 100, 1000, or 10000 trees:
test_paths <- paste0("/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/exp_1/", 
                     c("exp1_00001_0020_001_output_alignment.fa", "exp1_00010_0020_001_output_alignment.fa",
                       "exp1_00100_0020_001_output_alignment.fa", "exp1_01000_0020_001_output_alignment.fa"))

# here's paths for variables needed to test treelikeness metric functions
alignment_path <- al_tl_path
alignment_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/testing_metrics/testing_splitstree4/test.phy"
alignment_path <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/testing_metrics/testing_reticulation_index/exp1_00100_0020_001_output_alignment.fa"
sequence_format = "DNA"
substitution_model = "raw"
iqtree2_number_threads = "AUTO"
phylogemetric_number_of_threads = NA
number_scf_quartets = 100
number_of_taxa = 20

#### Treelikeness metric functionc ####

## Reticulation Index (Cai et. al. 2021)
reticulation_index <- function(alignment_path, reticulation_index_path, iqtree2_path, sequence_format = "DNA"){
  ## Reticulation index quantifies introgression at each node
  # Requires three inputs: 
  # 1. Rooted best-estimated species tree with branch lengths estimated in coalescent units, newick format
  #       - Obtained from running MPEST/ASTRAL on ML gene trees
  # 2. Rooted bootstrap species trees with branch lengths measured in coalescent units, newick format
  #       - Obtained from running MPEST/ASTRAL on bootstrap gene trees
  # 3. Rooted gene trees, newick format
  #       - Obtained from running IQ-Tree on each partition from the alignment file
  
  ## Determine file names
  partition_path <- gsub("output_alignment.fa", "partitions.nex", alignment_path)
  
  ## Estimating gene trees
  # Create a new folder to store all information about this alignment
  gene_folder <- paste0(dirname(alignment_path), "/gene_folder/")
  if (dir.exists(gene_folder) == FALSE){dir.create(gene_folder)}
  # Run function to separate alignment into genes using partition file
  gene_info <- genes.from.alignment(alignment_path, partition_path, gene_folder, sequence_format)
  # Estimate a maximum likelihood gene tree for each gene in IQ-Tree2
  estimate.iqtree2.gene.trees(gene_folder, iqtree2_path)
}



## Network tree-likeness test (Huson and Bryant 2006)
network.treelikeness.test <- function(alignment_path, splitstree_path, sequence_format = "DNA"){
  ## Uses Splitstree4.17.1 software to implement the Network Treelikeness Test described in Huson and Bryant (2006)
  # Software available from:
  # https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/
  
  ## Test steps:
  #   1. Infer a split network N from some data
  #         - Can be done in Splitstree using Networks -> NeighborNet
  #   2. Construct a confidence network for N with level 1 - alpha
  #         - Construct the set of splits with CIs excluding 0
  #         - Can be done in Splitstree4 using `Analysis -> Bootstrap...` then `Analysis -> Show Confidence Network`
  #   3. If the confidence network does not contain a tree, reject the null hypothesis that the data originated in a tree
  #         - Reject the null hypothesis if and only if this set is incompatible
  
  ## Construct and bootstrap a NeighborNet network
  # Convert fasta to nexus
  nexus_alignment_path <- convert.to.nexus(alignment_path, sequence_format = "DNA")
  # Name output path
  output_path <- paste0(alignment_path, "_Splitstree_output.nex")
  # Assemble the SplitsTree 4 command
  splitstree_command <- paste0(splitstree_path, " -g -x 'OPEN FILE=", nexus_alignment_path,"; ASSUME chartransform=Uncorrected_P HandleAmbiguousStates=Ignore Normalize=true; ASSUME disttransform=NeighbourNet; BOOTSTRAP RUNS=100; EXPORT FILE=", output_path," REPLACE=YES; QUIT'")
  # Call SplitsTree 4
  system(splitstree_command)
  
  ## Construct a confidence network using the bootstrap splits
  # Read in the nexus splits from the Splitstree output file
  splits <- read.nexus.splits(output_path)
  # Extract bootstrap confidence values from Splitstree4 output
  splitstree_text <- readLines(output_path)
  st_bootstraps_start <- grep("BEGIN st_bootstrap", splitstree_text)
  st_bootstraps_end <- grep("END; \\[st_Bootstrap\\]", splitstree_text)
  st_bootstraps_all <- splitstree_text[st_bootstraps_start:st_bootstraps_end]
  # Trim bootstrap splits to only splits present in network
  matrix_ind <- grep("MATRIX", st_bootstraps_all)
  line_ind <- grep("\\[--------------------------------------------------\\]", st_bootstraps_all)
  st_bootstraps <- st_bootstraps_all[(matrix_ind+1):(line_ind-1)]
  # Extract bootstrap confidence values from Splitstree4 output
  split_stbs <- unlist(strsplit(st_bootstraps, "\t"))
  CN_splits_df <- split_stbs[c(FALSE, TRUE, FALSE)]
  
}



## Q-residuals (Gray et. al. 2010)
q_residuals <- function(alignment_path, phylogemetric_path, sequence_format = "DNA", phylogemetric_number_of_threads = NA){
  ## Uses software`phylogemetric` (available at https://github.com/SimonGreenhill/phylogemetric) to calculate
  #     the Q-residuals of Gray et. al. (2010)
  # Outputs a value between 0 and 1, where 0 indicates the quartets are treelike
  
  ## Change fasta file (from alisim) into nexus file
  # Nexus format required for `phylogemetric` software
  f <- read.FASTA(file = alignment_path, type = sequence_format)
  nexus_path <- paste0(alignment_path, ".nex")
  if (sequence_format == "DNA"){
    nexus_format = "dna"
  } else if ((sequence_format == "AA") | (sequence_format == "Protein")){
    nexus_format = "protein"
  }
  write.nexus.data(f, nexus_path, format = nexus_format, interleaved = FALSE)
  # Remove the "INTERLEAVE = FALSE" section (so other programs aren't confused)
  txt <- readLines(nexus_path)
  ind <- grep("FORMAT", txt)
  txt[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=-;"
  write(txt, file = nexus_path, append = FALSE)
  
  ## Apply `phylogemetric` and save output as a vector
  if (is.na(phylogemetric_number_of_threads)){
    # Run basic command
    call <- paste0(phylogemetric_path, " qresidual ", nexus_path)
  } else {
    # Use multiple cores with the -w/--workers argument
    call <- paste0(phylogemetric_path, " -w ", phylogemetric_number_of_threads, " qresidual ", nexus_path)
  }
  program_output <- system(call, intern = TRUE)
  # Format the program output into a nice vector containing just the mean Q-residual value per taxa
  q_residuals <- as.numeric(unlist(strsplit(program_output, "\t"))[c(FALSE, TRUE)])
  
  ## Return mean Q-residual value
  mean_q_residual <- mean(q_residuals)
  return(mean_q_residual)
}


## TIGER (Cummins and McInerney 2011)
TIGER <- function(alignment_path, fast_TIGER_path, sequence_format = "DNA"){
  # Tree Independent Genertion of Evolutionary Rates
  # Function to calculate TIGER values from a multiple sequence alignment using the 
  #   pbfrandsen/fast_TIGER software (available here: https://github.com/pbfrandsen/fast_TIGER)
  
  ## Remove phylogenetically uninformative sites from the alignment
  # Uninformative sites can bias the TIGER values - see List (2021)
  alignment <- read.FASTA(alignment_path, type = sequence_format)
  # Detect informative sites (the sites to keep)
  is_inds <- pis(as.matrix(alignment), what = "index")
  # Index the alignment, to reduce it to just the informative sites
  informative_alignment <- as.matrix.DNAbin(alignment)[,c(is_inds)]
  # Write the informative alignment to disk
  informative_alignment_path <- paste0(alignment_path, "_ParsimonyInformativeSites_only.phy")
  write.phy(informative_alignment, file = informative_alignment_path)
  
  ## Run fast_TIGER
  # Change sequence format to either "dna" or "aa" (only allowable data type names)
  if (sequence_format == "DNA"){
    data_type = "dna"
  } else if ((sequence_format == "AA") | (sequence_format == "Protein")){
    data_type = "protein"
  }
  # Create system command and call fast_TIGER
  call <- paste0(fast_TIGER_path, " ", data_type, " ", informative_alignment_path)
  system(call)
  
  ## Open the results file and output mean TIGER value
  # TIGER values range from 0 to 1, with values closer to 1 indicating more stable/consistent sites
  TIGER_file <- paste0(informative_alignment_path, "_r8s.txt")
  TIGER_values <- as.numeric(readLines(TIGER_file))
  # Calculate the mean TIGER value
  mean_TIGER_value <- mean(TIGER_values)
  # Return the mean value
  return(mean_TIGER_value)
}



## Likelihood mapping (Strimmer and von Haeseler 1997)
likelihood.mapping <- function(alignment_path, iqtree2_path, iqtree2_number_threads = 1, number_of_taxa){
  # Function to call IQ-Tree and create a likelihood map for the alignment
  
  ## Create the likelihood map
  # Check whether likelihood mapping or IQ-Tree have run before. 
  # If one or both haven't run IQ-Tree to create the likelihood map
  iq_file <- paste0(alignment_path, ".iqtree")
  map_file <- paste0(alignment_path, ".lmap.eps")
  if ((file.exists(iq_file) == FALSE) | (file.exists(map_file) == FALSE)){
    number_of_quartets <- 25 * as.numeric(number_of_taxa)
    call <- paste0(iqtree2_path," -s ",alignment_path," -nt ", iqtree2_number_threads, " -lmap ",number_of_quartets," -redo -safe")
    system(call)
  }
  
  ## Extract results from likelihood map
  iq_log <- readLines(iq_file)
  ind <- grep("Number of fully resolved  quartets",iq_log)
  resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  ind <- grep("Number of partly resolved quartets",iq_log)
  partly_resolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  ind <- grep("Number of unresolved",iq_log)
  unresolved_q <- as.numeric(strsplit(strsplit(iq_log[ind],":")[[1]][2],"\\(")[[1]][1])
  total_q <- (resolved_q+partly_resolved_q+unresolved_q)
  prop_resolved <- resolved_q/total_q
  
  ## Collate results into a vector
  lm_results <- c(resolved_q, partly_resolved_q, unresolved_q, total_q, prop_resolved)
  names(lm_results) <- c("num_resolved_quartets", "num_partly_resolved_quartets", "num_unresolved_quartets",
                         "total_num_quartets", "proportion_resolved_quartets")
  
  ## Return results
  return(lm_results)
}



## Site concordance factors (Minh et. al. 2020)
scf <- function(alignment_path, iqtree2_path, iqtree2_number_threads = "AUTO", number_scf_quartets = 100, add.likelihood.map = FALSE, number_of_taxa = NA){
  # Function to calculate the site concordance factors for an alignment, given a maximum likelihood tree estimated in IQ-Tree
  # Optional: can perform likelihood map (so tree doesn't have to be estimated multiple times to do both likelihood mapping and sCF)
  
  ## Check that the treefile already exists: if it doesn't, run IQ-Tree and create it
  if (add.likelihood.map == FALSE){
    if (file.exists(paste0(alignment_path,".treefile")) == FALSE){
      # Given an alignment, estimate the maximum likelihood tree with IQ-Tree2
      call <- paste0(iqtree2_path," -s ",alignment_path," -nt ", iqtree2_number_threads, " -redo -safe")
      system(call)
    }
  } else if (add.likelihood.map == TRUE){
    if (((file.exists(paste0(alignment_path, ".iqtree")) == FALSE) | (file.exists(paste0(alignment_path, ".lmap.eps")) == FALSE)) &
        (is.na(number_of_taxa) == FALSE)) {
      number_of_quartets <- 25 * as.numeric(number_of_taxa)
      call <- paste0(iqtree2_path," -s ",alignment_path," -nt ", iqtree2_number_threads, " -lmap ",number_of_quartets," -redo -safe")
      system(call)
    }
  }
  ## Check if the site concordance factors have already been calculated: if they have not, calculate them
  if (file.exists(paste0(alignment_path,".treefile.cf.stat")) == FALSE){
    # Create the command and call it in the system
    # for sCF: iqtree -t concat.treefile -s ALN_FILE --scf 100 --prefix concord -nt 10
    treefile <- paste0(alignment_path,".treefile")
    call <- paste0(iqtree2_path," -t ",treefile," -s ",alignment_path," --scf ",number_scf_quartets," -nt 1 -redo -safe")
    system(call)
  }
  ## Retrieve the site concordance factors from the output table
  scf_table <- read.table(paste0(alignment_path,".treefile.cf.stat"), header = TRUE, sep = "\t")
  scf_results <- list(mean_scf = round(mean(scf_table$sCF), digits = 2), 
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



#### Tree estimation functions ####
estimate.iqtree2.gene.trees <- function(gene_folder, iqtree2_path){
  
}



#### Utility functions ####
convert.to.nexus <- function(alignment_path, sequence_format = "DNA"){
  ## Convert fasta file to nexus file (if there is no existing nexus file with the same name)
  
  ## Prepare parameters for file conversion
  # Name nexus file by simply appending ".nex" to end of existing file name
  nexus_alignment_path <- paste0(alignment_path,".nex")
  # Extract file type from alignment path
  suffix <- tail(strsplit(alignment_path,"\\.")[[1]],1)
  # Set format for output nexus file
  if ((sequence_format == "DNA") | (sequence_format == "dna")){
    nexus_format = "dna"
  } else if ((sequence_format == "Protein") | (sequence_format == "protein") | 
             (sequence_format == "AA") | (sequence_format == "aa")){
    nexus_format = "protein"
  }
  
  ## Convert to nexus using funcions based on suffix
  if (suffix == "fasta" |suffix == "fa" | suffix == "fna" | suffix == "ffn" | suffix == "faa" | suffix == "frn" | suffix == "fas"){
    ## If the file is a fasta file, convert it to nexus file format (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE){
      # Read in the fasta data
      data <- read.FASTA(alignment_path, type = sequence_format)
      # Write out the nexus data
      write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, interleaved = FALSE, datablock = FALSE) # write the output as a nexus file)
    }
  } else if (suffix == "phy" | suffix == "phylip"){
    ## If the file is a phy file, convert it to nexus file format (unless a nexus version already exists)
    if (file.exists(nexus_alignment_path) == FALSE){
      data <- read.phy(alignment_path)
      write.nexus.data(data, file = nexus_alignment_path,format = nexus_format, interleaved = FALSE, datablock = FALSE) # write the output as a nexus file)
    }
  }
  
  ## Open the nexus file and delete the interleave = YES or INTERLEAVE = NO part so IQ-TREE can read it
  nexus <- readLines(nexus_alignment_path)
  ind <- grep("BEGIN CHARACTERS",nexus)+2
  if ((sequence_format == "DNA") | (sequence_format == "dna")){
    nexus[ind] <- "  FORMAT DATATYPE=DNA MISSING=? GAP=-;"
  } else if ((sequence_format == "Protein") | (sequence_format == "protein") | 
             (sequence_format == "AA") | (sequence_format == "aa")){
    nexus[ind] <- "  FORMAT MISSING=? GAP=- DATATYPE=PROTEIN;"
  }
  # Write the edited nexus file out 
  writeLines(nexus,nexus_alignment_path)
  
  ## Output file name and path for nexus file
  return(nexus_alignment_path)
}



genes.from.alignment <- function(alignment_path, partition_path, gene_folder, sequence_format = "DNA", save.gene.info = FALSE){
  ## Function to take an alignment and a partition file, and save separate alignments for each gene (and a data frame with info about the genes)
  
  # Open the alignment as a matrix
  if (sequence_format == "DNA" | sequence_format == "dna"){
    aln_mat <- as.matrix(read.FASTA(alignment_path, type = "DNA"))
  } else if (sequence_format == "AA" | sequence_format == "aa" | sequence_format == "protein" | sequence_format == "Protein"){
    aln_mat <- as.matrix(read.FASTA(alignment_path, type = "AA"))
  }
  
  # Open partition file and extract genes
  partition_file <- readLines(partition_path)
  charpartition_inds <- grep("charset", partition_file)
  charpartition_lines <- partition_file[charpartition_inds]
  
  # Save each gene individually
  gene_info_list <- lapply(1:length(charpartition_lines), save.one.gene, charpartitions = charpartition_lines, alignment_matrix = aln_mat, output_folder = gene_folder)
  
  # If desired, save information about the genes
  if (save.gene.info == TRUE){
    gene_info_df <- do.call(rbind.data.frame, gene_info_list)
    names(gene_info_df) <- c("gene_name", "gene_start_position", "gene_end_position", "gene_length", "sequence_format")
    
    write.csv(gene_info_df, file = paste0(dirname(alignment_path), "/genes_from_alignment_info.csv"))
  }
}

save.one.gene <- function(index, charpartitions, alignment_matrix, output_folder, sequence_format = "DNA"){
  ## Small function to process a single charpartition row and save the associated gene
  
  # Extract the i^th row using the index and split the line to determine the gene name and location
  c_line <- charpartitions[index]
  c_line <- gsub(";", "", c_line)
  gene_name_split <- unlist(strsplit(c_line, "="))
  gene_name <- gsub("\tcharset ", "", gene_name_split[1])
  gene_name <- gsub(" ", "", gene_name)
  gene_range <- unlist(strsplit(c_line, ","))[2]
  gene_range <- gsub(" ", "", gene_range)
  gene_range_split <- strsplit(gene_range, "-")[[1]]
  gene_start <- as.numeric(gene_range_split[[1]])
  gene_end <- as.numeric(gene_range_split[[2]])
  gene_length = length(gene_start:gene_end)
  
  # Subset the matrix to include only the sites for gene 1 for all taxa
  # matrix[1:x, 1:y] where x is the number of taxa and y is the number of sites
  gene_alignment <- alignment_matrix[, gene_start:gene_end]
  
  # Name and output gene file
  gene_filepath <- paste0(output_folder, gene_name, ".fa")
  write.FASTA(gene_alignment, file = gene_filepath)
  
  # Assemble a vector of information about the gene
  gene_info <- c(gene_name, gene_start, gene_end, gene_length, sequence_format)
  names(gene_info) <- c("gene_name", "gene_start_position", "gene_end_position", "gene_length", "sequence_format")
  return(gene_info)
}



