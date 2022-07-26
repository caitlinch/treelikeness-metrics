# /caitlinch/treelikeness-metrics/00_prepare_empirical_data.R
# Caitlin Cherryh 2022

# This program takes in a partition file and separates an alignment based on the charsets in the partition

# Paper:
#   Oaks, J.R. (2011), A time-calibrated species tree of Crocodylia reveals a recent radiation of the true crocodiles. Evolution, 65: 3285-3297. https://doi.org/10.1111/j.1558-5646.2011.01373.x

# Data available from: 
#   Oaks, Jamie R (2011), Data from: A time-calibrated species tree of Crocodylia reveals a recent radiation of the true crocodiles, Dryad, Dataset, https://doi.org/10.5061/dryad.5k9s0

#### 1. Set parameters ####
# alignment_file          <- nexus alignment file
# partition_files         <- nexus partition files for nuclear and mtDNA
# output_directory        <- Directory for output alignments
# repo_directory          <- Location of caitlinch/treelikeness-metrics github repository (for access to functions)

alignment_file <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_alignment.nex"
partition_files <- c("/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_mtDNA.nex",
                     "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/oaks2011_crocodylia_nDNA.nex")
alignment_directory <- "/Users/caitlincherryh/Documents/C2_TreelikenessMetrics/00_data_oaks2011/"
repo_directory <- "/Users/caitlincherryh/Documents/Repositories/treelikeness-metrics/"



#### 2. Prepare analyses ####
# Open libraries
library(phylotools)
library(phangorn)

# Convert nexus to fasta
fasta_file <- gsub("\\.nex", "\\.fa", alignment_file)
# Open alignment file as matrix
dna_nex <- read.nexus.data(alignment_file)
write.dna(dna_nex, file = fasta_file, format = "fasta", colsep = "")
# Open fasta file
dna_mat <- as.matrix(read.dna(file = fasta_file, format = "fasta"))



#### 4. Save each gene individually ####
charset_type = c("mtDNA", "nDNA")
for (c in charset_type){
  # Make new output directory
  output_directory <- paste0(alignment_directory, c, "_gene_alignments/")
  if (dir.exists(output_directory) == FALSE){dir.create(output_directory)}
  
  # Select partition file
  c_partition_file <- grep(c, partition_files, value = TRUE)
  # Open partition file as lines
  partitions <- readLines(c_partition_file)
  # Extract only lines that contain a charset and remove extra formatting
  partitions <- grep("charset", partitions, value = TRUE)
  partitions <- gsub("\tcharset ", "", partitions)
  partitions <- gsub(";", "", partitions)
  
  ## Save each gene
  print('Saving individual genes:')
  all_gene_lengths <- c()
  all_gene_names <- c()
  for (i in 1:length(partitions)){
    # Get the ith gene from the list of genes
    gene_line <- partitions[i]
    
    # Split the line to determine the gene name
    gene_line <- gsub(";", "", gene_line)
    gene_line_split <- strsplit(gene_line, "=")[[1]]
    gene_name <- gsub(" ", "", gene_line_split[1])
    print(gene_name)
    
    # Check whether this line is for a single codon position or a whole gene
    get_codon_position <- grepl("\\\\3", gene_line)
    
    # Split the line to determine the gene start and end position
    gene_range <- gsub(" ", "", gene_line_split[2])
    gene_range_split <- strsplit(gene_range, "-")[[1]]
    if (length(gene_range_split) == 2){
      gene_start <- as.numeric(gene_range_split[[1]])
      if (get_codon_position == FALSE){
        # If getting whole gene, just take end of range
        gene_end <- as.numeric(gene_range_split[[2]])
        # Get the whole range of numbers as the gene_positions to extract
        gene_positions <- gene_start:gene_end
      } else if (get_codon_position == TRUE){
        # If getting only one codon position, get only last number (do not include the "\\3")
        gene_end <- as.numeric(gsub("\\\\3", "", gene_range_split[[2]]))
        # Get every third number starting from the gene_start and moving towards the gene_end
        gene_positions <- seq(gene_start, gene_end, 3)
      }
    } else if (length(gene_range_split) > 2){
      # If there are multiple components to the partition, extract all the sites for inclusion
      all_ranges <- strsplit(gene_line_split[2], " ")[[1]]
      all_ranges <- all_ranges[which(all_ranges != "")]
      gene_positions <- c()
      for (ind in 1:length(all_ranges)){
        j <- all_ranges[ind]
        j_get_codon_position <- grepl("\\\\3", j)
        if (j_get_codon_position == FALSE){
          j_split <- strsplit(j, "-")[[1]]
          j_start <- as.numeric(j_split[[1]])
          j_end <- as.numeric(j_split[[2]])
          j_positions <- c(j_start:j_end)
          gene_positions <- c(gene_positions, j_positions)
        } else if  (j_get_codon_position == TRUE){
          j_2 <- gsub("\\\\3", "", j)
          j_split <- strsplit(j_2, "-")[[1]]
          j_start <- as.numeric(j_split[[1]])
          j_end <- as.numeric(j_split[[2]])
          j_positions <- seq(j_start, j_end, 3)
          gene_positions <- c(gene_positions, j_positions)
        }
      }
    }
    
    # Length of gene will just be the number of sites included in the gene_positions
    gene_length <- length(gene_positions)
    
    # Save name and length of gene
    all_gene_lengths <- c(all_gene_lengths, gene_length)
    all_gene_names <- c(all_gene_names, gene_name)
    
    # Sort gene positions into numerical order
    # Subset the supermatrix to get the sites (columns) for this gene
    sorted_gene_positions <- sort(gene_positions, decreasing = FALSE)
    gene_mat <- dna_mat[, c(sorted_gene_positions)]
    
    # Assemble file name for gene
    gene_file <- paste0(output_directory, gene_name, ".fa")
    # Write gene to file
    write.FASTA(gene_mat, file = gene_file)
  }
  
  # Create dataframe with gene name and lengths
  df <- data.frame(dataset = "Oaks2011", gene_name = all_gene_names, gene_length = all_gene_lengths)
  write.csv(df, file = paste0(alignment_directory, "Oaks2011_", c,"_gene_lengths.csv"), row.names = FALSE)
}
